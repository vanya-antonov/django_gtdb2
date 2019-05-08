# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint
import logging
import os
import re

from Bio.Data import CodonTable
from django.conf import settings

from chelatase_db.lib.bio import run_tblastn_seq_vs_db, run_tblastn_file_vs_db
from chelatase_db.lib.config import PATHWAY_GENES_FAA
from chelatase_db.models.cof import ChelataseCof
from chelatase_db.models.feat import ChelataseFeat
from chelatase_db.models.fshift import ChelataseFshift
from chelatase_db.models.seq import ChelataseSeq
from gtdb2.models.org import Org


class ChelataseOrg(Org):
    class Meta:
        proxy = True

    # Modify values of some constant attributes
    FEAT_CLS_NAME = 'ChelataseFeat'
    FSHIFT_CLS_NAME = 'ChelataseFshift'

    # Merge two dicts: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(Org.PRM_INFO.items()) + list({
        'chel_num_chld': {'value_attr': 'num', 'type_fun': int},
    }.items()))

    @property
    def chld_feat_set(self):
        """Returns a QuerySet of ChelataseFeat objects corresponding to
        chlD or bchD genes.
        """
        # Follow the relationships by using the double underscores operator:
        # https://docs.djangoproject.com/en/2.1/topics/db/queries/#lookups-that-span-relationships
        return ChelataseFeat.objects.filter(
            seq__org=self,
            param_set__name='chel_gene_group',
            param_set__value='chlD_bchD')

    @property
    def chel_subunit_feat_set(self):
        """Returns a QuerySet of ChelataseFeat objects corresponding to
        the large, medium or small chelatase subunits (i.e. all the 
        features with the 'chel_subunit' prm.
        """
        return ChelataseFeat.objects.filter(
            seq__org=self,
            param_set__name='chel_subunit')

    def create_all_params(self):
        """Extend the parent's method to create additional params.
        """
        super().create_all_params()
        self._make_params_chel_statistics()

    def create_annotation(self):
        """Creates feats homologous to the cholorophyll and B12
        biosynthesis pathway genes.
        """
        user = self.gtdb.get_or_create_annotation_user()

        # Annotation is only created for orgs with chld gene(s)
        chld_feats = self._get_or_create_chld_feats(user)
        if len(chld_feats) == 0:
            return
        super().create_annotation()

        # Try to find cholorophyll and B12 pathway genes and create the feats
        db_path = self.gtdb.get_full_path_to(self.prm.blastdb_nucl_all)
        faa_fn = os.path.join(settings.BASE_DIR, PATHWAY_GENES_FAA)
        all_hsps = run_tblastn_file_vs_db(faa_fn, db_path, self.transl_table)
        for hsp in all_hsps:
            seq = ChelataseSeq.get_or_create_from_ext_id(
                user, self, hsp.sbjct_id)
            ChelataseFeat.get_or_create_from_gbk_annotation(
                user, seq, hsp.sbjct_start, hsp.sbjct_end, hsp.sbjct_strand)

    def get_or_create_feats_from_cof_by_tblastn(self, user, cof):
        """Returns a list of feats (with or without fshifts) identified
        in the given org sequences by using prot_seqs from the seed COF
        as queries to run tblastn.
        """
        all_feats = set()
        for q_fshift in cof.seed_fshifts:
            # Search all the org nt seqs for regions whose translation is
            # similar to the N- and C- parts of the query fs-prot
            blastdb_path = self.gtdb.get_full_path_to(
                self.prm['blastdb_nucl_all'])

            n_hsp_list = run_tblastn_seq_vs_db(
                q_fshift.prm['seq_prot_n'], blastdb_path, self.transl_table)
            c_hsp_list = run_tblastn_seq_vs_db(
                q_fshift.prm['seq_prot_c'], blastdb_path, self.transl_table)

            n_hits_dict = _make_hits_dict(n_hsp_list)
            c_hits_dict = _make_hits_dict(c_hsp_list)

            feat_set = self._get_or_create_feats_from_tblastn_hits(
                user, cof, n_hits_dict, c_hits_dict)
            all_feats.update(feat_set)

        return list(all_feats)

    def _make_params_chel_statistics(self):
        """Computes chelatase-related statistics and saves it as params:
         - chel_num_chld - total number of chlD/bchD genes
         - chel_num_M - total number of M chelatase subunits (cobT/chlD/bchD).
        """
        self.set_param('chel_num_chld', num=self.chld_feat_set.count())

    def _get_or_create_chld_feats(self, user):
        """Returns a list of ChelataseFeat objects."""
        chld_cof = ChelataseCof.get_or_create_chld_cof(self.user)
        chld_feats = self.get_or_create_feats_from_cof_by_tblastn(
            user, chld_cof)

        # Verify the automatic annotation of chlD genes
        for chld in chld_feats:
            if 'chel_subunit' not in chld.prm or chld.prm.chel_subunit != 'M':
                self.delete()
                raise ValueError("The identified chlD gene was not "
                                 "automatically annotated as medium subunit!")
        return chld_feats

    def _get_or_create_feats_from_tblastn_hits(
            self, user, cof, n_hits_dict, c_hits_dict):
        """Creates feats and fshifts (if needed) based on the tblastn hits.
        Returns a set of Feat objects.
        """
        feat_set = set()
        for chr_name in sorted(n_hits_dict.keys()):
            if chr_name not in c_hits_dict:
                continue

            seq = ChelataseSeq.get_or_create_from_ext_id(
                user, self, chr_name)
            all_fs_info = self._get_fshift_info_from_tblastn_hits_on_seq(
                seq, n_hits_dict[chr_name], c_hits_dict[chr_name])

            for fs in all_fs_info:
                if fs['coord'] is None:
                    # This is normal CDS -- no need to create fshift
                    feat = ChelataseFeat.get_or_create_from_gbk_annotation(
                        user, seq, fs['start'], fs['end'], fs['strand'])
                else:
                    # fsCDS needs a frameshift for full length translation
                    fshift = ChelataseFshift.get_or_create(
                        user=user, seq=seq, origin='tblastn',
                        start=fs['start'], end=fs['end'], strand=fs['strand'],
                        coord=fs['coord'], len=fs['len'])

                    # Create the parent feature that corresponds to the
                    # upstream part of fsCDS, i.e. where translation begins
                    parent = ChelataseFeat.get_or_create_from_gbk_annotation(
                        user, seq, fs['hsp_n'].sbjct_start,
                        fs['hsp_n'].sbjct_end, fs['hsp_n'].sbjct_strand)

                    # Finally, create the full-length fsCDS feat
                    feat = ChelataseFeat.get_or_create_fscds_from_parent(
                        user, parent, {fshift})

                if feat is not None:
                    feat_set.add(feat)

        return feat_set

    def _get_fshift_info_from_tblastn_hits_on_seq(
            self, seq, all_hsp_n, all_hsp_c):
        """Processes all the tblastn hits on the same seq and returns a
        list of dicts with info about fshifts to be created.
        """
        all_fs_info = []
        while len(all_hsp_n) > 0 and len(all_hsp_c) > 0:
            hsp_n, hsp_c = _get_closest_hsp_pair(all_hsp_n, all_hsp_c)
            if hsp_n is None or hsp_c is None:
                break

            # Each hsp can only be used for 1 fshift
            all_hsp_n.remove(hsp_n)
            all_hsp_c.remove(hsp_c)

            fs_dict = _get_fshift_info_from_two_hits(
                hsp_n, hsp_c, seq.seq, self.transl_table)
            if fs_dict is not None:
                fs_dict['hsp_n'] = hsp_n
                fs_dict['hsp_c'] = hsp_c
                all_fs_info.append(fs_dict)

        return all_fs_info

def _get_closest_hsp_pair(all_hsp_n, all_hsp_c):
    min_dist, min_hsp_n, min_hsp_c = None, None, None
    for hsp_n in all_hsp_n:
        for hsp_c in all_hsp_c:
            if hsp_n.sbjct_strand == hsp_c.sbjct_strand == 1:
                #          hsp_n                    hsp_c
                #    ----------------       --------------------
                #    start        end       start            end
                #   ==============================================>

                # the hits should be in the right order
                if hsp_n.sbjct_end >= hsp_c.sbjct_end:
                    continue
                dist = hsp_c.sbjct_start - hsp_n.sbjct_end
            elif hsp_n.sbjct_strand == hsp_c.sbjct_strand == -1:
                #          hsp_c                    hsp_n
                #    ----------------       --------------------
                #    start        end       start            end
                #   <=============================================

                # the hits should be in the right order
                if hsp_c.sbjct_end >= hsp_n.sbjct_end:
                    continue
                dist = hsp_n.sbjct_start - hsp_c.sbjct_end
            else:
                # The hits are on different strands
                continue

            if dist < 0:
                logging.warning("Distance is negative for hsp pair:\n"
                                "'%s'\n'%s'" % (vars(hsp_n), vars(hsp_c)))
                dist = 0

            if min_dist is None or dist < min_dist:
                min_dist, min_hsp_n, min_hsp_c = dist, hsp_n, hsp_c

    return min_hsp_n, min_hsp_c

def _get_fshift_info_from_two_hits(hsp_n, hsp_c, chr_seq, gencode):
    if hsp_n.sbjct_strand == hsp_c.sbjct_strand == 1:
        lh, rh = hsp_n, hsp_c
    elif hsp_n.sbjct_strand == hsp_c.sbjct_strand == -1:
        lh, rh = hsp_c, hsp_n
    else:
        raise ValueError("The HSPs are on different strands!")

    l_ss = _get_stop_stop_seq(lh.sbjct_start, lh.sbjct_end, lh.sbjct_strand,
                              chr_seq, gencode)
    r_ss = _get_stop_stop_seq(rh.sbjct_start, rh.sbjct_end, rh.sbjct_strand,
                              chr_seq, gencode)
    if l_ss is None or r_ss is None:
        return None
    elif l_ss['left'] == r_ss['left'] and l_ss['right'] == r_ss['right']:
        # left and right hits are in frame => fs-gene without FS
        return {'coord': None, 'len': 0,   # i.e. no frameshift
                'start': lh.sbjct_start, 'end': rh.sbjct_end,
                'strand': lh.sbjct_strand}

    ovlp_l, ovlp_r = get_overlap_region(l_ss['left'], l_ss['right'],
                                        r_ss['left'], r_ss['right'])
    if ovlp_l is None or ovlp_r is None:
        # No overlap  =>  No fs-gene
        return None
    elif (ovlp_r - ovlp_l) <= 0:
        logging.error('Something is wrong with the overlap coords (%i,%i)' %
                      (ovlp_l,ovlp_r))
        return None
    else:
        # True fs-gene
        gene_len = rh.sbjct_end - lh.sbjct_start
        if gene_len % 3 == 1:
            fs_len = 1   # i.e. +1
        elif gene_len % 3 == 2:
            fs_len = -1
        else:
            fs_len = 0   # i.e. in_frame_stop

        # putative fs-coord just before the middle stop codon
        fs_coord = (l_ss['right']-3) if lh.sbjct_strand == 1 else (r_ss['left']+3)
        return {'coord': fs_coord, 'len': fs_len, 'strand': lh.sbjct_strand,
                'start': lh.sbjct_start, 'end': rh.sbjct_end}

def get_overlap_region(s1,e1,s2,e2):
    """0-based system is used (like in Biopython):

    |            INPUT            |          RETURNS            |
    |-----------------------------|-----------------------------|
    |                             |                             |
    |  s1=3      e1=14            |                             |
    |     |----------|            |        [9,14]               |
    |           |---------|       |                             |
    |        s2=9     e2=19       |                             |
    |                             |                             |
    """
    if s1 > e1 or s2 > e2:
        raise Exception("Something is wrong with the intervals (%i,%i) and (%i,%i)" %
                        (s1,e1,s2,e2))

    if s1 <= s2 <= e1 and s1 <= e2 <= e1:
        # |----------------|
        #     |--------|
        return s2, e2
    elif s2 <= s1 <= e2 and s2 <= e1 <= e2:
        #     |--------|
        # |----------------|
        return s1, e1
    elif s1 <= s2 <= e1 and s2 <= e1 <= e2:
        # |------------|
        #       |-------------|
        return s2, e1
    elif s2 <= s1 <= e2 and  s1 <= e2 <= e1:
        #       |-------------|
        # |------------|
        return s1, e2
    else:
        return None, None

def _get_stop_stop_seq(left, right, strand, chr_seq, gencode):
    """Expand the given CDS part in both directions until stop codons
    on both sides (the stop codons are included as well).
    """
    # Validate the provided coordinates
    if (right - left) % 3 != 0:
        raise ValueError("The difference between the left and right "
                         "coordinates is not divisible by 3!!")

    chr_region = chr_seq[left:right]
    if not re.compile('^[ACGT]+$', re.IGNORECASE).match(str(chr_region)):
        logging.warning("Sequence %i-%i contains non-ACGT chars" %
                        (left, right))
        return None

    if strand == -1:
        chr_region = chr_region.reverse_complement()
    prot_seq = chr_region.translate(table=gencode)
    if '*' in prot_seq:
        logging.warning('Given region (%i,%i) already contains a stop codon' %
                        (left, right))
        return None

    ss_left = left
    while ss_left >= 3:
        cur_codon = chr_seq[ (ss_left-3) : ss_left ]
        codon_type = get_codon_type(cur_codon, gencode, strand)
        if codon_type == 'coding':
            ss_left -= 3  # expand the stop-stop seq
        elif codon_type == 'stop':
            ss_left -= 3  # stop-stop seq includes stops as well
            break
        else:
            break

    ss_right = right
    while ss_right < len(chr_seq)-3:
        cur_codon = chr_seq[ ss_right : (ss_right+3) ]
        codon_type = get_codon_type(cur_codon, gencode, strand)
        if codon_type == 'coding':
            ss_right += 3  # expand the stop-stop seq
        elif codon_type == 'stop':
            ss_right += 3  # stop-stop seq includes stops as well
            break
        else:
            break

    return {'left' : ss_left, 'right' : ss_right}

def get_codon_type(codon, gencode, strand=1):
    acgt_only_re = re.compile('^[ACGT]+$', re.IGNORECASE)
    if not acgt_only_re.match(str(codon)):
        return 'non_ACGT'

    if strand == -1:
        codon = codon.reverse_complement()

    if codon in CodonTable.unambiguous_dna_by_id[gencode].stop_codons:
        return 'stop'
    else:
        return 'coding'

def _make_hits_dict(hsp_list):
    """Returns a dict where the keys are the hit sequence IDs and
    the values are the lists of Bio.Blast.Record.HSP objects.
    """
    all_hits = {}
    for hsp in hsp_list:
        all_hits.setdefault(hsp.sbjct_id, []).append(hsp)
    return all_hits

