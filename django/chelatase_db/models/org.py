# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint
import logging
import os
import re

from django.conf import settings

from chelatase_db.lib.bio import run_tblastn_seq_vs_db, run_tblastn_file_vs_db, get_stop_stop_seq
from chelatase_db.lib.config import (
    PATHWAY_GENES_FAA, read_pathway_gene_info, read_kegg_orgs)
from chelatase_db.models.cof import ChelataseCof
from chelatase_db.models.feat import ChelataseFeat
from chelatase_db.models.fshift import ChelataseFshift
from chelatase_db.models.seq import ChelataseSeq
from gtdb2.lib.bio import get_overlapping_feats_from_list
from gtdb2.models.org import Org

#
# from chelatase_db.models import ChelataseFeat
# from dna_features_viewer import GraphicFeature, GraphicRecord
# from gtdb2.models import Org
# from gtdb2.models import Seq
# from gtdb2.models import Fshift
from io import StringIO



class ChelataseOrg(Org):
    class Meta:
        proxy = True

    # Modify values of some constant attributes
    FEAT_CLS_NAME = 'ChelataseFeat'
    FSHIFT_CLS_NAME = 'ChelataseFshift'

    # Merge two dicts: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(Org.PRM_INFO.items()) + list({
        'chel_num_chld': {'value_attr': 'num', 'type_fun': int},
        'chel_genotype_genes': {},
        'kegg_org_code': {},
        'chel_synthesis_chl': {'value_attr': 'num', 'type_fun': bool},
        'chel_synthesis_b12': {'value_attr' : 'num', 'type_fun' : bool},
        }.items()))

    # Image location: https://github.com/vanya-antonov/django_gtdb2/tree/master/frontend/src/assets/images
    CHLD_FUNCTIONS = {
        'M_Mg': {
            'image': 'chelTable0.png',
            'title': 'Medium subunit of Mg-chelatase',
            'descr': 'chlD gene encodes Medium subunit of Mg-chelatase'
        },
        'M_Mg_but_not_Co': {
            'image': 'chelTable1.png',
            'title': 'Medium subunit of Mg-chelatase',
            'descr': 'chlD gene encodes Medium subunit of Mg-chelatase while ' +
                     'Co-chelatase is encoded by a separate set of genes'
        },
        'M_Mg_and_M_Co': {
            'image': 'chelTable2.png',
            'title': 'Medium subunit of Mg- and Co-chelatase',
            'descr': 'chlD gene encodes both the Medium subunit of Mg-chelatase ' +
                     'and the Medium subunit of Co-chelatase'
        },
        'M_Co': {
            'image': 'chelTable3.png',
            'title': 'Medium subunit of Co-chelatase',
            'descr': 'chlD gene encodes Medium subunit of Co-chelatase'
        },
        'M_and_S_Co_fs': {
            'image': 'chelTable4.png',
            'title': 'Medium and Small subunits of Co-chelatase',
            'descr': 'Frameshifted chlD gene encodes both Medium and Small ' +
                     'subunits of Co-chelatase (without frameshifting Small ' +
                     'subunit is produced while frameshifting allows synthesis ' +
                     'of the full-length Medium subunit)'
        },
        'M_and_S_Co': {
            'image': 'chelTable5.png',
            'title': 'Medium and Small subunits of Co-chelatase',
            'descr': 'chlD gene encodes both Medium and Small ' +
                     'subunits of Co-chelatase (without frameshifting Large ' +
                     'subunit is synthesized while frameshifting changes the ' +
                     'reading frame and Small subunit is produced due to out-of-frame ' +
                     'prepature stop codon)'
        }
    }

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

    def get_full_pathway_gene_dict(self, pathway=None):
        """Returns a dict of lists of feats (genes) that have the
        'chel_pathway' and 'chel_gene_group' params. They keys are
        chel_gene_group values (e.g. 'chlI_bchI' or 'cobN').

        Arguments:
         - pathway - (string) include feats from a particular pathway
           only ('B12' or 'Chlorophyll' -- see pathway_genes.txt). By default
           it returns genes from all the pathways.
        """

        # Get statistics about all the B12/CHL genes
        info_dict = read_pathway_gene_info()
        all_gene_groups = {}
        for gene in info_dict.values():
            if pathway is not None and gene.get('chel_pathway', '') != pathway:
                continue

            # e.g. 'cobD_cobC'
            gene_group = gene['chel_gene_group']
            all_gene_groups[gene_group] = []

        # Add B12/CHL genes from the current org
        for feat in self.feat_set.all():
            if 'chel_pathway' not in feat.prm:
                continue
            if 'chel_gene_group' not in feat.prm:
                continue

            gene_group = feat.prm['chel_gene_group']
            if gene_group in all_gene_groups:
                all_gene_groups[gene_group].append(feat)

        return all_gene_groups

    def create_all_params(self):
        """Extend the parent's method to create additional params.
        """
        super().create_all_params()
        self._make_params_chel_statistics()
        self._make_genotype_str()
        self._make_params_kegg()
        self._make_params_synthesize_classification() # TODO проверить на корректность

    def create_annotation(self):
        """Creates feats homologous to the cholorophyll and B12
        biosynthesis pathway genes.
        """
        user = self.gtdb.get_or_create_annotation_user()

        # Annotation is only created for orgs with chld gene(s)
        logging.info("Searching for chlD genes...")
        chld_feats = self._get_or_create_chld_feats(user)
        logging.info("Number of identified chlD genes = %s" % len(chld_feats))
        if len(chld_feats) == 0:
            return
        super().create_annotation()

        # Try to find cholorophyll and B12 pathway genes and create the feats
        db_path = self.gtdb.get_full_path_to(self.prm.blastdb_nucl_all)
        faa_fn = os.path.join(settings.BASE_DIR, PATHWAY_GENES_FAA)
        all_hsps = run_tblastn_file_vs_db(faa_fn, db_path, self.transl_table)
        for hsp in all_hsps:
            logging.info(
                "Processing hit '%s' from query '%s' (evalue=%.1e)" %
                (hsp.sbjct_id, hsp.query_id, hsp.expect))
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

    def _make_params_synthesize_classification(self):
        """Computes biochemistry path :
         - chel_synthesis_chl - bool type if org synthesizes Chlorophyll
         - chel_synthesis_b12 - bool type if org synthesizes b12.
        """
        Chlorophyll_bool, B12_bool = self.get_prediction_text(output2text=False, from_db=False)

        self.set_param('chel_synthesis_chl', num=Chlorophyll_bool)
        self.set_param('chel_synthesis_b12', num=B12_bool)

    def _make_genotype_str(self):
        """Creates the 'chel_genotype_genes' param like '1xcobN, 1xfs-chlD'.
        """
        group2name = {
            'cobN': 'cobN',
            'cobT': 'cobT',
            'cobS': 'cobS',
            'chlH_bchH': 'chlH',
            'chlD_bchD': 'chlD',
            'chlI_bchI': 'chlI'}

        all_pathway_genes = self.get_full_pathway_gene_dict()
        chel_genes = {}
        for key, feats in all_pathway_genes.items():
            if key == 'chlD_bchD':
                # Separate the fs-chlD and chlD
                chel_genes['chlD'] = [f for f in feats if f.fshift is None]
                chel_genes['fs-chlD'] = [f for f in feats if f.fshift is not None]
            elif key == 'chlI_bchI':
                # Remove chlI that are the short products of fs-chlD genes
                chel_genes['chlI'] = [
                    f for f in feats if f.children_feat_set.count() == 0]
            elif key in group2name:
                # Change some names: 'chlD_bchD'  =>   'chlD'
                new_key = group2name[key]
                chel_genes[new_key] = feats

        genotype_names = [
            'cobN', 'cobT', 'cobS', 'chlH', 'chlD', 'fs-chlD', 'chlI']
        genotype_parts = []
        for name in genotype_names:
            num = len(chel_genes[name])
            if num > 0:
                genotype_parts.append(str(num) + 'x' + name)
        genotype_str = ', '.join(genotype_parts)

        self.set_param('chel_genotype_genes', genotype_str)

    def _make_params_kegg(self):
        """Creates params that will allow to show links to the KEGG
        database like: https://www.genome.jp/kegg-bin/show_pathway?dac00860
        """
        kegg_dict = read_kegg_orgs()
        if self.name in kegg_dict:
            self.set_param('kegg_org_code', kegg_dict[self.name]['org_code'])
            return

        # Couldn't find the exact name in the KEGG table -- try to use species
        my_species = _org_name_to_species(self.name)
        if my_species is None:
            return

        for kegg_org, kegg_info in kegg_dict.items():
            if my_species == _org_name_to_species(kegg_org):
                self.set_param('kegg_org_code', kegg_info['org_code'])
                return

    def _get_or_create_chld_feats(self, user):
        """Returns a list of ChelataseFeat objects."""
        chld_cof = ChelataseCof.get_or_create_chld_cof(self.user)
        chld_feats = self.get_or_create_feats_from_cof_by_tblastn(
            user, chld_cof)

        # Verify the automatic annotation of chlD genes
        for chld in chld_feats:
            if 'chel_subunit' not in chld.prm or chld.prm.chel_subunit != 'M':
                logging.warning(
                    "The identified chlD gene (%s) was not automatically "
                    "annotated as medium subunit and will be removed" %
                    chld.prm.location_str)
                if chld.fshift is not None:
                    chld.fshift.delete()
                chld.delete()
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
                elif fs['len'] is None or fs['len'] == 0:
                    # TODO: there is something strange with Rhodococcus
                    # kunmingensis (fs_len == 0  -- is it in-frame stop?)
                    logging.warning(
                        "Strange fs-gene (fs_coord = %s, but fs_len=%s) -- "
                        "skipping..." % (fs['coord'], fs['len']))
                    continue
                else:
                    # fsCDS needs a frameshift for full length translation
                    fshift = ChelataseFshift.get_or_create(
                        user=user, seq=seq, origin='tblastn',
                        start=fs['start'], end=fs['end'], strand=fs['strand'],
                        coord=fs['coord'], len=fs['len'])

                    # Create the parent feature that corresponds to the
                    # upstream part of fsCDS, i.e. where translation begins
                    parent = _get_or_create_parent_feat_from_tblastn_hits(
                        user, seq, fs['hsp_n'], fs['hsp_c'])

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

    @property
    def plot_gen_diagram(self):

        from chelatase_db.models import ChelataseFeat
        from dna_features_viewer import GraphicFeature, GraphicRecord
        from gtdb2.models import Org
        from gtdb2.models import Seq
        from gtdb2.models import Fshift

        id_org = self.id

        def calculate_label(x):
            if x // 10 ** 6 >= 1:
                return str(round(x / 10 ** 6, 2)) + 'M'
            else:
                return str(round(x / 10 ** 3)) + 'K'

        def get_color(gene):
            if 'fsh' in gene: return "#f60a01"

         #   if 'chl' in gene:
            if str(gene).startswith('bch') or str(gene).startswith('chl'):
                return '#93d14f'  # '#ffcccc'
            
            else:
                return "#8fb5e5"

        c_org = Org.objects.get(param_set__name='chel_genotype_genes', id=id_org)
        f_org = ChelataseFeat.objects.filter(seq__org=c_org).all()

        ch_set = ChelataseFeat.objects.filter(seq__org=id_org, param_set__name='chel_subunit').all()
        all_parents = [i.parent_id for i in ch_set if i.parent_id is not None]

        f_org = [i for i in ch_set if i.id not in all_parents]  # Извлекаем все гены
        fsh_set = ChelataseFshift.objects.filter(seq__org=c_org).all()  # Извлекаем все фрэйм-шифты

        dict_plot, locs, fshs = {'max_scale': 10, 'gen_len_max': 0, 'min_scale': 10 ** 8, 'min_loc': 10 ** 10}, [0], []

        for i in f_org:
            if max(i.end, i.start) > dict_plot['max_scale']: dict_plot['max_scale'] = max(i.end, i.start)
            if i.end - i.start < dict_plot['min_scale']: dict_plot['min_scale'] = i.end - i.start
            if i.end - i.start > dict_plot['gen_len_max']: dict_plot['gen_len_max'] = i.end - i.start
            if min(i.end, i.start) < dict_plot['min_loc']: dict_plot['min_loc'] = min(i.end, i.start)

        scale = 100  # Будет считать, что вся шкала это 100 единиц

        dict_plot['max_scale'] /= scale
        dict_plot['min_loc'] /= dict_plot['max_scale']
        dict_plot['gen_len_max'] /= dict_plot['max_scale']

        for i in f_org:  # Нормируем положения в генах
            loc = ((i.start + i.end) / 2) / dict_plot['max_scale']
            locs.append(loc)
            fshs.append(loc)
            dict_plot[i.name] = [{'len': (i.end - i.start) / dict_plot['max_scale'],
                                  'start': i.start / dict_plot['max_scale'], 'end': i.end / dict_plot['max_scale'],
                                  'loc': loc, 'label': i.prm.chel_gene, 'strand': i.strand, 'type': "gene"}]

        for i in fsh_set:  # нормируем положения в фрэймшифтах
            loc = ((i.start + i.end) / 2) / dict_plot['max_scale']
            locs.append(loc)
            dict_plot[i.name] = [{'len': (i.end - i.start) / dict_plot['max_scale'],
                                  'start': i.start / dict_plot['max_scale'], 'end': i.end / dict_plot['max_scale'],
                                  'loc': loc, 'type': "fsh"}]
        locs = sorted(locs)
        new_locs, tabs = {0: 0}, {}  # in tabs there are loc and distance
        left, right, shift_tabs, shift_genes = 0, 0, 0, 0

        min_gen_len = 10.10 * dict_plot['gen_len_max']
        min_dist = 0.20 * 100
        for i in locs:
            right = i
            if right - shift_tabs - left > min_dist and left != 0 and left != max(locs):
                tabs[(right + left) / 2] = calculate_label((right - left) * dict_plot['max_scale'])
                #last_tabs = (right + left) / 2 - shift_tabs - (right - left) / 4

            if right - left < min_gen_len and left != 0 and abs(right - left) != 0 and right in fshs and right != max(
                    locs):
                shift_genes = 6
                new_locs[right] = right - shift_tabs + shift_genes
                left = right - shift_tabs + shift_genes
            else:
                new_locs[right] = right - shift_tabs
                left = right

        features = []
        len_shift = 1 / 8  # Размер фреймшифта
        delta = 9  # Отступ от края
        for k, v in dict_plot.items():
            if k in ['max_scale', 'min_scale', 'gen_len_max', 'min_loc']: continue

            if dict_plot[k][0]['type'] == 'gene' and "_fs_" in k:
                features.append(GraphicFeature(start=new_locs[dict_plot[k][0]['loc']] - 4,
                                               end=new_locs[dict_plot[k][0]['loc']] - len_shift,
                                               strand=dict_plot[k][0]['strand'],
                                               color=get_color(dict_plot[k][0]['label']),
                                               label=dict_plot[k][0]['label']))
                features.append(GraphicFeature(start=new_locs[dict_plot[k][0]['loc']] + len_shift,
                                               end=new_locs[dict_plot[k][0]['loc']] + 2,
                                               strand=dict_plot[k][0]['strand'],
                                               color=get_color(dict_plot[k][0]['label']),
                                               label=dict_plot[k][0]['label']))

                features.append(GraphicFeature(start=new_locs[dict_plot[k][0]['loc']],  # Табличка
                                               end=new_locs[dict_plot[k][0]['loc']],
                                               strand=dict_plot[k][0]['strand'],
                                               color="#ccccff", label=k))

            elif dict_plot[k][0]['type'] == 'gene':
                features.append(
                    GraphicFeature(start=new_locs[dict_plot[k][0]['loc']], end=new_locs[dict_plot[k][0]['loc']],
                                   strand=dict_plot[k][0]['strand'], color=get_color(dict_plot[k][0]['label']),
                                   label=k))  # Табличка

                features.append(
                    GraphicFeature(start=new_locs[dict_plot[k][0]['loc']] - 4, end=new_locs[dict_plot[k][0]['loc']] + 2,
                                   strand=dict_plot[k][0]['strand'], color=get_color(dict_plot[k][0]['label']),
                                   label=dict_plot[k][0]['label']))
            if dict_plot[k][0]['type'] == 'fsh':
                features.append(GraphicFeature(start=new_locs[dict_plot[k][0]['loc']] - len_shift,
                                               end=new_locs[dict_plot[k][0]['loc']] + len_shift,
                                               color=get_color('fsh'),
                                               thickness=14, label='fsh'))

        for lo, dist in tabs.items():
            features.append(GraphicFeature(start=lo - 1 / 5, end=lo + 1 / 5, strand=0, color="#ffd700",
                                           label="about {} \n nucleotides".format(dist)))

        record = GraphicRecord(first_index=dict_plot['min_loc'] - delta,
                               sequence_length=115 - dict_plot['min_loc'] + delta - shift_tabs,
                               labels_spacing=-9, features=features, feature_level_height=2)
        record.plot(figure_width=15, figure_height=3)
        ax, plo = record.plot(figure_width=15, figure_height=3);
        svg_pic = StringIO()
        ax.figure.savefig(svg_pic, format="svg", bbox_inches='tight')

        return svg_pic.getvalue()

    def get_prediction_text(self, output2text=True, from_db=True):

        """
        return:
            True/False - org ability to synthesize Chlorophyll (feature_12)
            True/False - org ability to synthesize B12  (feature_34)
        """

        id_org = self.id
        # Magnesium = ['chlH_bchH', 'chlD_bchD', 'chlI_bchI']
        Chlorophyll = ['bchE', 'chlB_bchB', 'chlG_bchG', 'chlL_bchL', 'chlM_bchM', 'chlN_bchN']
        # Cobalt = ['cobN', 'cobT', 'cobS']
        Vitamin_B12 = ['cobD_cobC', 'cobO', 'cobP_cobU', 'cobQ', 'cobV_cobS', 'cysG_cobA']
        Vitamin_B12_genes, Chlorophyll_genes, all_genes_in_org = [], [], []

        c_org = Org.objects.get(param_set__name='chel_genotype_genes', id=id_org)


        if from_db:
            feature_12 = c_org.prm['chel_synthesis_chl']
            feature_34 = c_org.prm['chel_synthesis_b12']

        else:
            f_org = ChelataseFeat.objects.filter(seq__org=c_org).all()
            for i in f_org:
                if 'chel_pathway' in i.prm.keys() and i.prm.chel_evalue < 1e-50:
                    if i.prm.chel_pathway == "Chlorophyll" and i.prm.chel_gene_group in Chlorophyll:
                        Chlorophyll_genes.append(i.prm.chel_gene_group)
                    elif i.prm.chel_pathway == "B12" and i.prm.chel_gene_group in Vitamin_B12:
                        Vitamin_B12_genes.append(i.prm.chel_gene_group)

                    all_genes_in_org.append(i.prm.chel_gene_group)

            feature_1_chlH = len([gen for gen in all_genes_in_org if 'chlH' in gen]) >= 1
            feature_2_Chlorophyll = len(Chlorophyll_genes) / len(Chlorophyll) >= 0.5
            feature_3_cobN = len([gen for gen in all_genes_in_org if 'cobN' in gen]) >= 1
            feature_4_B12 = len(Vitamin_B12_genes) / len(Vitamin_B12) >= 0.5

            feature_12 = feature_1_chlH * feature_2_Chlorophyll
            feature_34 = feature_3_cobN * feature_4_B12

            if  'Archaea' in c_org.prm['taxonomy']:
                feature_12 *=0  # Archaea is not able to synthesize Chlorophyll

        def get_prediction(feature_chl, feature_b12):
            if feature_chl:
                feature_chl = "YES"
            else:
                feature_chl = 'NO'
            if feature_b12:
                feature_b12 = "YES"
            else:
                feature_b12 = 'NO'

            return 'Chlorophyll biosynthesis - {}   |   Vitamin B12 biosynthesis {}'.format(
                feature_chl, feature_b12)

        if output2text: return get_prediction(feature_12, feature_34)
        else: return feature_12, feature_34


    def classify_org(self):
        """Returns possible function of the chlD gene in the current org
        """
        c_org = Org.objects.get(param_set__name='chel_genotype_genes', id=self.id)

        num_fs_chlD = c_org.prm['num_fs_chlD']
        num_cobN = c_org.prm['num_cobN']
        num_cobT = c_org.prm['num_cobT']
        num_cobS = c_org.prm['num_cobS']

        num_chlH = c_org.prm['num_chlH']
        num_chlD = c_org.prm['num_chlD']
        num_chlI = c_org.prm['num_chlI']

        if num_fs_chlD == 0:  # there are no shifts
            if num_cobN * num_cobT * num_cobS > 0 and num_chlH * num_chlD * num_chlI > 0:  # M_Co_chel_and_M_Mg_chel
                return self.CHLD_FUNCTIONS['M_Mg_but_not_Co']
            elif num_cobN * num_cobT * num_cobS > 0 and num_chlH * num_chlD * num_chlI == 0:
                return self.CHLD_FUNCTIONS['M_Mg'] # 'COchel.png' # M_Co_chel
            elif num_cobN * num_cobT * num_cobS == 0 and num_chlH * num_chlD * num_chlI > 0:
                return self.CHLD_FUNCTIONS['M_Co']
            # TODO: find this case!
            #elif num_chlH * num_chlD * num_chlI > 0 and num_cobN > 0 and num_cobT * num_cobS == 0:
            #    return 'chelTable2.png'
            elif num_cobN > 0 and num_chlD > 0 and num_chlI > 0 and num_cobT * num_cobS * num_chlH == 0:
<<<<<<< HEAD
                return self.CHLD_FUNCTIONS['M_Co']
        else:
=======
                url = 'chelTable3.png'  # Co_chelatase
        else: # there are shifts
>>>>>>> 95bce1f3cef9e3a4fca47459664fdcdffa5f98f7
            if num_cobN > 0 and num_chlD > 0 and num_cobT * num_cobS * num_chlH * num_chlI == 0:
                return self.CHLD_FUNCTIONS['M_and_S_Co_fs']
            else:
                return {}  # other / unknown

        return {}

def _get_or_create_parent_feat_from_tblastn_hits(user, seq, hsp_n, hsp_c):
    """For a given pair of tBLASTn hits create the parent feature that
    corresponds to the upstream part of fsCDS (i.e. the shorter product).
    """
    half_len_n = (hsp_n.sbjct_end-hsp_n.sbjct_start)/2
    hsp_n_all_feats = get_overlapping_feats_from_list(
        seq.record.features, hsp_n.sbjct_start, hsp_n.sbjct_end,
        hsp_n.sbjct_strand, min_overlap=half_len_n, all_types=['CDS'])

    half_len_c = (hsp_c.sbjct_end-hsp_c.sbjct_start)/2
    hsp_c_all_feats = get_overlapping_feats_from_list(
        seq.record.features, hsp_c.sbjct_start, hsp_c.sbjct_end,
        hsp_c.sbjct_strand, min_overlap=half_len_c, all_types=['CDS'])

    parent = None
    if len(hsp_n_all_feats) == 1 and len(hsp_c_all_feats) == 1:
        hsp_n_feat = hsp_n_all_feats[0]
        hsp_c_feat = hsp_c_all_feats[0]
        if hsp_n_feat == hsp_c_feat:
            # The parent CDS has to be created manually because the
            # frameshifted gene annotated as a single
            # frameshifted gene without translation, i.e.
            #      hsp_n            hsp_c
            #   ---------         ------------
            # ====================================>
            #   hsp_n_feat          hsp_c_feat
            parent = ChelataseFeat.get_or_create_from_frameshifted_SeqFeature(
                user, seq, hsp_n_feat)
        else:
            # Parent CDS can be created from the annotated feature
            # because the frameshifted gene annotated as two genes, i.e.:
            #      hsp_n            hsp_c
            #   ---------         ------------
            # ==============>  ====================>
            #   hsp_n_feat          hsp_c_feat
            parent = ChelataseFeat.get_or_create_from_SeqFeature(
                user, seq, hsp_n_feat)
    else:
        logging.error("No annotated features for tBLASTn hit(s)!")

    return parent

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

    l_ss = get_stop_stop_seq(lh.sbjct_start, lh.sbjct_end, lh.sbjct_strand,
                              chr_seq, gencode)
    r_ss = get_stop_stop_seq(rh.sbjct_start, rh.sbjct_end, rh.sbjct_strand,
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

def _make_hits_dict(hsp_list):
    """Returns a dict where the keys are the hit sequence IDs and
    the values are the lists of Bio.Blast.Record.HSP objects.
    """
    all_hits = {}
    for hsp in hsp_list:
        all_hits.setdefault(hsp.sbjct_id, []).append(hsp)
    return all_hits

def _org_name_to_species(org_name):
    """Extracts species (the first two words) form the org name, e.g.:
    'Methanocaldococcus sp. FS406-22'  => 'Methanocaldococcus sp'
    """
    sprecies_re = re.compile(r'^[a-z]+ [a-z]+', re.IGNORECASE)
    mo = sprecies_re.search(org_name)
    if mo is not None:
        return mo.group()
    else:
        return None

