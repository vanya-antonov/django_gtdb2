# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint
import logging

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from django.db import models

from gtdb2.lib.bio import is_stop
from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.seq import Seq
from gtdb2.models.fshift import Fshift


class Feat(AbstractUnit):
    seq = models.ForeignKey(Seq, on_delete=models.CASCADE)
    type = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()
    origin = models.CharField(max_length=255, choices=Fshift.ORIGIN_CHOICES)
    fshift_set = models.ManyToManyField(Fshift, db_table='feat_fshifts')
    parent = models.ForeignKey('self', on_delete=models.CASCADE,
                               default=None, blank=True, null=True)

    class Meta:
        db_table = 'feats'

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(AbstractUnit.PRM_INFO.items()) + list({
        'gene': {},
        'gene_synonym': {},
        'location_str': {},
        'protein_id': {},
        'ribosomal_slippage': {},
        'seq_nt': {'value_attr': 'data'},
        'translation': {'value_attr': 'data'},
    }.items()))

    @property
    def location(self):
        """Returns a Bio.SeqFeature.FeatureLocation object
        https://github.com/biopython/biopython/blob/master/Bio/SeqFeature.py
        """
        if self.type == 'fsCDS':
            all_fshifts = list(self.fshift_set.all())
            stop_coord = self.end if self.strand == 1 else self.start
            loc = _make_compound_location(self.parent, all_fshifts, stop_coord)
        else:
            loc = FeatureLocation(self.start, self.end, self.strand)
        return loc

    @property
    def feature(self):
        "Retruns a Bio.SeqFeature.SeqFeature object."
        return SeqFeature(self.location, type=self.type)

    @classmethod
    def get_or_create_from_gbk_annotation(cls, user, seq, start, end, strand,
                                         f_type='CDS'):
        """Loads features from GenBank file that overlap with the given
        target region.
        """
        if f_type == 'CDS':
            all_types = ['CDS', 'fsCDS']
        else:
            all_types = [f_type]

        half_len = (end-start)/2

        # Check if a feature already exists in DB for the specified region
        db_feats = get_overlapping_feats_from_list(
            seq.feat_set.all(), start, end, strand,
            min_overlap=half_len, all_types=all_types)
        if len(db_feats) > 0:
            return db_feats[0]

        overlapping_feats = get_overlapping_feats_from_list(
            seq.record.features, start, end, strand,
            min_overlap=half_len, all_types=all_types)

        if len(overlapping_feats) == 0:
            return None
        f = overlapping_feats[0]
        return cls.get_or_create_from_SeqFeature(
            user, seq, f)

    @classmethod
    def get_or_create_from_SeqFeature(cls, user, seq, f, **kwargs):
        """The argument f is a Bio.SeqFeature.SeqFeature object.
        """
        # Make sure this feature is not present in the DB
        feat = cls.objects.filter(
            seq=seq, type=f.type, strand=f.location.strand,
            start=int(f.location.start), end=int(f.location.end)
        ).first()
        if feat is None:
            feat = cls.create_from_SeqFeature(user, seq, f, **kwargs)
        return feat

    @classmethod
    def create_from_SeqFeature(cls, user, seq, f, origin='annotation'):
        """The argument f is a Bio.SeqFeature.SeqFeature object.
        """
        self = cls(user=user, seq=seq, type=f.type,
                   start=int(f.location.start),
                   end=int(f.location.end),
                   strand=f.location.strand,
                   name=f.qualifiers.get('locus_tag', [None])[0],
                   descr=f.qualifiers.get('product', [None])[0],
                   origin=origin)
        self.save()

        self.create_all_params()

        return self

    @classmethod
    def get_or_create_fscds_from_parent(cls, user, parent_cds, all_fshifts):
        """Arguments:
         - parent_cds - Feat object corresponding to the CDS where
           translation initiates
         - all_fshifts - a list of all unique frameshifts that are required
           for the translation of fsCDS to be created.
        """
        # Check if this fsCDS already exists in DB using the
        # annotation approach: https://stackoverflow.com/a/8637972/310453
        fscds = (cls.objects
                    .filter(type='fsCDS', seq=parent_cds.seq)
                    .filter(parent=parent_cds)
                    .filter(fshift_set__in=all_fshifts)
                    .annotate(num_fshifts=models.Count('fshift_set'))
                    .filter(num_fshifts=len(all_fshifts))
                    .first())
        if fscds is not None:
            return fscds

        return cls.create_fscds_from_parent(user, parent_cds, all_fshifts)

    @classmethod
    def create_fscds_from_parent(cls, user, parent_cds, all_fshifts):
        "See get_or_create_fscds_from_parent() for more info."
        if len(all_fshifts) == 0:
            raise ValueError("fsCDS must have at least one fshift!")

        strand = parent_cds.strand
        start_codon = parent_cds.start if strand == 1 else parent_cds.end
        for fshift in all_fshifts:
            if fshift.seq != parent_cds.seq:
                raise ValueError("All fshifts must be on the same seq as "
                                "the parent CDS!")
            if fshift.strand != parent_cds.strand:
                raise ValueError("All fshifts must be on the same strand as "
                                "the parent CDS!")
            # TODO: implement is_downstream()
            #if not is_downstream(start_codon, fshift.coord, strand)
            #    raise ValueError("All fshifts must be located downstream "
            #                     "of the parent start codon!")

        c_loc = _make_compound_location(parent_cds, all_fshifts)

        # Make sure the fsCDS does not have internal stop codons
        fscds_f = SeqFeature(c_loc, type="CDS")
        fscds_nt = fscds_f.extract(parent_cds.seq.seq)
        fscds_aa = fscds_nt.translate(table=parent_cds.seq.transl_table)
        if '*' in fscds_aa.rstrip('*'):
            raise ValueError("fsCDS translation contains internal stop "
                             "codon(s):\n%s" % fscds_aa)

        fscds = cls(
            user=user, parent=parent_cds, origin='genetack',
            type='fsCDS', seq=parent_cds.seq,
            start=int(c_loc.start), end=int(c_loc.end), strand=c_loc.strand)
        fscds.save()

        fscds.fshift_set.set(all_fshifts)

        fscds.create_all_params()

        return fscds

    def create_all_params(self):
        """Create feature params depending on its type (CDS, fsCDS, rRNA etc)
        and origin (genbank, gtdb).
        """
        self._make_param_location_str()

        #TODO: if self.origin == 'gtdb':
        if self.origin == 'genetack':
            if self.type == 'fsCDS':
                self._make_param_gtdb_fscds()
            elif self.type == 'CDS':
                raise NotImplementedError("Parent GTDB CDS feat!!")
                self._make_param_gtdb_cds()
            else:
                raise ValueError("Unknonwn feature type = '%s'" % self.type)
        #TODO: if self.origin == 'genbank':
        elif self.origin == 'annotation':
            # f is a Bio.SeqFeature.SeqFeature object
            f = get_overlapping_feats_from_list(
                self.seq.record.features, self.start, self.end, self.strand,
                all_types=[self.type], min_overlap=len(self.feature)/2
            )[0]
            if self.type == 'CDS':
                self._make_param_gbk_cds(f)
            elif self.type == 'fsCDS':
                raise NotImplementedError("Annotated fsCDS feat!!")
                self._make_param_gbk_fscds(f)
            elif self.type == 'rRNA':
                self._make_param_gbk_rrna(f)
            else:
                raise ValueError("Unknonwn feature type = '%s'" % self.type)
        else:
            raise ValueError("Unknonwn feature origin = '%s'" % self.origin)

    def _make_param_location_str(self):
        """Creates location_str param like 'NC_000911.2:123-456(+)'.
        """
        strand_str = '+' if self.strand == 1 else '-'
        location_str = '%s:%d-%d(%s)' % (
            self.seq.id, self.start, self.end, strand_str)
        self.set_param('location_str', location_str)

    def _make_param_gtdb_fscds(self):
        """Creates params for predicted and not-annotated fsCDS feat.
        """
        self._make_param_fscds_name()
        self._make_param_fscds_seqs()

    def _make_param_fscds_name(self):
        "Creates fsCDS name like 'MEFER_RS06095_fs1122957'."
        if self.parent.name is None:
            all_parts = ['fsCDS']
        else:
            all_parts = [self.parent.name]

        # Sort frameshifts in the order of their appearance during translation
        all_fshifts = sorted(self.fshift_set.all(), key=lambda fs: fs.coord,
                             reverse=(self.strand == -1))
        for fs in all_fshifts:
            all_parts.append('fs' + str(fs.coord))

        self.name = '_'.join(all_parts)
        self.save()

    def _make_param_fscds_seqs(self):
        fscds_nt = self.feature.extract(self.seq.seq)
        fscds_aa = fscds_nt.translate(table=self.seq.transl_table)

        # Remove possible stop codon at the end
        fscds_aa = fscds_aa.rstrip('*')
        self.add_param('translation', data=fscds_aa, num=len(fscds_aa))
        self.add_param('seq_nt', data=fscds_nt, num=len(fscds_nt))

    def _make_param_gbk_rrna(self, f):
        """f is a Bio.SeqFeature.SeqFeature object.
        """
        self._make_param_seq_nt()

    def _make_param_gbk_cds(self, f):
        """f is a Bio.SeqFeature.SeqFeature object.
        """
        self._make_param_from_q(f, 'gene')
        self._make_param_from_q(f, 'gene_synonym')
        self._make_param_from_q(f, 'protein_id')
        self._make_param_from_q(f, 'ribosomal_slippage')

        self._make_param_seq_nt()
        self._make_param_translation(f)

    def _make_param_seq_nt(self):
        """Extracts the nt sequence of the feature and saves it in
        the params table.
        """
        # self.feature is a Bio.SeqFeature.SeqFeature object
        seq_nt = self.feature.extract(self.seq.seq).upper()
        self.add_param('seq_nt', data=seq_nt, num=len(seq_nt))

    def _make_param_translation(self, f):
        """f is a Bio.SeqFeature.SeqFeature object.
        """
        # Get the seqs
        cds_nt = f.extract(self.seq.seq).upper()

        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc25
        try:
            cds_aa = cds_nt.translate(
                table=self.seq.transl_table, to_stop=True, cds=True
            ).upper()
        except:
            # In-frame stop codon will cause the KeyError
            logging.error(
                "Can't translate CDS feat '%s' from seq '%s':\n%s" %
                (f.location, self.seq.id, cds_nt))
            cds_aa = None

        # Some checks
        if len(cds_nt) % 3 != 0:
            logging.error(
                "CDS sequence len is not divisible by 3 for feature:\n%s" % f)

        if 'translation' in f.qualifiers:
            translation = f.qualifiers['translation'][0].upper()
            if cds_aa is None:
                cds_aa = translation

            if translation != cds_aa:
                logging.error(
                    "The annotated and generated translations do not match "
                    "for feature '%s':\n\n%s\n%s\n\n" %
                    (f, translation, cds_aa))
        else:
            logging.error(
                "The feature doesn't have annotated translation:\n%s" % f)

        if cds_aa is not None:
            self.add_param('translation', data=cds_aa, num=len(cds_aa))

    def _make_param_from_q(self, f, name):
        "Create param based on gbk qualifier."
        for v in f.qualifiers.get(name, []):
            self.add_param(name, v)


class FeatParam(AbstractParam):
    parent = models.ForeignKey(Feat, on_delete=models.CASCADE,
                               related_name='param_set')

    class Meta:
        db_table = 'feat_params'


def _make_compound_location(parent_cds, all_fshifts, stop_coord=None):
    """The parent_cds and all_fshifts are assumed to be on the same seq
    and strand.
    Arguments:
     - stop_coord - the coordinate of the location end (for plus strand) or
       location start (for minus strand)
    """
    strand = parent_cds.strand

    # Sort frameshifts in the order of their appearance during translation
    all_fshifts = sorted(all_fshifts, key=lambda fs: fs.coord,
                         reverse=(strand == -1))

    # Get the translation initiation (start codon) coordinate
    current_coord = parent_cds.start if strand == 1 else parent_cds.end
    all_locs = []
    for fshift in all_fshifts:
        # Make the FeatureLocation for this region, shift the frame and
        # continue to the next fshift (if any)
        if strand == 1:
            # Moving from left to right
            loc = FeatureLocation(current_coord, fshift.coord, strand)
            current_coord = fshift.coord + fshift.len
        else:
            # Moving from right to left
            loc = FeatureLocation(fshift.coord, current_coord, strand)
            current_coord = fshift.coord - fshift.len
        all_locs.append(loc)

    # After all shifts of the reading frame, find the stop codon
    if stop_coord is None:
        stop_coord = _find_downstream_stop(
            parent_cds.seq.seq, current_coord, strand,
            parent_cds.seq.transl_table)

    if strand == 1:
        # Moving from left to right
        stop_loc = FeatureLocation(current_coord, stop_coord, strand)
    else:
        # Moving from right to left
        stop_loc = FeatureLocation(stop_coord, current_coord, strand)
    all_locs.append(stop_loc)

    return CompoundLocation(all_locs)

def _find_downstream_stop(seq, start_coord, strand, gcode):
    """Retruns coord AFTER the stop codon or start/end of sequence."""
    coord = start_coord
    if strand == 1:
        # Move to the right
        while coord <= len(seq)-3:
            codon = seq[coord:(coord+3)]
            coord += 3  # move to the next codon
            if is_stop(codon, strand, gcode):
                break
    else:
        # Move to the left
        while coord >= 3:
            codon = seq[(coord-3):coord]
            coord -= 3  # move to the next codon
            if is_stop(codon, strand, gcode):
                break

    return coord

def get_overlapping_feats_from_list(all_feats, start, end, strand,
                                    all_types=None, min_overlap=1):
    """Returns a list of objects from f_list that overlap with the specified
    genomic location. The returned objects will be sorted (descending) by the
    value of the '_overlap_len' attribute (int) that will be added to each
    object.

    Arguments:
     - all_feats - a list of feature-objects (e.g. Bio.SeqFeature.SeqFeature or
     gtdb2.models.feat.Feat) that have the .type (string) and the .location
     (Bio.SeqFeature.FeatureLocation) attributes.
    """
    target_loc = FeatureLocation(start, end, strand)
    overlapping_feats = []
    for f in all_feats:
        if all_types is not None and f.type not in all_types:
            continue

        if f.location.end < start or end < f.location.start:
            continue

        f._overlap_len = get_FeatureLocation_overlap_len(f.location, target_loc)
        if f._overlap_len >= min_overlap:
            overlapping_feats.append(f)

    # Sort by overlap_len: https://docs.python.org/3.6/howto/sorting.html
    return sorted(overlapping_feats, reverse=True,
                  key=lambda f: f._overlap_len)

def get_FeatureLocation_overlap_len(f1, f2):
    # https://github.com/biopython/biopython/issues/896
    if not isinstance(f1, (FeatureLocation, CompoundLocation)) or \
       not isinstance(f2, (FeatureLocation, CompoundLocation)):
        raise ValueError("Wrong feature types: %s / %s" %
                         (type(f1), type(f2)))

    if f1.strand is not None and \
       f2.strand is not None and \
       f1.strand != f2.strand:
        return 0

    return len(set(f1).intersection(set(f2)))

