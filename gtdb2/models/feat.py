# Copyright 2018 by Ivan Antonov. All rights reserved.

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
        'protein_id': {},
        'ribosomal_slippage': {},
        'translation': {'value_attr': 'data'},
    }.items()))

    @property
    def feature(self):
        "Retruns a Bio.SeqFeature.SeqFeature object."
        if self.type == 'fsCDS':
            all_fshifts = list(self.fshift_set.all())
            stop_coord = self.end if self.strand == 1 else self.start
            loc = _make_compound_location(self.parent, all_fshifts, stop_coord)
        else:
            loc = FeatureLocation(self.start, self.end, self.strand)
        return SeqFeature(loc, type=self.type)

    @classmethod
    def get_or_create_from_gbk_annotation(cls, user, seq, start, end, strand,
                                         f_type='CDS'):
        """Loads features from GenBank file that overlap with the given
        target region.
        """
        overlapping_feats = get_overlapping_feats_from_record(
            seq.record, start, end, strand, all_types=[f_type],
            min_overlap=1, max_feats=1)

        if len(overlapping_feats) == 0:
            return None
        f = overlapping_feats[0]

        # Check if this feature already present in DB
        feat = cls.objects.filter(
            seq=seq, type=f.type, strand=f.location.strand,
            start=int(f.location.start), end=int(f.location.end)
        ).first()
        if feat is not None:
            return feat

        return cls.create_from_SeqFeature(user, seq, f,
                                          origin='gbk_annotation')

    @classmethod
    def create_from_SeqFeature(cls, user, seq, f, origin):
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
        if self.type == 'fsCDS':
            self._make_fscds_name()
            self._make_param_fscds_translation()

        if self.origin == 'gbk_annotation':
            self._make_all_params_gbk()

    def _make_fscds_name(self):
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

    def _make_param_fscds_translation(self):
        fscds_nt = self.feature.extract(self.seq.seq)
        fscds_aa = fscds_nt.translate(table=self.seq.transl_table)

        # Remove possible stop codon at the end
        fscds_aa = fscds_aa.rstrip('*')
        self.add_param('translation', data=fscds_aa, num=len(fscds_aa))

    def _make_all_params_gbk(self):
        # f is a Bio.SeqFeature.SeqFeature object
        f = get_overlapping_feats_from_record(
            self.seq.record, self.start, self.end, self.strand,
            all_types=[self.type], min_overlap=1, max_feats=1)[0]

        self._make_param_from_q(f, 'gene')
        self._make_param_from_q(f, 'gene_synonym')

        if self.type == 'CDS':
            self._make_param_from_q(f, 'protein_id')
            self._make_param_from_q(f, 'ribosomal_slippage')

            for prot_seq in f.qualifiers.get('translation', []):
                self.add_param('translation', data=prot_seq, num=len(prot_seq))
        elif self.type == 'rRNA':
            # see add_16S_rRNA() from
            # ~/_my/Programming/python/scripts/frameshift/gtdb2_manager.py
            raise NotImplementedError("rRNA feat!!")

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

def get_overlapping_feats_from_record(
    record, start, end, strand, all_types=['CDS'], min_overlap=1, max_feats=None):
    """    _overlap_len attribute will be added to returned features
    """
    target_loc = FeatureLocation(start, end, strand)
    overlapping_feats = []
    for f in record.features:
        if f.type not in all_types:
            continue
        f._overlap_len = get_FeatureLocation_overlap_len(f.location, target_loc)
        if f._overlap_len >= min_overlap:
            overlapping_feats.append(f)

    if max_feats is not None and len(overlapping_feats) > max_feats:
        # Sort by overlap_len: https://docs.python.org/3.6/howto/sorting.html
        overlapping_feats = sorted(overlapping_feats, reverse=True,
                                   key=lambda f: f._overlap_len)
        # Take the longest feats only
        overlapping_feats = overlapping_feats[0:max_feats]

    return overlapping_feats

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

