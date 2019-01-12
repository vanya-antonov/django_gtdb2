# Copyright 2018 by Ivan Antonov. All rights reserved.

from django.db import models

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.seq import Seq
from gtdb2.models.fshift import Fshift

from Bio.SeqFeature import FeatureLocation, CompoundLocation


class Feat(AbstractUnit):
    seq = models.ForeignKey(Seq, on_delete=models.CASCADE)
    type = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()
    origin = models.CharField(max_length=255)
    fshift_set = models.ManyToManyField(Fshift, db_table='feat_fshifts')

    class Meta:
        db_table = 'feats'

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    prm_info = dict(list(AbstractUnit.prm_info.items()) + list({
        'gene': {},
        'gene_synonym': {},
        'protein_id': {},
        'ribosomal_slippage': {},
        'translation': {'value_attr': 'data'},
    }.items()))

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

        return cls.create_from_SeqFeature(user, seq, f)

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

        self.make_all_params(f)

        return self

    def make_all_params(self, f=None):
        if self.origin != 'annotation':
            raise NotImplementedError("Feat make_all_params!")

        if f is None:
            f = get_overlapping_feats_from_record(
                self.seq.record, self.start, self.end, self.strand,
                all_types=[self.type], min_overlap=1, max_feats=1)[0]

        self._make_param_from_q(f, 'gene')
        self._make_param_from_q(f, 'gene_synonym')

        if self.type == 'CDS':
            self._make_param_cds(f)
        elif self.type == 'rRNA':
            self._make_param_rrna(f)

    def _make_param_from_q(self, f, name):
        "Create param based on gbk qualifier."
        for v in f.qualifiers.get(name, []):
            self.add_param(name, v)

    def _make_param_cds(self, f):
        self._make_param_from_q(f, 'protein_id')
        self._make_param_from_q(f, 'ribosomal_slippage')

        for prot_seq in f.qualifiers.get('translation', []):
            self.add_param('translation', data=prot_seq, num=len(prot_seq))

    def _make_param_rrna(self, f):
        # see add_16S_rRNA() from
        # ~/_my/Programming/python/scripts/frameshift/gtdb2_manager.py
        raise NotImplementedError("rRNA feat!!")


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

class FeatParam(AbstractParam):
    parent = models.ForeignKey(Feat, on_delete=models.CASCADE,
                               related_name='param_set')

    class Meta:
        db_table = 'feat_params'

