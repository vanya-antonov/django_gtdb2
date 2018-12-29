# Copyright 2018 by Ivan Antonov. All rights reserved.

from django.db import models

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.seq import Seq
from gtdb2.models.fshift import Fshift


class Feat(AbstractUnit):
    seq = models.ForeignKey(Seq, on_delete=models.CASCADE)
    type = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()
    fshift_set = models.ManyToManyField(Fshift, db_table='feat_fshifts')

    class Meta:
        db_table = 'feats'

    @classmethod
    def get_or_create_from_locus_annotation(user, seq, start, end, strand):
        # see create_new_in_db_from_gbk_for_region() function in 
        # ~/_my/Programming/python/lib/mylib/genetackdb2.py
        raise NotImplementedError("Feat!")

    @classmethod
    def create_from_fshift(cls, user, fshift):
        raise NotImplementedError("Feat!")
        self = cls.get_or_create_from_locus_annotation(
            user, fshift.seq, fshift.start, fshift.end, fshift.strand)


class FeatParam(AbstractParam):
    parent = models.ForeignKey(Feat, on_delete=models.CASCADE,
                               related_name='param_set')

    class Meta:
        db_table = 'feat_params'

