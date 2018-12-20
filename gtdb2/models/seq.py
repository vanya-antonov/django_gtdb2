# Copyright 2018 by Ivan Antonov. All rights reserved.

import os

from django.db import models

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.org import Org


class Seq(AbstractUnit):
    id = models.CharField(max_length=255, primary_key=True)
    org = models.ForeignKey(Org, on_delete=models.CASCADE)
    type = models.CharField(max_length=255)
    len = models.IntegerField()

    class Meta:
        db_table = 'seqs'

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    prm_info = dict(list(AbstractUnit.prm_info.items()) + list({
        'gbk_fn': {},
    }.items()))

    def _get_record(self):
        if not hasattr(self, '_record'):
            gbk_path = self.gtdb.get_full_path_to(self.prm['gbk_fn'])
            self._record = SeqIO.read(gbk_path, "genbank")
        return self._record
    record = property(
        fget=_get_record,
        doc="Reads gbk_fn and returns SeqRecord object.")

    @classmethod
    def get_or_create_from_ext_id(cls, user, org, ext_id):
        "Returns existing or a newly created Seq object."
        seq = cls.objects.filter(pk=ext_id).first()
        if seq is None:
            seq = cls.create_from_ext_id(user, org, ext_id)
        return seq

    @classmethod
    def create_from_ext_id(cls, user, org, ext_id):
        "Reads gbk file from org_dir and creates seq from SeqRecord."
        gbk_fn = os.path.join(org.prm['dir_path'], 'seq_gbk', ext_id)
        record = SeqIO.read(cls.gtdb.get_full_path_to(gbk_fn), "genbank")
        seq = cls(user=user, org=org,
                  id=record.id, name=record.name, descr=record.description,
                  type=record.annotations['molecule_type'], len=len(record.seq))
        seq.save()
        seq.add_param('gbk_fn', gbk_fn)
        return seq


class SeqParam(AbstractParam):
    parent = models.ForeignKey(Seq, on_delete=models.CASCADE,
                               related_name='param_set')

    class Meta:
        db_table = 'seq_params'

