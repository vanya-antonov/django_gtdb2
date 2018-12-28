# Copyright 2018 by Ivan Antonov. All rights reserved.

from collections import Counter
import os

from django.db import models

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.org import Org


class Seq(AbstractUnit):
    # Overwrite AbstractUnit attribute - id is a string
    id = models.CharField(max_length=255, primary_key=True)

    # Overwrite AbstractUnit attribute - name must be unique and NOT NULL
    name = models.CharField(max_length=255, unique=True)

    org = models.ForeignKey(Org, on_delete=models.CASCADE)
    type = models.CharField(max_length=255)
    len = models.IntegerField()

    class Meta:
        db_table = 'seqs'

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    prm_info = dict(list(AbstractUnit.prm_info.items()) + list({
        'transl_table': {'type_fun': int},
    }.items()))

    @property
    def seq(self):
        "Returns Bio.Seq object."
        if not hasattr(self, '_seq'):
            if hasattr(self, '_record'):
                self._seq = self._record.seq
            else:
                # Read seq from fasta file
                self._seq = self.org.read_seq_file(
                    self.id, seq_dir='seq_fna').seq
        return self._seq

    @property
    def record(self):
        "Reads gbk_fn and returns Bio.SeqRecord object."
        if not hasattr(self, '_record'):
            self._record = self.org.read_seq_file(
                self.id, seq_dir='seq_gbk')
        return self._record

    @property
    def transl_table(self):
        "Returns transl_table ID (int)."
        if 'transl_table' in self.prm:
            return self.prm['transl_table']
        else:
            return self.org.transl_table

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
        record = org.read_seq_file(ext_id, seq_dir='seq_gbk')
        seq = cls(user=user, org=org,
                  id=record.id, name=record.name, descr=record.description,
                  type=record.annotations['molecule_type'], len=len(record.seq))
        seq.save()

        seq.make_all_params()

        return seq

    def make_all_params(self):
        self._make_param_transl_table()

    def _make_param_transl_table(self):
        "Analyzes all annotated CDSs and creates transl_table param."
        all_transl_f = list(filter(lambda f: 'transl_table' in f.qualifiers,
                                   self.record.features))

        if len(all_transl_f) == 0:
            logging.warning("Sequence '%s' doesn't have transl_table annotations" %
                            self.id)
            return None

        all_transl_tables = [f.qualifiers['transl_table'][0] for f in all_transl_f]

        # https://stackoverflow.com/a/10797913/310453
        transl_counts = Counter(all_transl_tables).most_common()
        if len(transl_counts) > 1:
            raise ValueError("Different genetic codes are annotated for seq "
                            "'%s': '%s'" % (self.id, transl_counts))
        return self.set_param('transl_table', transl_counts[0][0])

class SeqParam(AbstractParam):
    parent = models.ForeignKey(Seq, on_delete=models.CASCADE,
                               related_name='param_set')

    class Meta:
        db_table = 'seq_params'

