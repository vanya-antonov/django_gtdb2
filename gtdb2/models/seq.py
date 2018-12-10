# Copyright 2018 by Ivan Antonov. All rights reserved.

from django.db import models

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from .org import Org


class Seq(AbstractUnit):
    id = models.CharField(max_length=255, primary_key=True)
    org = models.ForeignKey(Org, on_delete=models.CASCADE)
    type = models.CharField(max_length=255)
    ext_id = models.CharField(max_length=255)
    len = models.IntegerField()

    class Meta:
        db_table = 'seqs'

    def as_SeqRecord(self, with_seq=False):
        "Returns SeqRecord object with seq from fasta file if requested."
        if with_seq:
            fasta_fn = self.org.get_full_path_to('seq_fna', self.id)
            seq = SeqIO.read(fasta_fn, "fasta").seq
        else:
            seq = Bio.Seq.UnknownSeq(self.len)
        return SeqRecord(
            seq, id=self.id, name=self.name, description=self.descr,
            annotations={'molecule_type': self.type})

    @classmethod
    def create_from_ext_id(cls, user, org, ext_id):
        "Reads gbk file from org_dir and creates seq from SeqRecord."
        gbk_fn = org.get_full_path_to('seq_gbk', ext_id)
        record = SeqIO.read(gbk_fn, "genbank")
        seq = cls(user=user, org=org,
                  id=record.id, name=record.name, descr=record.description,
                  type=record.annotations['molecule_type'], len=len(record.seq))
        seq.save()
        return seq


class SeqParam(AbstractParam):
    parent = models.ForeignKey(Seq, on_delete=models.CASCADE)

    class Meta:
        db_table = 'seq_params'

