from django.db import models

from .abstract import AbstractUnit, AbstractParam
from .user import User
from .org import Org


class Seq(AbstractUnit):
    id = models.CharField(max_length=255, primary_key=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    org = models.ForeignKey(Org, on_delete=models.CASCADE)
    type = models.CharField(max_length=255)
    ext_id = models.CharField(max_length=255)
    len = models.IntegerField()

    class Meta:
        db_table = 'seqs'


class SeqParam(AbstractParam):
    seq = models.ForeignKey(Seq, on_delete=models.CASCADE)

    class Meta:
        db_table = 'seq_params'

