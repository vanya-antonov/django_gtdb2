# Copyright 2018 by Ivan Antonov. All rights reserved.


from django.db import models

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.fshift import Fshift


class Cof(AbstractUnit):
    type = models.CharField(max_length=255)
    fshift_set = models.ManyToManyField(Fshift, db_table='cof_fshifts')

    class Meta:
        db_table = 'cofs'


class CofParam(AbstractParam):
    parent = models.ForeignKey(Cof, on_delete=models.CASCADE)

    class Meta:
        db_table = 'cof_params'

