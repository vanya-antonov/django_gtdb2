# Copyright 2018 by Ivan Antonov. All rights reserved.


from django.db import models
from django.db.models import signals
from django.dispatch import receiver

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.fshift import Fshift


class Cof(AbstractUnit):
    type = models.CharField(max_length=255)
    fshift_set = models.ManyToManyField(Fshift, db_table='cof_fshifts')

    class Meta:
        db_table = 'cofs'

    @classmethod
    def create_from_seed_fshifts(cls, user, seed_fshifts, name=None,
                                 descr=None, type='cof'):
        """Creates a new COF in db, associates it with the provided seed
        fshifts and sets the seed flag to them.
        """
        if len(seed_fshifts) == 0:
            raise ValueError("The list of seed fshifts is empty!")

        # Each fshift can belong to one cof only
        for fshift in seed_fshifts:
            if fshift.cof_set.count() > 0:
                raise ValueError("Fshift '%s' already belongs to a COF" %
                                 fshift.name)

        cof = cls(user=user, name=name, descr=descr, type=type)
        cof.save()
        cof.fshift_set.set(seed_fshifts)

        # Add the seed flag to each fshift
        for fshift in seed_fshifts:
            fshift.seed = True
            fshift.save()

        return cof


class CofParam(AbstractParam):
    parent = models.ForeignKey(Cof, on_delete=models.CASCADE,
                               related_name='param_set')

    class Meta:
        db_table = 'cof_params'


@receiver(signals.pre_delete, sender=Cof)
def on_cof_delete(sender, instance, using, **kwargs):
    """Make sure to remove the seed flag from all the cof seed fshifts.
    """
    for fshift in instance.fshift_set.all():
        fshift.seed = False
        fshift.save()

