# Copyright 2018 by Ivan Antonov. All rights reserved.

#from chelatase_db.models.cof import ChelataseCof
from gtdb2.models.fshift import Fshift


class ChelataseFshift(Fshift):
    class Meta:
        proxy = True

    # Add info about chelatase specific params
    # Fshift.prm_info.update({})

    # Overwrite the computed property to keep using the parent's param table
    param_set = property(fget=lambda self: self.orgparam_set)

