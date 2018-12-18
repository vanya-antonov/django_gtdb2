# Copyright 2018 by Ivan Antonov. All rights reserved.

from chelatase_db.models.cof import ChelataseCof
from gtdb2.models.org import Org


class ChelataseOrg(Org):
    class Meta:
        proxy = True

    # Merge two dicts: https://stackoverflow.com/a/38990/310453
    prm_info = dict(list(Org.prm_info.items()) + list({
        'num_chld_fshifts': {'prm_attr': 'num', 'type_fun': int},
        'num_chld_feats': {'prm_attr': 'num', 'type_fun': int},
    }.items()))

    # Overwrite the computed property to keep using the parent's param table
    param_set = property(fget=lambda self: self.orgparam_set)

    def create_chld_fshifts_and_feats(self, user):
        chld_cof = ChelataseCof.get_or_create(user)
        self.set_param('num_chld_feats', num_feats)

    def create_chelatase_feats(self, user):
        ChelataseFeat.get_or_create_small_subunits(user)
        ChelataseFeat.get_or_create_large_subunits(user)
        ChelataseFeat.get_or_create_chlorophyll_pathway(user)
        ChelataseFeat.get_or_create_b12_pathway(user)

