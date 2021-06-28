# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint

from gtdb2.models.feat import Feat


class TtaFeat(Feat):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(Feat.PRM_INFO.items()) + list({
        'tta__start_coord_tta': {'is_list': True, 'value_attr': 'num', 'type_fun': int},
    }.items()))

