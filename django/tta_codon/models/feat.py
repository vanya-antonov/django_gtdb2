# Copyright 2018 by Ivan Antonov. All rights reserved.

import os
from pprint import pprint
import subprocess

#from django.conf import settings

#from chelatase_db.lib.bio import run_blastp_seq_vs_file
#from chelatase_db.lib.config import PATHWAY_GENES_FAA, read_pathway_gene_info
from gtdb2.models.feat import Feat


class TtaFeat(Feat):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(Feat.PRM_INFO.items()) + list({
        'tta__start_coord_tta': {'value_attr': 'num', 'type_fun': int},
    }.items()))

