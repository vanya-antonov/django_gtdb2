# Copyright 2018 by Ivan Antonov. All rights reserved.

import csv
import json
import logging
import os

from django.conf import settings

from chelatase_db.lib.config import FEAT_PARAMS_BED
from gtdb2.lib.db import GeneTackDB


class ChelataseDB(GeneTackDB):

    def get_prm_from_bed_for_feat(self, feat):
        """Returns a content of the CHELATASE_SUBUNITS_BED file as dict of dicts.
        """
        # Chache the data obtained from the file
        if not hasattr(self, '_feat_params_bed'):
            self._feat_params_bed = read_feat_params_bed()
        # TODO: use bedtools intersect instead of exact match by key
        return self._feat_params_bed.get(
            feat.prm.location_str, [])


def read_feat_params_bed():
    """Returns a dict of lists where keys are location strings (e.g.
    'NC_00911.1:123-456(+)'). The duplicated keys are not allowed.
    """
    info_fn = os.path.join(settings.BASE_DIR, FEAT_PARAMS_BED)
    info_dict = {}
    with open(info_fn) as f:
        logging.info('Reading data from file "%s"...' % info_fn)
        reader = csv.DictReader(
            f, delimiter="\t",
            fieldnames=['chrom', 'start', 'end', 'name', 'score', 'strand'])
        for d in reader:
            d['start'] = int(d['start'])
            d['end'] = int(d['end'])
            location_str = '%s:%d-%d(%s)' % (
                d['chrom'], d['start'], d['end'], d['strand'])

            # d['name'] is JSON dict of list of dicts
            prm = json.loads(d['name'])
            if type(prm) == dict:
                prm = [prm]
            info_dict.setdefault(location_str, []).extend(prm)

    return info_dict


# Create global shared object for common use
CHEL_DB = ChelataseDB()
