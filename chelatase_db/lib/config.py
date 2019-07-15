# Copyright 2018 by Ivan Antonov. All rights reserved.

import csv
import logging
import os
#import pickle

from django.conf import settings


# For BLAST run
MAX_NUM_THREADS = 6
DEFAULT_EVALUE = 1e-6

PATHWAY_GENES_FAA = 'chelatase_db/data/pathway_genes.faa'
PATHWAY_GENES_TXT = 'chelatase_db/data/pathway_genes.txt'
CHEL_VIEWS_SQL = 'chelatase_db/data/chel_views.sql'
KEGG_ORGS_TXT = 'chelatase_db/data/kegg_orgs.txt'


def read_kegg_orgs():
    """Returns a dict of dicts where keys are KEGG org_names. The duplicated
    org_names are removed.
    """
    # The pickle doesn't exist -- read the text file and create the pickle
    info_fn = os.path.join(settings.BASE_DIR, KEGG_ORGS_TXT)
    info_dict = {}
    duplicated_names = set()
    with open(info_fn) as f:
        logging.info('Reading KEGG data from text file...')
        # https://stackoverflow.com/a/14158869/310453
        lines_wo_comments = filter(lambda row: row[0]!='#', f)
        for d in csv.DictReader(lines_wo_comments, delimiter="\t"):
            if d['org_name'] in info_dict:
                duplicated_names.add(d['org_name'])
            else:
                info_dict[d['org_name']] = d

    # Remove duplidated names
    for name in duplicated_names:
        # https://stackoverflow.com/a/11277439/310453
        info_dict.pop(name)

    return info_dict

def read_pathway_gene_info():
    """Returns a dict of dicts where keys of the first dict are the
    sequence IDs (the _id column of the .txt file).
    """
    # Read the info from file
    info_fn = os.path.join(settings.BASE_DIR, PATHWAY_GENES_TXT)
    info_dict = {}
    with open(info_fn) as f:
        # https://stackoverflow.com/a/14158869/310453
        lines_wo_comments = filter(lambda row: row[0]!='#', f)
        reader = csv.DictReader(lines_wo_comments, delimiter="\t")
        for row in reader:
            if row['_id'] in info_dict:
                raise ValueError(
                    "Sequence ID '%s' is duplicated in file '%s'" %
                    (row['_id'], fn))
            info_dict[row['_id']] = row

    # Postprocess the created dict of dicts
    for row in info_dict.values():
        # Remove any special key that does not have a value
        # About list(row.keys()): https://stackoverflow.com/a/11941855/310453
        for key in list(row.keys()):
            if key.startswith('_') and (row[key] == '' or  row[key] is None):
                del(row[key])

        # Convert some strings to integers, if available
        if '_min_len' in row:
            row['_min_len'] = int(row['_min_len'])
        if '_max_len' in row:
            row['_max_len'] = int(row['_max_len'])

    return info_dict

