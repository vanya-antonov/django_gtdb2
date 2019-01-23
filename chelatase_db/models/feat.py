# Copyright 2018 by Ivan Antonov. All rights reserved.

import csv
import os
from pprint import pprint
import subprocess

from django.conf import settings

from chelatase_db.lib.bio import run_blastp_vs_file
from gtdb2.models.feat import Feat


class ChelataseFeat(Feat):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(Feat.PRM_INFO.items()) + list({
        'chel_evalue': {'type_fun': float},
        'chel_gene': {},
        'chel_pathway': {},
        'chel_query': {},
        'chel_subunit': {},
    }.items()))


# TODO: Annotation of other pathway genes
# 
# Data file are in q_seq/group_seq and _q_groups.2_col located in
# ~/_my/DataLog/2018_yulia/0802.Ba.heatmap_with_several_queries
# 
# Additional info is in 0802.xlsx located on my Mac in 
# /Users/antonov/Projects/2018/Baranov/DataLog/0802.Ba.heatmap_with_several_queries
# 
# No code was written for this part of the work.
    PATHWAY_GENES_FAA = 'chelatase_db/data/pathway_genes.faa'
    PATHWAY_GENES_TXT = 'chelatase_db/data/pathway_genes.txt'

    @property
    def fshift(self):
        "Simulate one-to-one relationship between chelatase objects."
        num = self.fshift_set.count()
        if num == 0:
            return None
        elif num == 1:
            return self.fshift_set.first()
        else:
            raise ValueError("ChelataseFeat '%s' has '%s' fshifts" %
                             (self, num))

    def make_all_params(self):
        """Extend the parent's method to create additional params."""
        super().make_all_params()
        self._make_params_chel_pathway()

    def _make_params_chel_pathway(self):
        """Tries to determine if this feature (CDS or fsCDS) is homologous
        to some gene from the Chlorophyll or B12 biosynthesis pathway.
        """
        if 'translation' not in self.prm:
            return

        info_fn = os.path.join(settings.BASE_DIR, self.PATHWAY_GENES_TXT)
        info_dict = _read_pathway_gene_info(info_fn)

        faa_fn = os.path.join(settings.BASE_DIR, self.PATHWAY_GENES_FAA)
        hits_dict = run_blastp_vs_file(self.prm.translation, faa_fn)
        if len(hits_dict.keys()) == 0:
            return

        # Make a list of all the HSPs and sort it by evalue
        all_hsps = []
        for hit_id, hsp_list in hits_dict.items():
            for hsp in hsp_list:
                # add missing info
                hsp.hit_id =  hit_id
                all_hsps.append(hsp)
        all_hsps = sorted(all_hsps, key=lambda hsp: hsp.expect)

        aa_len = len(self.prm.translation)
        best_hsp = None
        for hsp in all_hsps:
            info = info_dict[hsp.hit_id]
            # Make sure the Feat satisfies the requirements
            if '_min_len' in info and aa_len < info['_min_len']:
                continue
            if '_max_len' in info and aa_len > info['_max_len']:
                continue
            best_hsp = hsp
            break

        if best_hsp is None:
            return

        self.set_param('chel_evalue', '%.2e' % best_hsp.expect,
                       num=best_hsp.expect)

        # Set params from info_dict: 'chel_gene', 'chel_subunit', etc
        for name, value in info_dict[best_hsp.hit_id].items():
            if not name.startswith('_'):
                self.set_param(name, value)


def _read_pathway_gene_info(fn):
    info_dict = {}
    with open(fn) as f:
        # https://stackoverflow.com/a/14158869/310453
        lines_wo_comments = filter(lambda row: row[0]!='#', f)
        reader = csv.DictReader(lines_wo_comments, delimiter="\t")
        for row in reader:
            if row['_id'] in info_dict:
                raise ValueError(
                    "Sequence ID '%s' is duplicated in file '%s'" %
                    (row['_id'], fn))
            if row['_min_len'].isdigit():
                row['_min_len'] = int(row['_min_len'])
            if row['_max_len'].isdigit():
                row['_max_len'] = int(row['_max_len'])
            info_dict[row['_id']] = row
    return info_dict

