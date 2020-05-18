# Copyright 2018 by Ivan Antonov. All rights reserved.

import os
from pprint import pprint
import subprocess

from django.conf import settings

from chelatase_db.lib.bio import run_blastp_seq_vs_file
from chelatase_db.lib.config import PATHWAY_GENES_FAA, read_pathway_gene_info
from gtdb2.models.feat import Feat


class ChelataseFeat(Feat):
    class Meta:
        proxy = True

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(Feat.PRM_INFO.items()) + list({
        'chel_evalue': {'type_fun': float},
        'chel_gene': {},
        'chel_gene_group': {},
        'chel_pathway': {},
        'chel_query': {},
        'chel_subunit': {},
    }.items()))

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

    def create_all_params(self):
        """Extend the parent's method to create additional params."""
        super().create_all_params()
        self._make_params_chel_pathway()

    def _make_params_chel_pathway(self):
        """Tries to determine if this feature (CDS or fsCDS) is homologous
        to some gene from the Chlorophyll or B12 biosynthesis pathway.
        """
        if 'translation' not in self.prm:
            return

        info_dict = read_pathway_gene_info()

        faa_fn = os.path.join(settings.BASE_DIR, PATHWAY_GENES_FAA)
        all_hsps = run_blastp_seq_vs_file(self.prm.translation, faa_fn)

        aa_len = len(self.prm.translation)
        best_hsp = None
        for hsp in all_hsps:
            info = info_dict[hsp.sbjct_id]
            # Make sure the Feat satisfies the requirements
            if '_min_len' in info and aa_len < info['_min_len']:
                continue
            if '_max_len' in info and aa_len > info['_max_len']:
                continue
            best_hsp = hsp
            break

        if best_hsp is None:
            return

        # Set params from info_dict: 'chel_gene', 'chel_subunit', etc
        for name, value in info_dict[best_hsp.sbjct_id].items():
            if value != '' and not name.startswith('_'):
                self.set_param(name, value)

        self.set_param('chel_evalue', '%.2e' % best_hsp.expect,
                       num=best_hsp.expect)


def _get_best_hit(hits_dict, info_dict, aa_len):
    # Make a list of all the HSPs and sort it by evalue
    all_hsps = []
    for hit_id, hsp_list in hits_dict.items():
        for hsp in hsp_list:
            # add missing info
            hsp.hit_id =  hit_id
            all_hsps.append(hsp)
    all_hsps = sorted(all_hsps, key=lambda hsp: hsp.expect)

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

    return best_hsp

