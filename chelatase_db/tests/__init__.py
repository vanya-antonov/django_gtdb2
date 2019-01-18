# Copyright 2018 by Ivan Antonov. All rights reserved.

import pickle

from chelatase_db.models.cof import ChelataseCof
from gtdb2.tests import GtdbTestCase


class ChelataseTestCase(GtdbTestCase):
    """This class provides some convenience methods."""

    # Overwrite dir path
    test_data_dir = 'chelatase_db/tests/data/'

    def create_chld_cof_from_pickles(self, seed_gtdb1_ids=None):
        """Creates chlD COF in the test database from the GTDB1 frameshifts
        stored in pickle files. Thus, any future calls of
        get_or_create_chld_cof() will return the existing cof.
        """
        if seed_gtdb1_ids is None:
            # Prepare a list of gtdb1 fs
            #seed_gtdb1_ids = ChelataseCof.chld_info['seed_fshifts'].keys()
            #seed_gtdb1_ids = [33540948, 946577251]
            seed_gtdb1_ids = [33540948]

        gtdb1_fs_list = []
        for gtdb1_id in seed_gtdb1_ids:
            fs_fn = self.get_full_path_to_test_file(
                'gtdb1_fs_' + str(gtdb1_id) + '.pickle')
            with open(fs_fn, 'rb') as handle:
                gtdb1_fs_list.append(pickle.load(handle))

        return ChelataseCof.create_from_gtdb1_seed_fshifts(
            self.user, gtdb1_fs_list)

