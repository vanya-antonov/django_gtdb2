# Copyright 2018 by Ivan Antonov. All rights reserved.

import pickle
from pprint import pprint

from chelatase_db.models.cof import ChelataseCof
from chelatase_db.tests import ChelataseTestCase
from gtdb2.lib.db import GeneTackDB
from gtdb2.models.org import Org


class ChelataseCofModelTests(ChelataseTestCase):
    def test_cof_create_from_seed_fshifts(self):
        # Prepare required objects
        #all_gtdb1_ids = ChelataseCof.info['seed_fshifts'].keys()
        all_gtdb1_ids = [33540948, 946577251]
        gtdb1_fs_list = []
        for gtdb1_id in all_gtdb1_ids:
            fs_fn = self.get_full_path_to_test_file(
                'gtdb1_fs_' + str(gtdb1_id) + '.pickle')
            with open(fs_fn, 'rb') as handle:
                gtdb1_fs_list.append(pickle.load(handle))

        cof = ChelataseCof.create_from_gtdb1_seed_fshifts(
            self.user, gtdb1_fs_list)

        self.assertEqual(cof.name, ChelataseCof.info['name'])
        self.assertEqual(cof.descr, ChelataseCof.info['descr'])
        self.assertEqual(len(cof.fshift_set.all()), len(all_gtdb1_ids))


