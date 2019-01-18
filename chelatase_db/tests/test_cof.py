# Copyright 2018 by Ivan Antonov. All rights reserved.

from chelatase_db.models.cof import ChelataseCof
from chelatase_db.tests import ChelataseTestCase


class ChelataseCofModelTests(ChelataseTestCase):

    def test_cof_create_from_gtdb1_seed_fshifts(self):
        cof = self.create_chld_cof_from_pickles()

        self.assertEqual(cof.name, 'chlD')
        self.assertTrue(cof.fshift_set.count() > 0)

