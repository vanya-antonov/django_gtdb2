# Copyright 2018 by Ivan Antonov. All rights reserved.

import pickle

from gtdb2.models.cof import Cof
from gtdb2.models.fshift import Fshift
from gtdb2.models.org import Org
from gtdb2.models.seq import Seq
from gtdb2.tests import GtdbTestCase


class CofModelTests(GtdbTestCase):
    def test_cof_create_from_seed_fshifts(self):
        """See the test_fshift.py for more info on how the gtdb1-fs pickle was
        created.
        """
        # Prepare required objects
        gbk_fn = self.get_full_path_to_test_file('M_fervens.gbk')
        org = Org.create_from_gbk(self.user, gbk_fn)
        seq = Seq.create_from_ext_id(self.user, org, 'NC_013156.1')
        fs_fn = self.get_full_path_to_test_file('gtdb1_fs_33540948.pickle')
        with open(fs_fn, 'rb') as handle:
            gtdb1_fs = pickle.load(handle)
        fshift = Fshift.create_from_gtdb1_fs(self.user, seq, gtdb1_fs)

        seed_fshifts = [fshift]

        # Initially fshifts do not have the seed flag
        for fshift in seed_fshifts:
            self.assertFalse(fshift.seed)

        cof_name = 'Test_cof'
        cof = Cof.create_from_seed_fshifts(
            self.user, seed_fshifts, name=cof_name)

        # The seed flag should be set after fshifts are added to the new cof
        for fshift in seed_fshifts:
            self.assertTrue(fshift.seed)

        # Check the seed_fshifts property
        self.assertEqual(set(cof.seed_fshifts), set(seed_fshifts))

        self.assertEqual(cof.name, cof_name)
        self.assertEqual(cof.fshift_set.count(), 1)
        self.assertEqual(fshift.cof_set.first(), cof)

        # Make sure that I can't create another COF from the same fshift
        with self.assertRaises(ValueError):
            Cof.create_from_seed_fshifts(self.user, [fshift], name=cof_name)

        # Make sure the cof can be properly deleted
        cof.delete()
        self.assertEqual(Cof.objects.count(), 0)

        # Make sure the seed flag was removed from fshifts upon cof deletion
        for fshift in Fshift.objects.all():
            self.assertFalse(fshift.seed)

