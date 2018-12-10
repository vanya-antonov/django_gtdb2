# Copyright 2018 by Ivan Antonov. All rights reserved.

import os

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.org import Org
from gtdb2.tests import GtdbTestCase


class OrgModelTests(GtdbTestCase):
    def test_org_create_from_gbk(self):
        gtdb = GeneTackDB()
        user = gtdb.get_default_user()
        fn = self.get_full_path_to_test_file('N_mexicana.gbk')

        # All the seq files will be created in a temporary folder
        #Org.subdir = gtdb.make_tmp_dir(keep=True)
        Org.subdir = gtdb.make_tmp_dir()
        org = Org.create_from_gbk(user, fn)

        # org.name is species name, so it must include genus, e.g.
        # org.genus = 'Nocardia'
        # org.name = 'Nocardia mexicana NBRC 108244'
        self.assertEqual(org.name, 'Nocardia mexicana NBRC 108244')
        self.assertIn(org.genus, org.name)

        # Make sure the fn was saved as one of the params
        self.assertEqual(org.orgparam_set.get(name='source_fn').value, fn)

        # Test AbstractUnit computed properties
        self.assertEqual(org.param_set.get(name='source_fn').value, fn)
        self.assertEqual(org.param_dict['source_fn'][0].value, fn)
        self.assertEqual(org.prm['source_fn'], fn)

        # Make sure org dir was created
        org_dir = gtdb.get_full_path_to(org.prm['dir_path'])
        self.assertTrue(os.path.isdir(org_dir))

        # Make sure the fna and gbk files were created
        seq_id = 'NZ_BDBV01000005.1'
        fna_fn = os.path.join(org_dir, 'seq_fna', seq_id)
        gbk_fn = os.path.join(org_dir, 'seq_gbk', seq_id)
        self.assertTrue(os.path.exists(fna_fn) and os.path.getsize(fna_fn) > 0)
        self.assertTrue(os.path.exists(gbk_fn) and os.path.getsize(gbk_fn) > 0)
        self.assertEqual(org.prm['num_seqs'], 3)

        # Check other generated params
        self.assertEqual(org.prm['short_name'], 'N. mexicana')
        self.assertEqual(org.prm['taxonomy'][0], 'Bacteria')
        self.assertEqual(org.prm['taxonomy'][-1], org.genus)

        # Make sure the org dir is removed if the org is deleted from db
        org.delete()
        self.assertFalse(os.path.isdir(org_dir))

