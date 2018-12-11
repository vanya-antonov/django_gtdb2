# Copyright 2018 by Ivan Antonov. All rights reserved.

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.org import Org
from gtdb2.models.seq import Seq
from gtdb2.tests import GtdbTestCase


class SeqModelTests(GtdbTestCase):
    def test_seq_create_from_ext_id(self):
        gtdb = GeneTackDB()
        user = gtdb.get_default_user()
        fn = self.get_full_path_to_test_file('N_mexicana.gbk')

        # All the seq files will be created in a temporary folder
        Org.subdir = gtdb.make_tmp_dir()
        org = Org.create_from_gbk(user, fn)

        seq = Seq.create_from_ext_id(user, org, 'NZ_BDBV01000005.1')

        # Check created params
        self.assertTrue('gbk_fn' in seq.prm)

        self.assertEqual(seq.id, 'NZ_BDBV01000005.1')
        self.assertEqual(seq.name, 'NZ_BDBV01000005')
        self.assertEqual(seq.type, 'DNA')
        self.assertEqual(seq.len, 34847)

        # Make sure seq can be read from file
        self.assertTrue(seq.record.seq.startswith('GAAGCGATAGCCAAGCAACTACA'))

