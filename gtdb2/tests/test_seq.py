# Copyright 2018 by Ivan Antonov. All rights reserved.

from datetime import datetime
import os

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.org import Org
from gtdb2.models.seq import Seq
from gtdb2.tests import GtdbTestCase


class SeqModelTests(GtdbTestCase):
    def test_seq_create_from_ext_id(self):
        fn = self.get_full_path_to_test_file('N_mexicana.gbk')
        org = Org.create_from_gbk(self.user, fn)

        seq = Seq.create_from_ext_id(self.user, org, 'NZ_BDBV01000005.1')

        # Check attributes
        self.assertEqual(seq.c_date.year, datetime.now().year)
        self.assertEqual(seq.id, 'NZ_BDBV01000005.1')
        self.assertEqual(seq.name, 'NZ_BDBV01000005')
        self.assertEqual(seq.type, 'DNA')
        self.assertEqual(seq.len, 34847)

        # Check created params
        self.assertTrue('gbk_fn' in seq.prm)

        # We must save relative path to the file
        rel_path = 'orgs/Nocardia_mexicana_NBRC_108244/seq_gbk/NZ_BDBV01000005.1'
        self.assertEqual(seq.prm['gbk_fn'], rel_path)

        # Make sure the file exists
        gbk_fn = self.gtdb.get_full_path_to(seq.prm['gbk_fn'])
        self.assertTrue(os.path.exists(gbk_fn) and os.path.getsize(gbk_fn) > 0)

        # Make sure seq can be read from file
        self.assertTrue(seq.record.seq.startswith('GAAGCGATAGCCAAGCAACTACA'))

