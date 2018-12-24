# Copyright 2018 by Ivan Antonov. All rights reserved.

from datetime import datetime
import os

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.org import Org
from gtdb2.models.seq import Seq
from gtdb2.tests import GtdbTestCase


class SeqModelTests(GtdbTestCase):

    def setUp(self):
        super().setUp()
        fn = self.get_full_path_to_test_file('N_mexicana.gbk')
        self.org = Org.create_from_gbk(self.user, fn)

    def test_seq_create_from_ext_id(self):
        seq = Seq.create_from_ext_id(self.user, self.org, 'NZ_BDBV01000005.1')

        # Check attributes
        self.assertEqual(seq.c_date.year, datetime.now().year)
        self.assertEqual(seq.id, 'NZ_BDBV01000005.1')
        self.assertEqual(seq.name, 'NZ_BDBV01000005')
        self.assertEqual(seq.type, 'DNA')
        self.assertEqual(seq.len, 34847)

        # Make sure seq can be read from file
        self.assertTrue(seq.seq.startswith('GAAGCGATAGCCAAGCAACTACA'))
        self.assertTrue(seq.record.seq.startswith('GAAGCGATAGCCAAGCAACTACA'))

    def test_seq_make_all_params(self):
        seq = Seq.create_from_ext_id(self.user, self.org, 'NZ_BDBV01000005.1')
        seq.make_all_params()

        # Check 'transl_table' prm
        self.assertTrue('transl_table' in seq.prm)
        self.assertEqual(seq.prm['transl_table'], 11)

