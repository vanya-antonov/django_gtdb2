# Copyright 2018 by Ivan Antonov. All rights reserved.

from gtdb2.models.org import Org
from gtdb2.models.seq import Seq
from gtdb2.models.feat import Feat
from gtdb2.tests import GtdbTestCase


class SeqModelTests(GtdbTestCase):

    def setUp(self):
        super().setUp()
        fn = self.get_full_path_to_test_file('N_mexicana.gbk')
        self.org = Org.create_from_gbk(self.user, fn)
        self.seq = Seq.create_from_ext_id(
            self.user, self.org, 'NZ_BDBV01000005.1')

    def test_feat_create_from_gbk_annotation(self):
        feat = Feat.get_or_create_from_gbk_annotation(
            self.user, self.seq, 14612-1, 16618, -1)

        # Check attributes
        self.assertEqual(feat.name, 'NM1_RS02650')
        self.assertEqual(feat.descr, 'magnesium chelatase')
        self.assertEqual(feat.type, 'CDS')
        self.assertEqual(feat.origin, 'annotation')

        # Check params
        self.assertEqual(feat.prm['protein_id'], 'WP_068013468.1')
        self.assertTrue(feat.prm['translation'].startswith('MAPRTDTTTSARD'))

