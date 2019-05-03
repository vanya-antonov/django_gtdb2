# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint

from chelatase_db.models.feat import ChelataseFeat
from chelatase_db.tests import ChelataseTestCase
from gtdb2.models.org import Org
from gtdb2.models.seq import Seq


class ChelataseFeatModelTests(ChelataseTestCase):

    def setUp(self):
        super().setUp()

        fn = self.get_full_path_to_test_file('N_mexicana.gbk')
        org = Org.create_from_gbk(self.user, fn)
        self.seq = Seq.create_from_ext_id(
            self.user, org, 'NZ_BDBV01000005.1')

    def test_feat_make_all_params_chld(self):
        """Make sure the chel_* params are automatically added to
        a new chlD feat.
        """
        # make_all_params() is called inside the create method
        start, end = 14612-1, 16618   # convert to 0-based coordinates
        feat = ChelataseFeat.get_or_create_from_gbk_annotation(
            self.user, self.seq, start, end, -1)

        # This is the medium chelatase subunit
        self.assertEqual(feat.descr, "magnesium chelatase")

        # ... and it should be automatically annotated
        self.assertEqual(feat.prm.chel_subunit, 'M')
        self.assertEqual(feat.prm.chel_gene_group, 'chlD_bchD')
        self.assertEqual(feat.prm.chel_pathway, 'Chlorophyll')
        self.assertTrue(feat.prm.chel_evalue < 1e-60)

    def test_feat_make_all_params_b12(self):
        """Make sure that some of the chel_* params are automatically added to
        a new feat corresponding to a B12-pathway gene.
        """
        # make_all_params() is called inside the create method
        start, end = 24455-1, 26180   # convert to 0-based coordinates
        feat = ChelataseFeat.get_or_create_from_gbk_annotation(
            self.user, self.seq, start, end, -1)

        # This is one of the B12 pathway genes
        self.assertEqual(feat.descr, "cobyric acid synthase")

        # ... and it should be automatically annotated
        self.assertEqual(feat.prm.chel_gene, 'cobQ')
        self.assertEqual(feat.prm.chel_pathway, 'B12')
        self.assertTrue(feat.prm.chel_evalue < 1e-100)

        # ... but not applicable params should not be added
        self.assertTrue('chel_subunit' not in feat.prm)

