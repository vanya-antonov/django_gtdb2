# Copyright 2018 by Ivan Antonov. All rights reserved.

from chelatase_db.models.feat import ChelataseFeat
from chelatase_db.tests import ChelataseTestCase
from gtdb2.models.org import Org
from gtdb2.models.seq import Seq


class ChelataseFeatModelTests(ChelataseTestCase):

    def test_feat_make_all_params(self):
        fn = self.get_full_path_to_test_file('N_mexicana.gbk')
        org = Org.create_from_gbk(self.user, fn)
        seq = Seq.create_from_ext_id(
            self.user, org, 'NZ_BDBV01000005.1')

        # make_all_params() is called inside the create method
        start, end = 14612-1, 16618   # convert to 0-based coordinates
        feat = ChelataseFeat.get_or_create_from_gbk_annotation(
            self.user, seq, start, end, -1)

        # This is the medium chelatase subunit
        self.assertEqual(feat.descr, "magnesium chelatase")

        # ... and it should be automatically annotated
        self.assertEqual(feat.prm.chel_subunit, 'M')
        self.assertEqual(feat.prm.chel_gene, 'bchD')
        self.assertEqual(feat.prm.chel_pathway, 'Chlorophyll')
        self.assertTrue(feat.prm.chel_blastp_evalue < 1e-80 )

