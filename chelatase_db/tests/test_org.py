# Copyright 2018 by Ivan Antonov. All rights reserved.

import os
from pprint import pprint

from django.conf import settings

from chelatase_db.models.cof import ChelataseCof
from chelatase_db.models.feat import ChelataseFeat
from chelatase_db.models.fshift import ChelataseFshift
from chelatase_db.models.org import ChelataseOrg
from chelatase_db.models.seq import ChelataseSeq
from chelatase_db.tests import ChelataseTestCase


class ChelataseOrgModelTests(ChelataseTestCase):

    def setUp(self):
        super().setUp()

        mf_gtdb1_id = 33540948    # mf = Methanocaldococcus fervens

        # Create chld COF with 1 fshift only
        self.chld_cof = self.create_chld_cof_from_pickles(
            seed_gtdb1_ids = [mf_gtdb1_id])

        # Make sure that org, seq and fshift were created
        self.assertEqual(ChelataseOrg.objects.count(), 1)
        self.assertEqual(ChelataseSeq.objects.count(), 1)
        self.assertEqual(ChelataseFshift.objects.count(), 1)

    def test_org_get_or_create_feats_from_cof_by_tblastn(self):
        # M. fervens org should be created in setUp()
        mf_org = ChelataseOrg.objects.filter(
            name='Methanocaldococcus fervens AG86'
        ).first()

        # Create params required for tblastn run
        mf_org._make_param_blastdb()
        mf_org._make_param_transl_table()
        mf_chld_feats = mf_org.get_or_create_feats_from_cof_by_tblastn(
            self.user, self.chld_cof)

        # Methanocaldococcus fervens contains 1 chlD gene only
        self.assertEqual(len(mf_chld_feats), 1)

        # This chlD gene contains a single -1 frameshift
        mf_chld = mf_chld_feats[0]
        self.assertEqual(mf_chld.fshift.len, -1)
        self.assertEqual(mf_chld.fshift.name, 'NC_013156.1:1122957:-1')
        self.assertEqual(mf_chld.prm.chel_subunit, 'M')

        # Make sure the parent CDS was created as well
        mf_chli = mf_chld.parent
        self.assertEqual(mf_chli.prm.chel_subunit, 'S')

    def test_org_make_all_params(self):
        """Creates chlD and other pathway feats by tBLASTn."""

        # make_all_params() is called inside the 'create' method
        gbk_fn = self.get_full_path_to_test_file('N_mexicana.gbk')
        org = ChelataseOrg.create_from_gbk(self.user, gbk_fn)

        # The org has only one chlD gene
        self.assertEqual(len(org.chld_feats), 1)
        chld_feat = org.chld_feats[0]

        self.assertEqual(chld_feat.descr, 'magnesium chelatase')
        self.assertEqual(chld_feat.prm.chel_subunit, 'M')

        # The chlD gene does not have a frameshift
        self.assertEqual(chld_feat.fshift_set.count(), 0)
        self.assertEqual(chld_feat.fshift, None)

        # ... and therefore does not have a parent
        self.assertEqual(chld_feat.parent, None)

        # The chlD gene is more similar to bchD
        self.assertEqual(chld_feat.prm.chel_pathway, 'Chlorophyll')
        self.assertEqual(chld_feat.prm.chel_subunit, 'M')
        self.assertEqual(chld_feat.prm.chel_gene, 'bchD')
        self.assertTrue(chld_feat.prm.chel_evalue < 1e-80)

        # In addition to the medium subunit, it has the large subunit as well
        all_feats = list(ChelataseFeat.objects.filter(seq__org=org).all())
        all_large = [f for f in all_feats if f.prm.chel_subunit == 'L']
        self.assertEqual(len(all_large), 1)
        cobn_feat = all_large[0]

        # Verify the annotation of the cobN feat
        self.assertEqual(cobn_feat.prm.chel_pathway, 'Cobalamin')
        self.assertEqual(cobn_feat.prm.chel_subunit, 'L')
        self.assertEqual(cobn_feat.prm.chel_gene, 'cobN')

