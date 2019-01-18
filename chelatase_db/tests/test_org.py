# Copyright 2018 by Ivan Antonov. All rights reserved.

import os
from pprint import pprint

from django.conf import settings

from chelatase_db.models.cof import ChelataseCof as CCof
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

    def test_org_make_all_params(self):
        "Creates chlD fsCDS in D.acidovorans genome by tBLASTn."
        da_gtdb1_id = 286671156   # da = Delftia acodovorans

        # make_all_params() is called inside the 'create' method
        info = CCof.chld_info['seed_fshifts'][da_gtdb1_id]
        gbk_fn = os.path.join(settings.BASE_DIR, info['org_gbk'])
        org = ChelataseOrg.create_from_gbk(self.user, gbk_fn)

        #len(org.chld_feats)

        pprint(vars(org))


#        org_id = chld_cof.fshift_set.first().org.id
#        org = ChelataseOrg.objects.get(pk=org_id)

        #org.get_or_create_chld_feats(self.user)
#        org.make_all_params(self.user)
#        print("TEST IS FINISHED!")

