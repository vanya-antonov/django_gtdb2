# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint

from chelatase_db.models.fshift import ChelataseFshift
from chelatase_db.models.org import ChelataseOrg
from chelatase_db.models.seq import ChelataseSeq
from chelatase_db.tests import ChelataseTestCase


class ChelataseOrgModelTests(ChelataseTestCase):

    def test_org_get_or_create_feats_from_cof_by_tblastn(self):
        da_gtdb1_id = 286671156   # da = Delftia acodovorans
        mf_gtdb1_id = 33540948    # mf = Methanocaldococcus fervens

        # Create chld COF with 1 fshift only
        chld_cof = self.create_chld_colf_from_pickles(
            seed_gtdb1_ids = [mf_gtdb1_id])

        # Make sure that org, seq and fshift were created
        self.assertEqual(ChelataseOrg.objects.count(), 1)
        self.assertEqual(ChelataseSeq.objects.count(), 1)
        self.assertEqual(ChelataseFshift.objects.count(), 1)

        mf_org = ChelataseOrg.objects.first()
        self.assertEqual(mf_org.name, 'Methanocaldococcus fervens AG86')

        mf_fshift = ChelataseFshift.objects.first()
        self.assertEqual(mf_fshift.name, 'NC_013156.1:1122957:-1')

        # Create params required for tblastn run
        mf_org._make_param_blastdb()
        mf_org._make_param_transl_table()
        mf_chld_feats = mf_org.get_or_create_feats_from_cof_by_tblastn(
            self.user, chld_cof)

        # Make sure no new objects were created
        self.assertEqual(ChelataseOrg.objects.count(), 1)
        self.assertEqual(ChelataseSeq.objects.count(), 1)
        self.assertEqual(ChelataseFshift.objects.count(), 1)

        # Methanocaldococcus fervens contains 1 chlD gene only
        self.assertEqual(len(mf_chld_feats), 1)

        # This chlD contains -1 frameshift
        mf_chld = mf_chld_feats[0]
        self.assertEqual(mf_chld.fshift.len, -1)


#        org_id = chld_cof.fshift_set.first().org.id
#        org = ChelataseOrg.objects.get(pk=org_id)

        #org.get_or_create_chld_feats(self.user)
#        org.make_all_params(self.user)
#        print("TEST IS FINISHED!")

