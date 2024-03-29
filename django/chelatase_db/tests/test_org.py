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

        # Created chlD COF from pickles so that
        # ChelataseCof.get_or_create_chld_cof() will return an existing
        # COF. To speed up the process we create COF with 1 fshift only.
        self.chld_cof = self.create_chld_cof_from_pickles(
            [self.MF_GTDB1_ID])

    def test_org_get_or_create_feats_from_cof_by_tblastn(self):
        # M. fervens org should be created together with the COF
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

        # ...and this is the only fshift in this genome
        self.assertEqual(mf_org.fshift_set.count(), 1)

        # Make sure the parent CDS was created as well
        mf_chli = mf_chld.parent
        self.assertEqual(mf_chli.prm.chel_subunit, 'S')



    def test_org_get_or_create_feats_from_cof_by_tblastn_2(self):
        """Test creation of fshifts located on the minus strand."""
        # Create chld COF with 1 fshift only
        chld_cof = self.create_chld_cof_from_pickles(
            seed_gtdb1_ids = [self.DA_GTDB1_ID])

        # The org was created together with the COF
        org = ChelataseOrg.objects.filter(
            name='Delftia acidovorans SPH-1'
        ).first()

        # Create params required for tblastn run
        org._make_param_blastdb()
        org._make_param_transl_table()
        all_feats = org.get_or_create_feats_from_cof_by_tblastn(
            self.user, chld_cof)

        # There should be 1 chlD gene only
        self.assertEqual(len(all_feats), 1)

        # This chlD gene contains a single -1 frameshift
        feat = all_feats[0]
        self.assertEqual(feat.fshift.len, -1)
        self.assertEqual(feat.prm.chel_subunit, 'M')

        # Make sure the parent CDS was created as well
        self.assertTrue(feat.parent is not None)
        self.assertEqual(feat.parent.prm.chel_subunit, 'S')


    def test_synthesis_classification_empty(self):

        gbk_fn = self.get_full_path_to_test_file('P_aeruginosa.gbk')
        org = ChelataseOrg.create_from_gbk(self.user, gbk_fn)   
        

        pprint(org.prm)
        


    def test_org_create_from_gbk(self):
        """Creates chlD and other pathway feats by tBLASTn."""

        # make_all_params() is called inside the 'create' method
        gbk_fn = self.get_full_path_to_test_file('N_mexicana.gbk')
        org = ChelataseOrg.create_from_gbk(self.user, gbk_fn)

        # The org has only one chlD gene
        self.assertEqual(org.prm.chel_num_chld, 1)
        self.assertEqual(org.chld_feat_set.count(), 1)
        chld_feat = org.chld_feat_set.first()

        self.assertEqual(chld_feat.descr, 'magnesium chelatase')
        self.assertEqual(chld_feat.prm.chel_subunit, 'M')

        # The chlD gene does not have a frameshift
        self.assertEqual(chld_feat.fshift_set.count(), 0)
        self.assertEqual(chld_feat.fshift, None)

        # ... and therefore does not have a parent
        self.assertEqual(chld_feat.parent, None)

        # The chlD gene should be automatically annotated
        self.assertEqual(chld_feat.prm.chel_pathway, 'Chlorophyll')
        self.assertEqual(chld_feat.prm.chel_subunit, 'M')
        self.assertTrue(chld_feat.prm.chel_gene in ['chlD', 'bchD'])
        self.assertTrue(chld_feat.prm.chel_evalue < 1e-60)

        # In addition to the medium subunit, it has the large subunit as well
        large_feat_set = ChelataseFeat.objects.filter(
            seq__org=org,
            param_set__name='chel_subunit',
            param_set__value='L')
        self.assertEqual(large_feat_set.count(), 1)
        cobn_feat = large_feat_set.first()

        # Verify the annotation of the cobN feat
        self.assertEqual(cobn_feat.prm.chel_pathway, 'B12')
        self.assertEqual(cobn_feat.prm.chel_subunit, 'L')
        self.assertEqual(cobn_feat.prm.chel_gene, 'cobN')

        # Verify the genotype string
        self.assertEqual(org.prm.chel_genotype_genes, '1xcobN, 1xchlD')

        # The org also has some other genes from the B12 pathway
        chel_feats = [f for f in org.feat_set.all() if 'chel_gene_group' in f.prm]
        cobQ_feats = [f for f in chel_feats if f.prm.chel_gene_group == 'cobQ']
        self.assertEqual(len(cobQ_feats), 1)

        cobO_feats = [f for f in chel_feats if f.prm.chel_gene_group == 'cobO']
        self.assertEqual(len(cobO_feats), 1)

        cobA_feats = [f for f in chel_feats if f.prm.chel_gene_group == 'cysG_cobA']
        self.assertEqual(len(cobA_feats), 1)

        # Test the get_full_pathway_gene_dict()
        b12_genes = org.get_full_pathway_gene_dict(pathway='B12')
        self.assertTrue('cobN' in b12_genes)
        self.assertTrue('cobP_cobU' in b12_genes)
        self.assertTrue('chlD_bchD' not in b12_genes)
        self.assertEqual(len(b12_genes['cobN']), 1)
        self.assertEqual(len(b12_genes['cobT']), 0)

        chl_genes = org.get_full_pathway_gene_dict(pathway='Chlorophyll')
        self.assertTrue('cobN' not in chl_genes)
        self.assertTrue('cobP_cobU' not in chl_genes)
        self.assertTrue('chlD_bchD' in chl_genes)
        self.assertTrue('chlH_bchH' in chl_genes)
        self.assertEqual(len(chl_genes['chlD_bchD']), 1)
        self.assertEqual(len(chl_genes['chlH_bchH']), 0)

    def test_org_create_from_gbk_witout_fs(self):
        """Make sure that fs-chlD gene is NOT created, because
        this bacteria already has full-length chlD gene and
        frameshifting is not needed (it was a bug to create
        fs-chlD gene).
        """

        # Number of fshifts before we create the new org
        num_fs_before = ChelataseFshift.objects.count()

        gbk_fn = self.get_full_path_to_test_file('R_sulfidophilum.gbk')
        org = ChelataseOrg.create_from_gbk(self.user, gbk_fn)

        # We should not create a new frameshift for this org
        num_fs_after = ChelataseFshift.objects.count()
        self.assertEqual(num_fs_before, num_fs_after)

    def test_org_make_params_kegg(self):
        # M. fervens org should be created together with the COF
        mf_org = ChelataseOrg.objects.filter(
            name='Methanocaldococcus fervens AG86'
        ).first()
        mf_org._make_params_kegg()
        self.assertEqual(mf_org.prm.kegg_org_code, 'mfe')

        # Create chld COF with 1 fshift only
        self.create_chld_cof_from_pickles(
            seed_gtdb1_ids = [self.MS_GTDB1_ID])
        ms_org = ChelataseOrg.objects.filter(
            name='Methanocaldococcus sp. FS406-22'
        ).first()
        ms_org._make_params_kegg()
        self.assertEqual(ms_org.prm.kegg_org_code, 'mfs')

    def test_synthesis_classification(self):
        # P_aeruginosa org should have classfication params
        gbk_fn = self.get_full_path_to_test_file('P_aeruginosa.gbk')
        org = ChelataseOrg.create_from_gbk(self.user, gbk_fn)
         
        print(org.prm['taxonomy']) 
        self.assertEqual(org.prm.chel_synthesis_chl, 0)
        self.assertEqual(org.prm.chel_synthesis_b12, 1)


