# Copyright 2018 by Ivan Antonov. All rights reserved.

import os
from pprint import pprint
import shutil

from Bio import SeqIO

from gtdb2.models.org import Org, OrgParam
from gtdb2.tests import GtdbTestCase


class OrgModelTests(GtdbTestCase):

    def setUp(self):
        super().setUp()
        self.fn = self.get_full_path_to_test_file('N_mexicana.gbk')

    def test_org_get_all_seq_ids(self):
        org = Org.create_from_gbk(self.user, self.fn)
        all_ids = org.get_all_seq_ids()

        # Make sure all IDs are present
        true_ids = ['NZ_BDBV01000172.1', 'NZ_BDBV01000005.1',
                    'NZ_BDBV01000039.1']
        self.assertEqual(set(all_ids), set(true_ids))

        # Make sure full paths are returned if requested (.gbk files)
        gbk_paths = org.get_all_seq_ids(fullpath=True)
        self.assertTrue(os.path.exists(gbk_paths[0]))
        self.assertTrue(os.path.getsize(gbk_paths[0]) > 0)
        self.assertIsNotNone(next(SeqIO.parse(gbk_paths[0], "genbank")))
        with self.assertRaises(StopIteration):
            next(SeqIO.parse(gbk_paths[0], "fasta"))

        # Make sure full paths are returned if requested (.fasta files)
        fna_paths = org.get_all_seq_ids(seq_dir='seq_fna', fullpath=True)
        self.assertTrue(os.path.exists(fna_paths[0]))
        self.assertTrue(os.path.getsize(fna_paths[0]) > 0)
        self.assertIsNotNone(next(SeqIO.parse(fna_paths[0], "fasta")))
        with self.assertRaises(StopIteration):
            next(SeqIO.parse(fna_paths[0], "genbank"))

    def test_org_get_or_create_from_gbk(self):
        # Org does not exist yet
        self.assertIsNone(Org.get_by_gbk(self.fn))

        org = Org.get_or_create_from_gbk(self.user, self.fn)
        self.assertIsNotNone(org)

        # Make sure the same org is returned and no new orgs were created
        org2 = Org.get_or_create_from_gbk(self.user, self.fn)
        self.assertEqual(org, org2)
        self.assertEqual(Org.objects.count(), 1)

    def test_org_create_from_gbk(self):
        org = Org.create_from_gbk(self.user, self.fn)

        # org.name is species name, so it must include genus, e.g.
        # org.genus = 'Nocardia'
        # org.name = 'Nocardia mexicana NBRC 108244'
        self.assertEqual(org.name, 'Nocardia mexicana NBRC 108244')
        self.assertIn(org.genus, org.name)

        # Make sure org dir was created
        org_dir = org.get_full_path_to_subdir()
        self.assertTrue(os.path.isdir(org_dir))

        # Make sure the fna and gbk files were created
        seq_id = 'NZ_BDBV01000005.1'
        fna_fn = os.path.join(org_dir, 'seq_fna', seq_id)
        gbk_fn = os.path.join(org_dir, 'seq_gbk', seq_id)
        self.assertTrue(os.path.exists(fna_fn) and os.path.getsize(fna_fn) > 0)
        self.assertTrue(os.path.exists(gbk_fn) and os.path.getsize(gbk_fn) > 0)

        self._validate_params_Nocardia_mexicana(org)

        # Make sure the org dir is removed if the org is deleted from db
        org.delete()
        self.assertFalse(os.path.isdir(org_dir))

    def test_org_run_genetack_gm(self):
        gbk_fn = self.get_full_path_to_test_file('S_griseus.50kb.gbk')
        org = Org.create_from_gbk(self.user, gbk_fn)

        # Use pre-calculated models save execution time
        gm_mod_fn = self.get_full_path_to_test_file('S_griseus.gm_mod.txt')
        fs_mod_fn = self.get_full_path_to_test_file('S_griseus.fs_mod.txt')
        gtgm_out_fn = org.run_genetack_gm(
            gm_mod_fn=gm_mod_fn, fs_mod_fn=fs_mod_fn)

        # Make sure the output file exists
        self.assertTrue(os.path.exists(gtgm_out_fn))

    def test_org_create_genetack_fshifts_from_file(self):
        gbk_fn = self.get_full_path_to_test_file('S_griseus.50kb.gbk')
        org = Org.create_from_gbk(self.user, gbk_fn)

        gtgm_fn = self.get_full_path_to_test_file('S_griseus.50kb.genetackgm')
        fsgenes_fna_fn = self.get_full_path_to_test_file(
            'S_griseus.50kb.fsgene_seqs.fna')
        all_fshifts = org.create_fshifts_from_genetackgm_file(
            self.user, gtgm_fn, fsgenes_fna_fn)

    def test_org_make_all_params_2(self):
        """Make sure all params are preserved and no duplicates are created
        after another function call.
        """
        org = Org.create_from_gbk(self.user, self.fn)

        num_params_1 = OrgParam.objects.count()
        org.create_all_params()
        num_params_2 = OrgParam.objects.count()

        self.assertEqual(num_params_1, num_params_2)

        self._validate_params_Nocardia_mexicana(org)

    def _validate_params_Nocardia_mexicana(self, org):
        # Test the computed attributes & Make sure the fn was saved
        self.assertEqual(org.param_set.get(name='source_fn').value, self.fn)
        self.assertEqual(org.prm_dict_of_lists['source_fn'][0].value, self.fn)
        self.assertEqual(org.prm_dict['source_fn'], self.fn)
        self.assertEqual(org.prm['source_fn'], self.fn)
        self.assertEqual(org.prm.source_fn, self.fn)   # AttrDict style

        # Check the param values
        self.assertEqual(org.prm['num_seqs'], 3)
        self.assertEqual(org.prm['short_name'], 'N. mexicana')
        self.assertEqual(org.prm['transl_table'], [11])
        self.assertEqual(org.transl_table, 11)

        true_taxonomy = ['Bacteria', 'Actinobacteria', 'Corynebacteriales',
                         'Nocardiaceae', 'Nocardia']
        self.assertEqual(org.prm['taxonomy'], true_taxonomy)

        # Make sure external IDs were added
        self.assertEqual(org, OrgParam.get_parent_by_xref(
            'BioProject', 'PRJNA224116'))
        self.assertEqual(org, OrgParam.get_parent_by_xref(
            'BioSample', 'SAMD00018777'))
        self.assertEqual(org, OrgParam.get_parent_by_xref(
            'Assembly', 'GCF_001613165.1'))
        self.assertEqual(org, OrgParam.get_parent_by_xref(
            'taxon', 1210089))

        # Verify the created features - one 16S rRNA
        all_rrna_feats = org.feat_set.filter(type='rRNA').all()
        self.assertEqual(len(all_rrna_feats), 1)
        rrna_feat = all_rrna_feats[0]
        self.assertEqual(rrna_feat.name, 'NM1_RS40620')
        self.assertEqual(rrna_feat.descr, '16S ribosomal RNA')

        # Make sure the seq was loaded
        self.assertTrue(rrna_feat.prm.seq_nt.startswith('TTCAACGGAGAGTT'))
        self.assertTrue(rrna_feat.prm.seq_nt.endswith('GGATCACCTCCTTTCTAA'))

        # This rRNA should also be assigned as org rRNA
        self.assertTrue('seq_rrna_16s' in org.prm)
        self.assertEqual(org.prm.seq_rrna_16s, rrna_feat.prm.seq_nt)

        # Make sure blastdbs were created
        db_path = self.gtdb.get_full_path_to(org.prm['blastdb_nucl_all'])
        db_files = [db_path + '.nsq', db_path + '.nhr', db_path + '.nin']
        for fn in db_files:
            self.assertTrue(os.path.exists(fn) and os.path.getsize(fn) > 0)

