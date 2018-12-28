# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint
import pickle

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.fshift import Fshift, FshiftParam
from gtdb2.models.org import Org
from gtdb2.models.seq import Seq
from gtdb2.tests import GtdbTestCase


class FshiftModelTests(GtdbTestCase):

    def setUp(self):
        super().setUp()

        # Prepare required objects
        gbk_fn = self.get_full_path_to_test_file('M_fervens.gbk')
        org = Org.create_from_gbk(self.user, gbk_fn)
        self.seq = Seq.create_from_ext_id(self.user, org, 'NC_013156.1')

    def test_fshift_simple_create(self):
        fshift = Fshift(user=self.user, seq=self.seq, type='genetack',
                        start=1121859, coord=1122957, end=1123904,
                        len=-1, strand=1)
        fshift.save()

        # Make sure the fshift is present in DB
        self.assertEqual(Fshift.objects.count(), 1)
        self.assertEqual(Fshift.objects.first(), fshift)

    def test_fshift_make_all_params(self):
        fshift = Fshift(user=self.user, seq=self.seq, type='genetack',
                        start=1121859, coord=1122957, end=1123904,
                        len=-1, strand=1)
        fshift.save()

        # name is generated inside make_all_params()
        self.assertIsNone(fshift.name)

        fshift.make_all_params()

        self.assertEqual(fshift.name, 'NC_013156.1:1122957:-1')

        self._validate_params_1122957(fshift)

    def test_fshift_make_all_params_2(self):
        """Make sure all params are preserved and no duplicates are created
        after another function call.
        """
        fshift = Fshift(user=self.user, seq=self.seq, type='genetack',
                        start=1121859, coord=1122957, end=1123904,
                        len=-1, strand=1)
        fshift.save()
        fshift.make_all_params()

        num_params_1 = FshiftParam.objects.count()
        fshift.make_all_params()
        num_params_2 = FshiftParam.objects.count()

        self.assertEqual(num_params_1, num_params_2)

        self._validate_params_1122957(fshift)

    def test_fshift_create_from_gtdb1(self):
        """While testing I can't get access to the data in any production
        database. So, the gtdb1_fs object was stored in a pickle file that
        was created as follows:
        cd ~/_my/github/django_gtdb2
        python3 manage.py shell
        >>> import pickle
        >>> from gtdb1.models import GtFs
        >>> gtdb1_id = 33540948
        >>> gtdb1_fs = GtFs.objects.using('gtdb1').get(pk=gtdb1_id)
        >>> job = gtdb1_fs.job
        >>> out_fn = 'gtdb2/tests/data/gtdb1_fs_%s.pickle' % gtdb1_id
        >>> with open(out_fn, 'wb') as handle:
        >>>     pickle.dump(gtdb1_fs, handle, protocol=pickle.HIGHEST_PROTOCOL)
        """
        gtdb1_fs_id = 33540948
        pickle_fn = 'gtdb1_fs_' + str(gtdb1_fs_id) + '.pickle'
        pickle_path = self.get_full_path_to_test_file(pickle_fn)
        with open(pickle_path, 'rb') as handle:
            gtdb1_fs = pickle.load(handle)
        fshift = Fshift.create_from_gtdb1_fs(
            self.user, self.seq, gtdb1_fs)

        # Make sure the old date was preserved
        self.assertTrue(fshift.c_date.year == 2010)

        # Validate fs attributes
        self.assertEqual(fshift.name, 'NC_013156.1:1122957:-1')
        self.assertEqual(fshift.type, 'genetack')
        self.assertEqual(fshift.coord, 1122957)
        self.assertEqual(fshift.len, -1)
        self.assertEqual(fshift.strand, 1)

        # Make sure the gtdb1 ID was added to xrefs
        fshift_by_ext_id = FshiftParam.get_parent_by_xref('gtdb1', 33540948)
        self.assertIsNotNone(fshift_by_ext_id)
        self.assertEqual(fshift, fshift_by_ext_id)

        self._validate_params_1122957(fshift)

    def _validate_params_1122957(self, fshift):
        "Validates params of fshift with coord 1122957."

        # Seq part before fshift is lowercase, and after -- is uppercase
        self.assertTrue(fshift.prm['seq_nt'].startswith('atg'))
        self.assertFalse(fshift.prm['seq_nt'].startswith('ATG'))
        self.assertTrue(fshift.prm['seq_nt'].endswith('TAG'))
        self.assertFalse(fshift.prm['seq_nt'].endswith('tag'))

        # Check seq_prot_n and seq_prot_c
        seq_prot = fshift.prm['seq_prot'].upper()
        self.assertTrue(seq_prot.startswith(fshift.prm['seq_prot_n']))
        self.assertTrue(seq_prot.endswith(fshift.prm['seq_prot_c']))

        # Corrected DNA seq does not have inner stop codons
        nt_corr = fshift.get_prm_bio_seq('seq_nt_corr')
        prot_corr = nt_corr.translate(table=fshift.seq.transl_table)
        self.assertTrue(prot_corr.endswith('*'))
        prot_corr = prot_corr.rstrip('*')
        self.assertFalse('*' in prot_corr)

        # Translation of the original (frameshifted) seq produces stop codons
        nt_intact = fshift.get_prm_bio_seq('seq_nt')
        prot_intact = nt_intact.translate(table=fshift.seq.transl_table)
        prot_intact = prot_intact.rstrip('*')
        self.assertTrue('*' in prot_intact)

