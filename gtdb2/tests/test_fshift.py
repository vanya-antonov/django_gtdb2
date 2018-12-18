# Copyright 2018 by Ivan Antonov. All rights reserved.

import pickle

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.fshift import Fshift, FshiftParam
from gtdb2.models.org import Org
from gtdb2.models.seq import Seq
from gtdb2.tests import GtdbTestCase


class FshiftModelTests(GtdbTestCase):
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
        gtdb = GeneTackDB()
        user = gtdb.get_default_user()

        # All the seq files will be created in a temporary folder
        Org.subdir = gtdb.make_tmp_dir()

        # Prepare required objects
        gbk_fn = self.get_full_path_to_test_file('M_fervens.gbk')
        org = Org.create_from_gbk(user, gbk_fn)
        seq = Seq.create_from_ext_id(user, org, 'NC_013156.1')
        fs_fn = self.get_full_path_to_test_file('gtdb1_fs_33540948.pickle')
        with open(fs_fn, 'rb') as handle:
            gtdb1_fs = pickle.load(handle)

        fshift = Fshift.create_from_gtdb1_fs(user, seq, gtdb1_fs)

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

