# Copyright 2018 by Ivan Antonov. All rights reserved.

import os

from django.conf import settings
from django.test import TestCase

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.org import Org


class GtdbTestCase(TestCase):
    """This class add some commonly used attributes and provides some
    convenience methods.
    """

    test_data_dir = 'gtdb2/tests/data/'

    def setUp(self):
        self.gtdb = GeneTackDB()
        self.user = self.gtdb.get_or_create_default_user()

        # All the seq files (if any) will be created in a temporary folder
        #Org.subdir = gtdb.make_tmp_dir(keep=True)
        Org.subdir = self.gtdb.make_tmp_dir()

    def get_full_path_to_test_file(self, fn):
        return os.path.join(settings.BASE_DIR, self.test_data_dir, fn)

