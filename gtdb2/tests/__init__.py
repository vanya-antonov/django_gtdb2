# Copyright 2018 by Ivan Antonov. All rights reserved.

import logging
#logging.basicConfig(level=logging.DEBUG)
import os
import shutil
import tempfile

from django.conf import settings
from django.test import TestCase, override_settings

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.org import Org

# All the files (if any) will be created in a tmp dir
# https://docs.djangoproject.com/en/2.1/topics/testing/tools/#overriding-settings
@override_settings(GTDB_DIR=tempfile.mkdtemp(prefix='__gtdb.'))

class GtdbTestCase(TestCase):
    """This class add some commonly used attributes and provides some
    convenience methods.
    """

    test_data_dir = 'gtdb2/tests/data/'

    def setUp(self):
        self.gtdb = GeneTackDB()
        self.user = self.gtdb.get_or_create_default_user()

    def tearDown(self):
        "Deletes tmp dir."
        # Make sure that we delete tmp folder and not the real GTDB_DIR
        if '__gtdb.' in settings.GTDB_DIR:
            shutil.rmtree(settings.GTDB_DIR)
            logging.debug("Tmp folder '%s' has been removed." %
                          settings.GTDB_DIR)
        else:
            logging.warning("GTDB_DIR is NOT a tmp folder '%s' in testing!" %
                            settings.GTDB_DIR)

    def get_full_path_to_test_file(self, fn):
        return os.path.join(settings.BASE_DIR, self.test_data_dir, fn)

