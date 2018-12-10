# Copyright 2018 by Ivan Antonov. All rights reserved.

import os

from django.conf import settings
from django.test import TestCase


class GtdbTestCase(TestCase):
    """This class provides some convenience methods."""

    def get_full_path_to_test_file(self, fn):
        return os.path.join(settings.BASE_DIR, 'gtdb2/tests/data', fn)


