# Copyright 2018 by Ivan Antonov. All rights reserved.

import csv
import json
import logging
import os

from django.conf import settings

from chelatase_db.lib.config import FEAT_PARAMS_BED
from gtdb2.lib.db import GeneTackDB


class ChelataseDB(GeneTackDB):
    pass


# Create global shared object for common use
CHEL_DB = ChelataseDB()
