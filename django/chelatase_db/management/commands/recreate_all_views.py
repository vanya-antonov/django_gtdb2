# Copyright 2018 by Ivan Antonov. All rights reserved.

import os
from pprint import pprint
import subprocess

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import connection

from chelatase_db.lib.config import CHEL_VIEWS_SQL


class Command(BaseCommand):
    help = """Re-create all the views (i.e. tables with the _v suffix)."""

    def add_arguments(self, parser):
        parser.add_argument('run', metavar='__RUN__')

    def handle(self, *args, **options):
        if options['run'] != '__RUN__':
            return

        sql_fn = os.path.join(settings.BASE_DIR, CHEL_VIEWS_SQL)

        # https://docs.djangoproject.com/en/2.2/topics/db/sql/#executing-custom-sql-directly
        with connection.cursor() as cursor:
            query_list = []
            with open(sql_fn) as f:
                for line in f:
                    line = line.strip()
                    if line != '':
                        query_list.append(line)

                    if line.endswith(';'):
                        query_str = ' '.join(query_list)
                        cursor.execute(query_str)
                        query_list = []

        self.stdout.write(self.style.SUCCESS('ALL DONE!'))

