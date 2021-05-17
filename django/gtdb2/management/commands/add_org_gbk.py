from django.core.management.base import BaseCommand, CommandError

from gtdb2.models.org import Org

from gtdb2.lib.command import AbstractCommand
from gtdb2.lib.db import GeneTackDB


class Command(AbstractCommand):
    help = "create new org in db and copy all seqs to gtdb dir"

    def add_arguments(self, parser):
        parser.add_argument(
            'fn', metavar='ALL_SEQS.gbk',
            help=".gbk file with all intronless org sequences")

    def handle(self, *args, **options):
        gbk_fn = options['fn']
        gtdb = GeneTackDB()
        user = gtdb.get_or_create_default_user()

        org = Org.get_or_create_from_gbk(user, gbk_fn)
        self.stdout.write("Successfully created org '%s' with id '%s'" %
            (org.name, org.id))

