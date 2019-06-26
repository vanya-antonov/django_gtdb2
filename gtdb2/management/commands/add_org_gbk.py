from django.core.management.base import BaseCommand, CommandError

from gtdb2.models.org import Org

from gtdb2.lib.command import AbstractCommand


class Command(AbstractCommand):
    help = "create new org in db and copy all seqs to gtdb dir"

    def add_arguments(self, parser):
        parser.add_argument(
            'fn', metavar='ALL_SEQS.gbk',
            help=".gbk file with all intronless org sequences")

    def handle(self, *args, **options):
        org = Org.create_from_gbk(options['fn'])
        self.print_success("Successfully created org '%s' with id '%s'" %
                           (org.name, org.id))

