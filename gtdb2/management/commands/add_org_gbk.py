from django.core.management.base import BaseCommand, CommandError

from gtdb2.lib.manager import add_org_gbk


class Command(BaseCommand):
    help = "create new org in db and copy all seqs to gtdb dir"

    def add_arguments(self, parser):
        parser.add_argument(
            'fn', metavar='ALL_SEQS.gbk',
            help=".gbk file with all intronless org sequences")

    def handle(self, *args, **options):
        org = add_org_gbk(options['fn'])
        self.stdout.write(self.style.SUCCESS(
            "Successfully created org '%s' with id '%s'" %
            (org.name, org.id)))

