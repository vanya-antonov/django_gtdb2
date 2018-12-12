# Copyright 2018 by Ivan Antonov. All rights reserved.

from django.core.management.base import BaseCommand, CommandError

from gtdb1.models import GtFs

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.fshift import FShift


class Command(BaseCommand):
    help = """Copy validated chelatase fs-genes from GTDB1 into GTDB2 fshifts
    and download/add any corresponding orgs/seqs if needed."""

#    def add_arguments(self, parser):
#        parser.add_argument(
#            'fn', metavar='ALL_SEQS.gbk',
#            help=".gbk file with all intronless org sequences")

    def handle(self, *args, **options):
        
        user = GeneTackDB().get_default_user()
        fshift = FShift.create_from_gtdb1(user, gtdb1_id)
#        org = Org.create_from_gbk(options['fn'])
#        self.print_success("Successfully created org '%s' with id '%s'" %
#                           (org.name, org.id))

