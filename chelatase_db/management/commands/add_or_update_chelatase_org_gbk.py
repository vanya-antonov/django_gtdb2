# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint

from Bio import SeqIO

from django.core.management.base import BaseCommand, CommandError

#from chelatase_db.models.cof import ChelataseCof
from chelatase_db.models.org import ChelataseOrg
from gtdb2.lib.db import GeneTackDB


class Command(BaseCommand):
    help = """Creates a new org in db (or updates seqs if the org already
    exists) and generates all the chelatase params. It also creates
    chelatase seed fshifts and cof if needed."""

    def add_arguments(self, parser):
        parser.add_argument(
            'fn', metavar='ORG.gbk',
            help=".gbk file with all intronless org sequences")

    def handle(self, *args, **options):
        gbk_fn = options['fn']
        gtdb = GeneTackDB()
        user = gtdb.get_or_create_default_user()

        org = ChelataseOrg.get_by_gbk(gbk_fn)
        if org is None:
            org = ChelataseOrg.create_from_gbk(user, gbk_fn)
            self.stdout.write(self.style.SUCCESS(
                "The new org '%s' has been created with '%s' chlD genes" %
                (org.name, org.prm['num_chld_feats'])))
        else:
            new_seqs = org.update_seqs_with_gbk(gbk_fn)
            if len(new_seqs) > 0:
                org.make_all_params(user)
                self.stdout.write(self.style.SUCCESS(
                    "'%s' seqs have been created/updated in org '%s'" %
                    (len(new_seqs), org.name)
            else:
                self.stdout.write(self.style.SUCCESS(
                    "All '%s' sequences are up-to-date" % org.name))

