# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint

from Bio import SeqIO

from django.core.management.base import BaseCommand, CommandError

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

        org = ChelataseOrg.get_or_create_from_gbk(user, gbk_fn)
        upd_info = org.update_seqs_with_gbk(gbk_fn)
        if(upd_info['n_new'] == 0 and upd_info['n_updated'] == 0 and
           'num_chld_feats' in org.prm):
            self.stdout.write(self.style.SUCCESS(
                "The '%s' genome is up-to-date (it has '%s' chlD genes)" %
                (org.name, org.prm['num_chld_feats'])))
            return

        org.create_chld_fshifts_and_feats(user)

        pprint(vars(org))
        exit()
        if org.prm['num_chld_feats'] == 0:
            self.stdout.write(
                "The '%s' genome does not have chlD genes!" % org.name)
            return
        org.create_chelatase_feats(user)
        org.make_all_params()

