import sys

from Bio import SeqIO

from django.db import models

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.abstract import AbstractUnit, AbstractParam


class Org(AbstractUnit):
    genus = models.CharField(max_length=255)
    phylum = models.CharField(max_length=255)
    kingdom = models.CharField(max_length=255)

    class Meta:
        db_table = 'orgs'

    @classmethod
    def _create_from_SeqRecord(cls, record, user=None):
        if user is None:
            user = cls.db.get_default_user()
        print(user)
        print(user.name)
        sys.exit()

    @classmethod
    def create_from_gbk(cls, fn, user=None):
        org = None
        for record in SeqIO.parse(fn, "genbank"):
            if org is None:
                # This is the 1st record in the file
                org = cls._create_from_SeqRecord(record, user)
            elif record.annotations['organism'] != org.name:
                raise Exception('')

            gtdb.save_SeqRecord(record, 'fasta')
            gtdb.save_SeqRecord(record, 'genbank')
        org.update_org_params()
        org.make_genetack_model()
        org.make_blastdb()
        gtdb.make_global_blastdb()
        """
                org = Org.objects. ... .first()
                if org is not None:


            # Create new org, if needed
            org_id = Org.get_db_id_by_SeqRecord(gtdb, record)
            if org_id is None:
                org_id = Org.create_new_in_db_from_SeqRecord(gtdb, record, args.user_id)
                org = Org(gtdb, org_id)
                org.add_param('source_fn', fn)
                logging.info("New org '{}' has been created ({})".format(record.annotations['organism'], org_id))
           
            # Create new seq, if needed
            seq_id = Seq.get_db_id_by_ext_id(gtdb, record.id)
            if seq_id is None:
                org = Org(gtdb, org_id)
                seq_id = Seq.create_new_in_db_from_SeqRecord(gtdb, record, org, args.user_id)
                logging.info("New seq '{}' has been created ({})".format(record.id, seq_id))
                """

class OrgParam(AbstractParam):
    org = models.ForeignKey(Org, on_delete=models.CASCADE)

    class Meta:
        db_table = 'org_params'

