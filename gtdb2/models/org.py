import os
import sys
import re

from Bio import SeqIO

from ivanya.biopython.genbank import get_genus_phylum_kingdom

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
    def create_from_gbk(cls, user, fn):
        org = None
        for record in SeqIO.parse(fn, "genbank"):
            if org is None:
                # This is the 1st record in the file
                org = cls._create_from_SeqRecord(user, record)
                org.add_param('source_fn', fn)
                return org
            elif record.annotations['organism'] != org.name:
                # Make sure that all seqs corresponds to the same org
                raise Exception("Seq '%s' does NOT belong to '%s' in file '%s'" %
                                (record.id, org.name, fn))
            #gtdb.save_SeqRecord(record, 'fasta')
            #gtdb.save_SeqRecord(record, 'genbank')
        org.update_org_params()
        #org.make_genetack_model()
        #org.make_blastdb()
        #gtdb.make_global_blastdb()
        return org

    @classmethod
    def _create_from_SeqRecord(cls, user, record):
        genus, phylum, kingdom = get_genus_phylum_kingdom(record)
        org = Org(user=user,
                  name=record.annotations['organism'],
                  genus=genus,
                  phylum=phylum,
                  kingdom=kingdom)
        org.save()
        return org
        org_dir = cls._create_org_dir(name)

        db_id = gtdb.get_random_db_id('orgs', 'id')
        gtdb.exec_sql_in(
            "INSERT INTO orgs (id, user_id, name, genus, phylum, kingdom, dir_path)",
            db_id, user_id, name, genus, phylum, kingdom, dir_path)
        org = Org(gtdb, db_id)
        
        # ORG_PARAMS
        #dbx_dict = db_xref_list_to_dict(record.dbxrefs)
        #for (name,value) in dbx_dict.items():
        #    org.add_param(name, value)
        
        #num = 0
        #for value in record.annotations.get('taxonomy',[]):
        #    org.add_param('taxonomy', value, num)
        #    num += 1
        
        return db_id
        print(user)
        print(user.name)
        exit()

    @classmethod
    def _create_org_dir(cls, org_name):
        """The function returns new path RELATIVE to the GTDB root dir."""
        # 'Natranaerobius thermophilus JW/NM-WN-LF'  =>  'Natranaerobius_thermophilus_JW_NM_WN_LF'
        folder = re.compile('[^\w\d]+').sub('_', org_name).strip('_')
        subdir = os.path.join(cls.gtdb.root_dir, 'orgs')
        print(folder)
        exit()
        full_path = make_new_fullpath_for_basename(subdir, folder)
        os.makedirs(full_path)
        return os.path.relpath(full_path, gtdb.gtdb_dir)


class OrgParam(AbstractParam):
    parent = models.ForeignKey(Org, on_delete=models.CASCADE)

    class Meta:
        db_table = 'org_params'

