# Copyright 2018 by Ivan Antonov. All rights reserved.

import logging
import os
import shutil
import sys
import re

from Bio import SeqIO

from django.db import models
from django.db.models import signals
from django.dispatch import receiver

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.abstract import AbstractUnit, AbstractParam


class Org(AbstractUnit):
    genus = models.CharField(max_length=255)
    phylum = models.CharField(max_length=255)
    kingdom = models.CharField(max_length=255)

    class Meta:
        db_table = 'orgs'

    subdir = 'orgs'
    prm_info = {
        'dir_path': {'prm_attr': 'value'},
        'num_seqs': {'prm_attr': 'num', 'type_fun': int},
        'short_name': {'prm_attr': 'value'},
        'source_fn': {'prm_attr': 'value'},
        'taxonomy': {'prm_attr': 'value', 'is_list': True, 'sort_attr': 'num'},
        'transl_table': {'prm_attr': 'value', 'type_fun': int, 'is_list': True,
                         'sort_attr': 'num', 'reverse': True},
    }

    @classmethod
    def create_from_gbk(cls, user, fn):
        """Creates a new org in db from a genbank file containing all org
        sequences. All the sequences are saved to GTDB dir, but only sequences
        with important annotated features (such as rRNAs or annotated
        programmed frameshifts) are added to db - see make_all_params().
        """
        org = None
        for record in SeqIO.parse(fn, "genbank"):
            if org is None:
                # This is the first record in the file
                org = cls._create_from_SeqRecord(user, record)
                org.add_param('source_fn', fn)
            elif record.annotations['organism'] != org.name:
                # Make sure that all seqs correspond to the same org
                raise Exception("Seq '%s' does NOT belong to '%s' in file '%s'" %
                                (record.id, org.name, fn))
            org._create_seq_file(record, 'fasta-2line')
            org._create_seq_file(record, 'genbank')
        org.make_all_params()
        #org.make_genetack_model()
        #org.make_blastdb()
        #gtdb.make_global_blastdb()
        return org

    @classmethod
    def _create_from_SeqRecord(cls, user, record):
        "Creates a new org in db, creates org dir and adds some basic params."
        sp, ge, ph, ki = _get_species_genus_phylum_kingdom(record)
        if None in (sp, ge, ph, ki):
            raise Exception("Can't determine species, genus, phylum or "
                            "kingdom: %s, %s, %s, %s" % (sp, ge, ph, ki))
        org = Org(user=user, name=sp, genus=ge, phylum=ph, kingdom=ki)
        org.save()

        org._create_org_dir()

        return org

    def _create_org_dir(self):
        "Creates dir for new org and saves its path in the params table."
        # 'Natranaerobius thermophilus JW/NM-WN-LF'  =>
        # 'Natranaerobius_thermophilus_JW_NM_WN_LF'
        folder = re.compile('[^\w\d]+').sub('_', self.name).strip('_')
        rel_path = os.path.join(self.subdir, folder)
        os.makedirs(self.gtdb.get_full_path_to(rel_path))
        self.add_param('dir_path', rel_path)

    def _create_seq_file(self, record, fmt):
        "Saves sequence in the org dir in the format specificed as 'fmt'."
        fmt2dir = {'fasta-2line': 'seq_fna', 'genbank': 'seq_gbk'}

        # Create dir like 'orgs/Delftia_acidovorans_SPH_1/seq_fna' if needed
        org_dir = self.gtdb.get_full_path_to(self.prm['dir_path'])
        full_path = os.path.join(org_dir, fmt2dir[fmt])
        os.makedirs(full_path, exist_ok=True)

        file_path = os.path.join(full_path, record.id)
        SeqIO.write(record, file_path, fmt)

    def get_full_path_to(self, *args):
        "Returns a full path for a path relative to ORG dir."
        return self.gtdb.get_full_path_to(self.prm['dir_path'], *args)

    def make_all_params(self):
        "Generates/updates the majority of org params."
        # Get all org gbk files with full paths
        all_gbk = os.listdir(self.get_full_path_to('seq_gbk'))
        all_gbk = [self.get_full_path_to('seq_gbk', fn) for fn in all_gbk]

        self.set_param('num_seqs', num=len(all_gbk))

        record = SeqIO.read(all_gbk[0], "genbank")
        self._make_param_taxonomy(record)

        #print(record.dbxrefs)
        #dbx_dict = db_xref_list_to_dict(record.dbxrefs)
        #for (name,value) in dbx_dict.items():
        #    org.add_param(name, value)

        self._make_param_short_name()

    def _make_param_taxonomy(self, record):
        "Creates the 'taxonomy' param from the given SeqRecord."
        self.delete_param('taxonomy')
        num = 0
        for value in record.annotations.get('taxonomy', []):
            self.add_param('taxonomy', value, num)
            num += 1

    def _make_param_short_name(self):
        "e.g. 'Mycobacterium tuberculosis H37Rv'  => 'M. tuberculosis'"
        self.delete_param('short_name')
        short_name_re = re.compile(r'^(\w)\w*\s+(\w+).*$')
        if short_name_re.match(self.name):
            short_name = short_name_re.sub(r'\1. \2', self.name)
            self.add_param('short_name', short_name)
        else:
            logging.warning("Can't generate short name from '%s'" % self.name)


@receiver(signals.pre_delete, sender=Org)
def on_org_delete(sender, instance, using, **kwargs):
    """Make sure to remove org dir if the org is deleted. This is done using
    Django signals: https://stackoverflow.com/a/12678428/310453
    """
    shutil.rmtree(instance.get_full_path_to())
    print('Deleting org dir: ' + instance.get_full_path_to())

def _get_species_genus_phylum_kingdom(record):
    taxa_l = record.annotations.get('taxonomy')
    if taxa_l is None or len(taxa_l) < 3:
        return None, None, None, None

    kingdom = taxa_l[0]
    if kingdom not in ['Bacteria', 'Archaea', 'Eukaryota']:
        logging.warning("Unknown kingdom '%s'" % kingdom)
        return None, None, None, None

    species = record.annotations['organism']
    if taxa_l[-1] in species:
        genus = taxa_l[-1]
    elif taxa_l[-2] in species:
        genus = taxa_l[-2]
    else:
        logging.warning("Can't determine genus for species name '%s' and taxonomy %s" %
                        (species, taxa_l))
        genus = None

    if genus is not None and ' ' in genus:
        # Genus must be a single word!
        logging.warning("Wrong genus '%s' for species name '%s'" % (genus, species))
        genus = None

    phylum = taxa_l[1]
    return species, genus, phylum, kingdom


class OrgParam(AbstractParam):
    parent = models.ForeignKey(Org, on_delete=models.CASCADE)

    class Meta:
        db_table = 'org_params'

