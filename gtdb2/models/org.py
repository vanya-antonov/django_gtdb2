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

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    prm_info = dict(list(AbstractUnit.prm_info.items()) + list({
        'num_seqs': {'value_attr': 'num', 'type_fun': int},
        'taxonomy': {'value_attr': 'value', 'is_list': True, 'sort_attr': 'num'},
        'transl_table': {'value_attr': 'value', 'type_fun': int, 'is_list': True,
                         'sort_attr': 'num', 'reverse': True},
    }.items()))

    dir_path = property(
        fget=lambda self: self.gtdb.get_full_path_to(self.prm['dir_path']),
        doc="Returns a full path to the org dir.")

    def get_all_seq_ids(self, seq_dir='seq_gbk', fullpath=False):
        """Returns a list of all seq ids that are the names of the files in
        the seq_gbk/seq_fna folder.

        Arguments:
         - seq_dir - name of the folder with sequence files
         - fullpath - if True appends full paths to the sequence files.
        """
        seq_dir = os.path.join(self.dir_path, seq_dir)
        all_ids = os.listdir(seq_dir)
        if fullpath:
            all_ids = [os.path.join(seq_dir, fn) for fn in all_ids]
        return all_ids

    def make_all_params(self):
        "Generates/updates the majority of org params."
        self._make_param_short_name()

        # Get all org gbk files with full paths
        all_gbk = self.get_all_seq_ids(seq_dir='seq_gbk', fullpath=True)
        self.set_param('num_seqs', num=len(all_gbk))

        record = SeqIO.read(all_gbk[0], "genbank")
        self._make_param_taxonomy(record)

        self._make_param_xref(record)

    def update_seqs_with_gbk(self, gbk_fn):
        """Compares current list of org seqs with the seqs from the gbk file.
        Removes seqs that are not present in gbk. Adds seqs that are present
        in gbk, but not in org_dir and replaces seqs with records from the
        file if they have newer versions. Returns a tuple with info about
        the update statistics.
        """
        # gbk file may be very large, so get a list of IDs first
        gbk_seq_ids = [r.id for r in SeqIO.parse(gbk_fn, "genbank")]
        org_seq_ids = self.get_all_seq_ids()
        if set(org_seq_ids) == set(gbk_seq_ids):
            return {"n_new": 0, "n_updated": 0, "n_deleted": 0}
        raise NotImplementedError("Some org seqs should be updated!")

    @classmethod
    def get_or_create_from_gbk(cls, user, gbk_fn):
        "Returns existing or a newly created Org object."
        record = next(SeqIO.parse(gbk_fn, "genbank"))
        org = cls.get_by_SeqRecord(record)
        if org is None:
            org = cls.create_from_gbk(user, gbk_fn)
        return org

    @classmethod
    def get_by_SeqRecord(cls, record):
        "Returns existing Org object or None."
        org_name = record.annotations['organism']
        org = cls.objects.filter(name=org_name).first()

        # TODO: validate org and record xrefs!
        return org

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
                raise ValueError("Sequence '%s' does NOT belong to org '%s'" %
                                 (record.id, org.name))
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
            raise ValueError("Can't determine species, genus, phylum or "
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

    def _make_param_xref(self, record):
        "Creates the 'xref' params from the given SeqRecord."
        all_xrefs = record.dbxrefs
        all_xrefs += _get_source_feature_xrefs(record)
        for xref in all_xrefs:
            self.add_gbk_xref_param(xref)

    def _make_param_taxonomy(self, record):
        "Creates the 'taxonomy' params from the given SeqRecord."
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


class OrgParam(AbstractParam):
    parent = models.ForeignKey(Org, on_delete=models.CASCADE,
                               related_name='param_set')

    class Meta:
        db_table = 'org_params'


@receiver(signals.pre_delete, sender=Org)
def on_org_delete(sender, instance, using, **kwargs):
    """Make sure to remove org dir if the org is deleted. This is done using
    Django signals: https://stackoverflow.com/a/12678428/310453
    """
    shutil.rmtree(instance.dir_path)

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

def _get_source_feature_xrefs(record):
    "Feature with type='source' has an important xref = 'taxon:1210089'."
    for f in record.features:
        if f.type == 'source':
            return f.qualifiers['db_xref']
    return []

