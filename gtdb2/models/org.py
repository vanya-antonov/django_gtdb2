# Copyright 2018 by Ivan Antonov. All rights reserved.

import logging
import os
import shutil
import subprocess
import sys
import re

from Bio import SeqIO

from django.apps import apps
from django.db import models
from django.db.models import signals
from django.dispatch import receiver

from gtdb2.lib.db import GeneTackDB
from gtdb2.models.abstract import AbstractUnit, AbstractParam


class Org(AbstractUnit):
    # Overwrite AbstractUnit attribute - name must be unique and NOT NULL
    name = models.CharField(max_length=255, unique=True)

    genus = models.CharField(max_length=255)
    phylum = models.CharField(max_length=255)
    kingdom = models.CharField(max_length=255)

    class Meta:
        db_table = 'orgs'

    # Define these attributes so they can be overwritten by the derived classes
    FEAT_CLS_NAME = 'Feat'
    FSHIFT_CLS_NAME = 'Fshift'

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(AbstractUnit.PRM_INFO.items()) + list({
        'blastdb_nucl_all': {},
        'dir_name': {},
        'num_seqs': {'value_attr': 'num', 'type_fun': int},
        'short_name': {},
        'source_fn': {},
        'taxonomy': {'is_list': True, 'sort_attr': 'num'},
        'transl_table': {'type_fun': int, 'is_list': True,
                         'sort_attr': 'num', 'reverse': True},
    }.items()))

    # The '%'s will be substituted with the value of 'dir_name' prm
    SUBDIR_INFO = {
        'main': 'orgs/%s/',
        'seq_fna': 'orgs/%s/seq_fna/',
        'seq_gbk': 'orgs/%s/seq_gbk/',
        'blastdb': 'orgs/%s/blastdb/',}

    @property
    def feat_set(self):
        """Retruns a QuerySet of all the Feat objects (or its subclasses)
        that belong to the current org.
        """
        return self._get_grandchildren_set(self.FEAT_CLS_NAME)

    @property
    def fshift_set(self):
        """Retruns a QuerySet of all Fshift objects (or its subclasses)
        that belong to the current org.
        """
        return self._get_grandchildren_set(self.FSHIFT_CLS_NAME)

    @property
    def transl_table(self):
        """Returns a translation table most typical to the org. For example,
        for human it should return the standard genetic code and not the
        mitochondrial genetic code.
        """
        return self.prm['transl_table'][0]

    @classmethod
    def get_or_create_from_gbk(cls, user, gbk_fn):
        "Returns existing or a newly created Org object."
        org = cls.get_by_gbk(gbk_fn)
        if org is None:
            org = cls.create_from_gbk(user, gbk_fn)
        return org

    @classmethod
    def get_by_gbk(cls, gbk_fn):
        "Returns existing Org object or None."
        record = next(SeqIO.parse(gbk_fn, "genbank"))
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
            elif record.annotations['organism'] != org.name:
                # Make sure that all seqs correspond to the same org
                raise ValueError("Sequence '%s' does NOT belong to org '%s'" %
                                 (record.id, org.name))
            org._create_seq_file(record, 'fasta-2line', 'seq_fna')
            org._create_seq_file(record, 'genbank', 'seq_gbk')

        org.set_param('source_fn', fn)

        org.create_all_params()
        org.create_annotation()

        return org

    @classmethod
    def _create_from_SeqRecord(cls, user, record):
        "Creates a new org in db, creates org dir and adds some basic params."
        sp, ge, ph, ki = _get_species_genus_phylum_kingdom(record)
        if None in (sp, ge, ph, ki):
            raise ValueError("Can't determine species, genus, phylum or "
                             "kingdom: %s, %s, %s, %s" % (sp, ge, ph, ki))
        org = cls(user=user, name=sp, genus=ge, phylum=ph, kingdom=ki)
        org.save()

        org._create_org_dir()

        return org

    def get_full_path_to_subdir(self, name='main'):
        "Returns a full path to subdir by its name or alias."
        subdir = self.SUBDIR_INFO[name] % self.prm['dir_name']
        return self.gtdb.get_full_path_to(subdir)

    def get_all_seq_ids(self, seq_dir='seq_gbk', fullpath=False):
        """Returns a list of all seq ids that are the names of the files in
        the seq_gbk/seq_fna folder.

        Arguments:
         - seq_dir - name of the folder with sequence files
         - fullpath - if True appends full paths to the sequence files.
        """
        seq_dir = self.get_full_path_to_subdir(seq_dir)
        all_ids = os.listdir(seq_dir)
        if fullpath:
            all_ids = [os.path.join(seq_dir, fn) for fn in all_ids]
        return all_ids

    def read_seq_file(self, ext_id, seq_dir='seq_gbk'):
        "Returns a SeqRecord object."
        dir2fmt = {'seq_gbk': 'genbank', 'seq_fna': 'fasta'}
        dir_path = self.get_full_path_to_subdir(seq_dir)
        seq_path = os.path.join(dir_path, ext_id)
        return SeqIO.read(seq_path, dir2fmt[seq_dir])

    def create_all_params(self):
        "Generates/updates the majority of params."
        self._make_param_short_name()

        # Get all org gbk files with full paths
        all_gbk = self.get_all_seq_ids(seq_dir='seq_gbk', fullpath=True)
        self.set_param('num_seqs', num=len(all_gbk))

        record = SeqIO.read(all_gbk[0], "genbank")
        self._make_param_taxonomy(record)
        self._make_param_xref(record)

        self._make_param_transl_table()

        self._make_param_blastdb()
        #org.make_genetack_model()

    def create_annotation(self):
        # create annotated fshifts
        # create genetack fshifts
        # create 16S rRNA feats
        pass

    def update_seqs_with_gbk(self, gbk_fn):
        """Compares current list of org seqs with the seqs from the gbk file.
        Removes seqs that are not present in gbk. Adds seqs that are present
        in gbk, but not in org_dir and replaces seqs with records from the
        file if they have newer versions. Returns a list of new seqs (i.e.
        the ones that were updated or just created).
        """
        # gbk file may be very large, so get a list of IDs first
        gbk_seq_ids = [r.id for r in SeqIO.parse(gbk_fn, "genbank")]
        org_seq_ids = self.get_all_seq_ids()
        if set(org_seq_ids) == set(gbk_seq_ids):
            return []
            #return {"n_new": 0, "n_updated": 0, "n_deleted": 0}
        raise NotImplementedError("Some org seqs should be updated!")

    def _get_grandchildren_set(self, grandchildren_cls_name):
        """Retruns a QuerySet of all grandchildren of particular type
        (assuming that seqs are children of the org).
        For example, if grandchildren_cls_name = 'Feat' it will
        return a QuerySet of Feat objects, i.e.
        ORG (self)  ==>  ALL_SEQS  ==>  ALL_FEATS.

        Arguments:
         - grandchildren_cls_name - a string corresponding to grandchildren
         class name (e.g. 'Feat' or 'Fshift')
        """
        # Get the app name, i.e. 'gtdb2':
        # https://stackoverflow.com/a/2742722/310453
        app_label = self._meta.app_label

        # Get the class object for the db model:
        # https://stackoverflow.com/a/36234846/310453
        grandchildren_cls = apps.get_model(app_label, grandchildren_cls_name)

        # Create the QuerySet
        return grandchildren_cls.objects.filter(seq__org=self)

    def _create_org_dir(self):
        "Creates dir for the new org and saves its name in the params table."
        # 'Natranaerobius thermophilus JW/NM-WN-LF'  =>
        # 'Natranaerobius_thermophilus_JW_NM_WN_LF'
        dir_name = re.compile('[^\w\d]+').sub('_', self.name).strip('_')
        self.set_param('dir_name', dir_name)

        dir_path = self.get_full_path_to_subdir()
        os.makedirs(dir_path)

    def _create_seq_file(self, record, fmt, subdir):
        "Saves sequence in the org dir in the format specificed as 'fmt'."
        # Create dir like 'orgs/Delftia_acidovorans_SPH_1/seq_fna' if needed
        full_path = self.get_full_path_to_subdir(subdir)
        os.makedirs(full_path, exist_ok=True)

        file_path = os.path.join(full_path, record.id)
        SeqIO.write(record, file_path, fmt)

    def _make_param_xref(self, record):
        "Creates the 'xref' params from the given SeqRecord."
        all_xrefs = record.dbxrefs
        all_xrefs += _get_source_feature_xrefs(record)
        for xref in all_xrefs:
            # xref is a string like 'Assembly:GCF_001613165.1'
            parts = xref.split(':')
            if len(parts) != 2:
                raise ValueError("Wrong gbk_xref='%s'" % xref)
            self.set_param_xref(parts[0], parts[1])

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
            self.set_param('short_name', short_name)
        else:
            logging.warning("Can't generate short name from '%s'" % self.name)

    def _make_param_transl_table(self):
        """Analyzes all features from all org seqs and creates transl_table
        param with statistics about all genetic codes of the org.
        """
        all_gcodes = {}
        for gbk_fn in self.get_all_seq_ids(seq_dir='seq_gbk', fullpath=True):
            record = SeqIO.read(gbk_fn, "genbank")
            for f in record.features:
                gcode = f.qualifiers.get('transl_table', [None])[0]
                if gcode is not None:
                    all_gcodes.setdefault(gcode, 0)
                    all_gcodes[gcode]+= 1

        self.delete_param('transl_table')
        for gcode, num_f in all_gcodes.items():
            self.add_param('transl_table', value=gcode, num=num_f)

    def _make_param_blastdb(self):
        "Creates all blastdbs and saves them as org params."
        # Remove any existing DBs
        blastdb_dir = self.get_full_path_to_subdir('blastdb')
        if os.path.exists(blastdb_dir):
            # About ignore_errors: https://stackoverflow.com/a/303225/310453
            shutil.rmtree(blastdb_dir, ignore_errors=True)
        os.makedirs(blastdb_dir)   # Create empty dir

        # Create blastdbs inside it
        self._make_param_blastdb_nucl_all(blastdb_dir)
        #self._make_param_blastdb_nt(blastdb_dir)

    def _make_param_blastdb_nucl_all(self, blastdb_dir):
        """Creates a blastdb from all org seqs without considering their
        genetic codes. It is equivalent to genome for prokaryotes and to
        transcriptome for eukaryotes.
        """
        db_path = os.path.join(blastdb_dir, 'nucl_all')
        cmd_str = 'makeblastdb -dbtype nucl -title "%s" -out %s' % (
            self.name, db_path)
        p = subprocess.Popen(
            cmd_str, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, shell=True, universal_newlines=True)

        all_ext_ids = self.get_all_seq_ids()
        for ext_id in all_ext_ids:
            record = self.read_seq_file(ext_id)
            p.stdin.write('>' + record.id + '\n' +
                          str(record.seq) + '\n')
        p.stdin.close()

        # Save relative path to db as param
        rel_path = os.path.relpath(db_path, self.gtdb.root_dir)
        return self.set_param('blastdb_nucl_all', rel_path,
                              num=len(all_ext_ids))


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
    dir_path = instance.get_full_path_to_subdir()
    if os.path.exists(dir_path):
        shutil.rmtree(instance.get_full_path_to_subdir())
    else:
        logging.error("Directory doesn't exist: %s" % dir_path)

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

