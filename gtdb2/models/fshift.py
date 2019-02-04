# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint
import re

import Bio.Seq
from Bio.Alphabet import generic_dna, generic_protein
from django.db import models

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.seq import Seq


class Fshift(AbstractUnit):
    ORIGIN_CHOICES = (
        ('gbk_annotation', 'GenBank annotation'),
        ('genetack', 'GeneTack prediction'),
        ('tblastn', 'tBLASTn prediction'),)
    seq = models.ForeignKey(Seq, on_delete=models.CASCADE)
    origin = models.CharField(max_length=255, choices=ORIGIN_CHOICES)
    coord = models.IntegerField()
    len = models.IntegerField()
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()
    seed = models.BooleanField(default=False)

    class Meta:
        db_table = 'fshifts'

    # Merge with parent prm_info: https://stackoverflow.com/a/38990/310453
    PRM_INFO = dict(list(AbstractUnit.PRM_INFO.items()) + list({
        'seq_nt': {'value_attr': 'data'},
        'seq_nt_corr': {'value_attr': 'data'},
        'seq_nt_n': {'value_attr': 'data'},
        'seq_nt_c': {'value_attr': 'data'},
        'seq_prot': {'value_attr': 'data'},
        'seq_prot_n': {'value_attr': 'data'},
        'seq_prot_c': {'value_attr': 'data'},
    }.items()))

    @property
    def org(self):
        return self.seq.org

    @classmethod
    def get_or_create(cls, *args, **kwargs):
        """Tries to find similar fshift in the DB and, if fails, creates a
        new one.
        """
        all_fshifts = cls.objects.filter(
            seq=kwargs['seq'], strand=kwargs['strand'], len=kwargs['len'],
            coord__gt=kwargs['start'], coord__lt=kwargs['end'],
        ).all()
        if len(all_fshifts) == 0:
            # Create a new fshift
            self = cls(*args, **kwargs)
            self.create_all_params()
            return self
        elif len(all_fshifts) == 1:
            # Return existing fshift
            return all_fshifts[0]
        else:
            raise NotImplemented("Fshift with the closest coord should "
                                 "be returned!")

    @classmethod
    def create_from_gtdb1_fs(cls, user, seq, gtdb1_fs):
        "Creates Fshift based on the info from GTDB1 fs."
        fs = _validate_gtdb1_fs(gtdb1_fs, seq)

        self = cls(
            user=user, seq=seq, c_date=fs.job.c_date, origin='genetack',
            descr=fs.descr, coord=fs.fs_coord, len=int(fs.type),
            start=fs.start, end=fs.end, strand=fs.strand)
        self.save()

        # add gtdb1 ID to xrefs
        self.set_param_xref('gtdb1', fs.fs_id)

        self.create_all_params()

        return self

    def get_prm_bio_seq(self, name):
        "Returns Bio.Seq.Seq object for a given seq param."
        if name.startswith('seq_nt'):
            alphabet = generic_dna
        elif name.startswith('seq_prot'):
            alphabet = generic_protein
        else:
            raise ValueError("Can't determine alphabet for prm '%s'" % name)
        return Bio.Seq.Seq(self.prm[name], alphabet)

    def set_seq_param(self, name, seq):
        "Saves seq as param-data and len(seq) as param-num."
        self.set_param(name, data=seq, num=len(seq))

    def create_all_params(self):
        "Generates/updates name attribute and the majority of params."
        self._make_fshift_name()
        self._make_param_seq_nt_and_prot()

    def _make_fshift_name(self):
        "Generates name like: 'NC_010002.1:1922457:-1'."
        self.name = '%s:%d:%+d' % (self.seq.id, self.coord, self.len)
        self.save()

    def _make_param_seq_nt_and_prot(self):
        "Creates seq_nt* and seq_prot* params."
        start, end, strand = self.start, self.end, self.strand
        fs_coord, fs_type = self.coord, self.len
        seq = self.seq.seq

        up_len_nt = fs_coord - start
        down_len_nt = end - fs_coord
        if strand == -1:
            up_len_nt, down_len_nt = down_len_nt, up_len_nt

        # Backward frameshiting increases the C-terminal length of fs-protein
        down_len_nt -= fs_type

        if up_len_nt % 3 != 0 or down_len_nt % 3 != 0:
            raise ValueError(
                "The length of the upstream / downstream fsgene part "
                "(%d / %d) is not divisible by 3: start=%d, end=%d, "
                "fs_coord=%d, fs_type=%+d, strand=%d" %
                (up_len_nt, down_len_nt, start, end, fs_coord, fs_type, strand))

        chunk_nt = seq[start:end]
        if strand == -1:
            chunk_nt = chunk_nt.reverse_complement()

        up_chunk_nt = chunk_nt[:up_len_nt]
        down_chunk_nt = chunk_nt[-down_len_nt:]
        prot_seq_n = up_chunk_nt.translate(table = self.seq.transl_table)
        prot_seq_c = down_chunk_nt.translate(table = self.seq.transl_table)
        prot_seq_c = prot_seq_c.rstrip('*')  # remove possible stop codon at the end

        if '*' in prot_seq_n + prot_seq_c:
            raise ValueError("FS-prot seq contains in-frame stop codon:\n"
                             "%s\n\n%s" % (prot_seq_n, prot_seq_c))

        self.set_seq_param('seq_prot_n', prot_seq_n.upper())
        self.set_seq_param('seq_prot_c', prot_seq_c.upper())
        self.set_seq_param('seq_nt_n', up_chunk_nt.upper())
        self.set_seq_param('seq_nt_c', down_chunk_nt.upper())

        seq_prot = prot_seq_n.lower() + prot_seq_c.upper()
        self.set_seq_param('seq_prot', seq_prot)

        seq_nt = chunk_nt[:up_len_nt].lower() + chunk_nt[up_len_nt:].upper()
        self.set_seq_param('seq_nt', seq_nt)

        seq_nt_corr = up_chunk_nt.lower() + down_chunk_nt.upper()
        self.set_seq_param('seq_nt_corr', seq_nt_corr)


class FshiftParam(AbstractParam):
    parent = models.ForeignKey(Fshift, on_delete=models.CASCADE,
                               related_name='param_set')

    class Meta:
        db_table = 'fshift_params'


def _validate_gtdb1_fs(fs, seq):
    """Makes sure that GTDB1 fs matches GTDB2 seq and modifies/adds some
    attributes to the fs object.
    """
    if fs.job.species not in seq.org.name:
        raise ValueError("Provided gtdb1_fs.species and seq.org do not match:"
                         "\n%s\n%s" %(fs.job.species, seq.org.name))

    fs.strand = -1 if fs.strand == '-' else 1

    if fs.strand == -1:
        fs.fs_coord -= 1   # switch to zero-based coordinate system

    up_match = re.compile('^[a-z]+').search(fs.init_gene_seq)
    down_match = re.compile('[A-Z]+$').search(fs.init_gene_seq)
    if up_match is None or down_match is None:
        raise ValueError("Wrong fsgene sequence '%s'" % fs.init_gene_seq)
    fs.up_len_nt = up_match.end() - up_match.start()
    fs.down_len_nt = down_match.end() - down_match.start()

    if fs.strand == 1:
        fs.start = fs.fs_coord - fs.up_len_nt
        fs.end = fs.fs_coord + fs.down_len_nt
    else:
        fs.start = fs.fs_coord - fs.down_len_nt
        fs.end = fs.fs_coord + fs.up_len_nt

    fsgene_seq = seq.record.seq[fs.start:fs.end]
    if fs.strand == -1:
        fsgene_seq = fsgene_seq.reverse_complement()

    if fsgene_seq.upper() != fs.init_gene_seq.upper():
        raise ValueError("GT_FS and fsgene sequences do not match:\n%s\n%s" %
                         (fs['init_gene_seq'], fsgene_seq))

    # If up_len_nt is not divisible by 3 => the FS was predicted in a
    # middle of a codon => Move the fs-coord upstream to avoid stop
    # codon in cases like 'aaa_tgA_CCC'.
    adjust_len = fs.up_len_nt % 3
    if fs.strand == 1:
        fs.fs_coord -= adjust_len
    else:
        fs.fs_coord += adjust_len

    return fs

