# Copyright 2018 by Ivan Antonov. All rights reserved.

import re

from django.db import models

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.seq import Seq


class Fshift(AbstractUnit):
    seq = models.ForeignKey(Seq, on_delete=models.CASCADE)
    type = models.CharField(max_length=255)
    coord = models.IntegerField()
    len = models.IntegerField()
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()

    class Meta:
        db_table = 'fshifts'

    @classmethod
    def create_from_gtdb1_fs(cls, user, seq, gtdb1_fs):
        "Creates Fshift based on the info from GTDB1 fs."
        fs = _validate_gtdb1_fs(gtdb1_fs, seq)

        fshift = Fshift(
            user=user, seq=seq, c_date=fs.job.c_date, type='genetack',
            name=_get_fshift_name(seq, fs.fs_coord, int(fs.type)),
            descr=fs.descr, coord=fs.fs_coord, len=int(fs.type),
            start=fs.start, end=fs.end, strand=fs.strand)
        fshift.save()

        # add gtdb1 ID to xrefs
        fshift.add_xref_param('gtdb1', fs.fs_id)

        return fshift


class FshiftParam(AbstractParam):
    parent = models.ForeignKey(Fshift, on_delete=models.CASCADE)

    class Meta:
        db_table = 'fshift_params'


def _get_fshift_name(seq, fs_coord, fs_len):
    "Generates name like: 'NC_010002.1:1922457:-1'."
    if fs_len == 0:
        # to avoid 'NC_005085.1:1684916:+0'
        return '%s:%d:%d' % (seq.id, fs_coord, fs_len)
    else:
        return '%s:%d:%+d' % (seq.id, fs_coord, fs_len)

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

