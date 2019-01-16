# Copyright 2018 by Ivan Antonov. All rights reserved.

import re

from Bio.Data.CodonTable import unambiguous_dna_by_id


def is_stop(codon, strand=1, gcode=1):
    """Returns True if the codon is a stop codon and False if it is
    a coding codon.
    """
    if len(codon) != 3:
        raise ValueError("The codon length must be 3!")

    codon = codon.upper()
    acgt_only_re = re.compile('^[ACGT]+$', re.IGNORECASE)
    if not acgt_only_re.match(str(codon)):
        raise ValueError("The codon '%s' contains non-ACGT letters!" % codon)

    if strand == -1:
        codon = codon.reverse_complement()

    return codon in unambiguous_dna_by_id[gcode].stop_codons

