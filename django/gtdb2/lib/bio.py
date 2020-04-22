# Copyright 2018 by Ivan Antonov. All rights reserved.

import re

from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqFeature import FeatureLocation, CompoundLocation


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

def get_overlapping_feats_from_list(all_feats, start, end, strand,
                                    all_types=None, min_overlap=1):
    """Returns a list of objects from f_list that overlap with the specified
    genomic location. The returned objects will be sorted (descending) by the
    value of the '_overlap_len' attribute (int) that will be added to each
    object.

    Arguments:
     - all_feats - a list of feature-objects (e.g. Bio.SeqFeature.SeqFeature or
     gtdb2.models.feat.Feat) that have the .type (string) and the .location
     (Bio.SeqFeature.FeatureLocation) attributes.
    """
    target_loc = FeatureLocation(start, end, strand)
    overlapping_feats = []
    for f in all_feats:
        if all_types is not None and f.type not in all_types:
            continue

        if f.location.end < start or end < f.location.start:
            continue

        f._overlap_len = get_FeatureLocation_overlap_len(f.location, target_loc)
        if f._overlap_len >= min_overlap:
            overlapping_feats.append(f)

    # Sort by overlap_len: https://docs.python.org/3.6/howto/sorting.html
    return sorted(overlapping_feats, reverse=True,
                  key=lambda f: f._overlap_len)

def get_FeatureLocation_overlap_len(f1, f2):
    # https://github.com/biopython/biopython/issues/896
    if not isinstance(f1, (FeatureLocation, CompoundLocation)) or \
       not isinstance(f2, (FeatureLocation, CompoundLocation)):
        raise ValueError("Wrong feature types: %s / %s" %
                         (type(f1), type(f2)))

    if f1.strand is not None and \
       f2.strand is not None and \
       f1.strand != f2.strand:
        return 0

    return len(set(f1).intersection(set(f2)))

