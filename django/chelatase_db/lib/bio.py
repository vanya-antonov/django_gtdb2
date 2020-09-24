# Copyright 2018 by Ivan Antonov. All rights reserved.

import io
import logging
import re
import subprocess

from Bio.Blast import NCBIXML
from Bio.Data import CodonTable

from .config import MAX_NUM_THREADS, DEFAULT_EVALUE


def run_blastp_seq_vs_file(query_seq, target_fn, evalue=DEFAULT_EVALUE):
    """Arguments:
     - query_seq - (str) query protein sequence
     - target_fn - (str) full path to a file with other protein seqs
    """
    # Use the fixed -dbsize to make results (evalues) reproducible
    cmd_str = ' '.join([
        'blastp  -outfmt 5  -dbsize 1  -max_target_seqs 999999',
        '-evalue', str(evalue),
        '-subject', target_fn])
    xml_str = _run_blast_cmd_str_with_query_seq(cmd_str, query_seq)
    return _parse_blast_xml_str(xml_str)

def run_tblastn_seq_vs_db(query_seq, target_db, gcode,
                evalue=DEFAULT_EVALUE):
    """Arguments:
     - query_seq - (str) query protein sequence
     - target_db - (str) full path to the blast db
     - gcode - (int) target db genetic code
    """
    # Use the fixed -dbsize to make results (evalues) reproducible
    cmd_str = ' '.join([
        'tblastn  -outfmt 5  -dbsize 1  -max_target_seqs 999999',
        '-evalue', str(evalue),
        '-num_threads', str(MAX_NUM_THREADS),
        '-db_gencode', str(gcode),
        '-db', target_db])
    xml_str = _run_blast_cmd_str_with_query_seq(cmd_str, query_seq)
    return _parse_blast_xml_str(xml_str)

def run_tblastn_file_vs_db(query_fn, target_db, gcode,
                           evalue=DEFAULT_EVALUE):
    """Arguments:
     - query_fn - (str) full path to the fn with query sequence(s)
     - target_db - (str) full path to the blast db
     - gcode - (int) target db genetic code
    """
    # Use the fixed -dbsize to make results (evalues) reproducible
    cmd_str = ' '.join([
        'tblastn  -outfmt 5  -dbsize 1  -max_target_seqs 999999',
        '-evalue', str(evalue),
        '-num_threads', str(MAX_NUM_THREADS),
        '-db_gencode', str(gcode),
        '-query', query_fn,
        '-db', target_db])
    xml_str = subprocess.getoutput(cmd_str)
    return _parse_blast_xml_str(xml_str)

def _run_blast_cmd_str_with_query_seq(cmd_str, query_seq):
    """Returns XML string produced by BLAST run.
    """
    fasta_txt = '>query\n' + query_seq + '\n'
    return subprocess.run(
        cmd_str, input=fasta_txt, stdout=subprocess.PIPE,
        shell=True, universal_newlines=True
    ).stdout

def _parse_blast_xml_str(xml_str):
    """Returns a sorted list of Bio.Blast.Record.HSP objects with some extra
    attributes -- query_id, sbjct_id, sbjct_strand. Also the coordinate system
    of all the HSP objects is converted to 0-based.
    https://github.com/biopython/biopython/blob/master/Bio/Blast/Record.py
    """
    if xml_str == '':
        return []

    f = io.StringIO(xml_str)

    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc102
    all_hsps = []
    for blast_record in NCBIXML.parse(f):
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                # Hsp objects use 1-based coordinates - let's fix it
                hsp.query_start -= 1
                hsp.sbjct_start -= 1

                # Add the query and the hit IDs
                hsp.query_id = blast_record.query
                hsp.sbjct_id = alignment.hit_def

                # Add a more convenient attribute for the hit strand
                hsp.sbjct_strand = 1 if hsp.frame[1] >= 0 else -1

                all_hsps.append(hsp)

    f.close()

    return sorted(all_hsps, key = lambda hsp: hsp.expect)

def get_stop_stop_seq(left, right, strand, chr_seq, gencode):
    """Expand the given CDS part in both directions until stop codons
    on both sides (the stop codons are included as well).
    """
    # Validate the provided coordinates
    if (right - left) % 3 != 0:
        raise ValueError("The difference between the left and right "
                         "coordinates is not divisible by 3!!")

    chr_region = chr_seq[left:right]
    if not re.compile('^[ACGT]+$', re.IGNORECASE).match(str(chr_region)):
        logging.warning("Sequence %i-%i contains non-ACGT chars" %
                        (left, right))
        return None

    if strand == -1:
        chr_region = chr_region.reverse_complement()
    prot_seq = chr_region.translate(table=gencode)
    if '*' in prot_seq:
        logging.warning('Given region (%i,%i) already contains a stop codon' %
                        (left, right))
        return None

    ss_left = left
    while ss_left >= 3:
        cur_codon = chr_seq[ (ss_left-3) : ss_left ]
        codon_type = get_codon_type(cur_codon, gencode, strand)
        if codon_type == 'coding':
            ss_left -= 3  # expand the stop-stop seq
        elif codon_type == 'stop':
            ss_left -= 3  # stop-stop seq includes stops as well
            break
        else:
            break

    ss_right = right
    while ss_right < len(chr_seq)-3:
        cur_codon = chr_seq[ ss_right : (ss_right+3) ]
        codon_type = get_codon_type(cur_codon, gencode, strand)
        if codon_type == 'coding':
            ss_right += 3  # expand the stop-stop seq
        elif codon_type == 'stop':
            ss_right += 3  # stop-stop seq includes stops as well
            break
        else:
            break

    return {'left' : ss_left, 'right' : ss_right}

def get_codon_type(codon, gencode, strand=1):
    acgt_only_re = re.compile('^[ACGT]+$', re.IGNORECASE)
    if not acgt_only_re.match(str(codon)):
        return 'non_ACGT'

    if strand == -1:
        codon = codon.reverse_complement()

    if codon.upper() in CodonTable.unambiguous_dna_by_id[gencode].stop_codons:
        return 'stop'
    else:
        return 'coding'
