# Copyright 2018 by Ivan Antonov. All rights reserved.

import io
import subprocess

from Bio.Blast import NCBIXML

from .config import MAX_NUM_THREADS, DEFAULT_EVALUE


def run_blastp_vs_file(prot_seq, file_path, evalue=DEFAULT_EVALUE):
    # Use the fixed -dbsize to make results (evalues) reproducible
    cmd_str = ' '.join([
        'blastp  -outfmt 5  -dbsize 1  -max_target_seqs 999999',
        '-evalue', str(evalue),
        '-subject', file_path])
    return _run_blast_cmd_str(prot_seq, cmd_str)

def run_tblastn(prot_seq, blastdb_path, gcode,
                evalue=DEFAULT_EVALUE,
                num_threads=MAX_NUM_THREADS):
    # Use the fixed -dbsize to make results (evalues) reproducible
    cmd_str = ' '.join([
        'tblastn  -outfmt 5  -dbsize 1  -max_target_seqs 999999',
        '-evalue', str(evalue),
        '-num_threads', str(num_threads),
        '-db_gencode', str(gcode),
        '-db', blastdb_path])

    hits_dict = _run_blast_cmd_str(prot_seq, cmd_str)

    # Overwrite the .strand attribute becuase it is None anyway
    for all_hsps in hits_dict.values():
        for hsp in all_hsps:
            hsp.strand = 1 if hsp.frame[1] >= 0 else -1

    return hits_dict

def _run_blast_cmd_str(query_seq, cmd_str):
    """Returns a dict where the keys are the hit sequence IDs and
    the values are the lists of Bio.Blast.Record.HSP objects
    (https://github.com/biopython/biopython/blob/master/Bio/Blast/Record.py).
    The coordinate system of the HSP objects is converted to 0-based.
    """
    fasta_txt = '>query\n' + query_seq + '\n'
    proc = subprocess.run(
        cmd_str, input=fasta_txt, stdout=subprocess.PIPE,
        shell=True, universal_newlines=True)

    f = io.StringIO(proc.stdout)
    all_res = list(NCBIXML.parse(f))
    f.close()
    if(len(all_res) != 1):
        raise ValueError('BLAST XML must contain results for a single query!')

    all_hits = {}
    for alignment in all_res[0].alignments:
        for hsp in alignment.hsps:
            # Hsp objects use 1-based coordinates - let's fix it
            hsp.query_start -= 1
            hsp.sbjct_start -= 1

            all_hits.setdefault(alignment.hit_def, []).append(hsp)

    return all_hits

