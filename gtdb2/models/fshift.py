# Copyright 2018 by Ivan Antonov. All rights reserved.

from pprint import pprint

from django.db import models

from gtdb2.models.abstract import AbstractUnit, AbstractParam
from gtdb2.models.seq import Seq


class FShift(AbstractUnit):
    seq = models.ForeignKey(Seq, on_delete=models.CASCADE)
    coord = models.IntegerField()
    type = models.IntegerField()
    origin = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()

    class Meta:
        db_table = 'fshifts'

    @classmethod
    def create_from_gtdb1_fs(cls, user, org, gtdb1_fs):
        if gtdb1_fs.job.species not in org.name:
            raise ValueError("Provided gtdb1_fs and org do not match!")

#        pprint(vars(gtdb1_fs))
#        pprint(vars(gtdb1_fs.job))
#        pprint(vars(org))

'''
def create_new_in_db_from_GT_FS(gtdb, fs_id, source='genetack'):
    fs, seq = _get_info_about_GT_FS(gtdb, fs_id)
    
    # If up_len_nt is not divisible by 3 => the FS was predicted in a middle of codon.
    # Move the fs-coord upstream to avoid stop codon in cases like 'aaa_tgA_CCC'.
    adjust_len = fs['up_len_nt'] % 3
    if fs['strand'] == 1:
        fs['fs_coord'] -= adjust_len
    else:
        fs['fs_coord'] += adjust_len
    
    fsgene = FSGene.create_new_in_db(
        gtdb, seq.id, fs['user_id'], fs['start'], fs['end'], fs['strand'], fs['fs_coord'], fs['type'],
        source=source, cof_id=fs['cof_id'], c_date=fs['c_date'], db_id=fs_id)
    fsgene.make_all_params()
    
def _get_info_about_GT_FS(gtdb, fs_id):
    res = gtdb.exec_sql_ar("""
        select f.fs_id, j.user_id, j.c_date, s.id AS seq_id, s.ext_id AS seq_ext_id,
               f.fs_coord, f.type, f.init_gene_seq,
               IF(f.strand = "+", 1, -1) AS strand,
               (select cof_id from cof_gtfs cg where cg.fs_id = f.fs_id limit 1) AS cof_id
        from gt_fs f, jobs j, seqs s
        where j.job_id=f.job_id
        and s.name = j.name
        and f.fs_id = %s
    """, fs_id
    )
    if len(res) == 0:
        logging.warning("Can't get info about FS_ID = '%s'" % fs_id)
        return None
    elif not res[0]['seq_ext_id'].endswith('.1'):
        logging.warning("Sequence '%s' has wrong version!" % res[0]['seq_ext_id'])
        return None
    else:
        fs = res[0]
        fs['type'] = int(fs['type'])
    
    if fs['strand'] == -1:
        fs['fs_coord'] -= 1   # switch to zero-based coordinate system?
    
    up_match = re.compile('^[a-z]+').search(fs['init_gene_seq'])
    down_match = re.compile('[A-Z]+$').search(fs['init_gene_seq'])
    if up_match is None or down_match is None:
        logging.warning("Wrong fsgene sequence '%s'" % fs['init_gene_seq'])
        return None
    fs['up_len_nt'] = up_match.end() - up_match.start()
    fs['down_len_nt'] = down_match.end() - down_match.start()
    
    if fs['strand'] == 1:
        fs['start'] = fs['fs_coord'] - fs['up_len_nt']
        fs['end'] = fs['fs_coord'] + fs['down_len_nt']
    else:
        fs['start'] = fs['fs_coord'] - fs['down_len_nt']
        fs['end'] = fs['fs_coord'] + fs['up_len_nt']
    
    seq = SeqWithSeq(gtdb, fs['seq_id'])
    fsgene_seq = seq[fs['start']:fs['end']]
    if fs['strand'] == -1:
        fsgene_seq = fsgene_seq.reverse_complement()
    
    if fsgene_seq.upper() != fs['init_gene_seq'].upper():
        logging.warning("GT_FS and fsgene sequence do not match:\n%s\n%s" % (fs['init_gene_seq'], fsgene_seq))
        return None
    
    return fs, seq
'''


class FShiftParam(AbstractParam):
    parent = models.ForeignKey(FShift, on_delete=models.CASCADE)

    class Meta:
        db_table = 'fshift_params'

