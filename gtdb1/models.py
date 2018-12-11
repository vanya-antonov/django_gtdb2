# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey has `on_delete` set to the desired behavior.
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from django.db import models


class CofGtfs(models.Model):
    cof = models.ForeignKey('Cofs', models.DO_NOTHING)
    fs = models.ForeignKey('GtFs', models.DO_NOTHING)
    is_core = models.IntegerField()

    class Meta:
        managed = False
        db_table = 'cof_gtfs'
        unique_together = (('cof', 'fs'),)

"""
class CofGtfsPairs(models.Model):
    cof = models.ForeignKey('Cofs', models.DO_NOTHING)
    fs_id_1 = models.ForeignKey('GtFs', models.DO_NOTHING, db_column='fs_id_1')
    fs_id_2 = models.ForeignKey('GtFs', models.DO_NOTHING, db_column='fs_id_2')
    type = models.CharField(max_length=128)
    ma = models.IntegerField(db_column='Ma', blank=True, null=True)  # Field name made lowercase.
    ms = models.IntegerField(db_column='Ms', blank=True, null=True)  # Field name made lowercase.
    na = models.IntegerField(db_column='Na', blank=True, null=True)  # Field name made lowercase.
    ns = models.IntegerField(db_column='Ns', blank=True, null=True)  # Field name made lowercase.
    ka = models.FloatField(db_column='Ka', blank=True, null=True)  # Field name made lowercase.
    ks = models.FloatField(db_column='Ks', blank=True, null=True)  # Field name made lowercase.
    kaks = models.FloatField(db_column='KaKs', blank=True, null=True)  # Field name made lowercase.
    mut_cods = models.IntegerField(blank=True, null=True)
    num_codons_wo_gaps = models.IntegerField(blank=True, null=True)
    mut_count = models.IntegerField(blank=True, null=True)
    mut_gap_free_len = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'cof_gtfs_pairs'
        unique_together = (('cof', 'fs_id_1', 'fs_id_2', 'type'),)
"""

class CofHits(models.Model):
    q_fs_id = models.IntegerField()
    q_start = models.IntegerField()
    q_end = models.IntegerField()
    q_fs_ali_coord = models.IntegerField()
    h_fs_id = models.IntegerField()
    h_start = models.IntegerField()
    h_end = models.IntegerField()
    h_fs_ali_coord = models.IntegerField()
    ali_len = models.IntegerField()
    score = models.IntegerField()
    evalue = models.CharField(max_length=255)
    identity = models.FloatField()

    class Meta:
        managed = False
        db_table = 'cof_hits'


class CofLabels(models.Model):
    cof = models.ForeignKey('Cofs', models.DO_NOTHING)
    label_id = models.IntegerField()

    class Meta:
        managed = False
        db_table = 'cof_labels'
        unique_together = (('cof', 'label_id'),)


class CofWords(models.Model):
    cof_id = models.CharField(max_length=255)
    word = models.CharField(max_length=255)
    obs_num = models.FloatField(blank=True, null=True)
    weight = models.FloatField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'cof_words'


class Cofs(models.Model):
    cof_id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=255, blank=True, null=True)
    c_date = models.DateTimeField()
    num_b_fs = models.IntegerField(blank=True, null=True)
    num_p_fs = models.IntegerField(blank=True, null=True)
    num_r_fs = models.IntegerField(blank=True, null=True)
    num_bp_fs = models.IntegerField(blank=True, null=True)
    dir = models.CharField(max_length=255, blank=True, null=True)
    has_paml = models.IntegerField(blank=True, null=True)
    num_fs = models.IntegerField(blank=True, null=True)
    num_genus = models.IntegerField(blank=True, null=True)
    kingdoms = models.CharField(max_length=32, blank=True, null=True)
    fs_dna_num_left_gaps = models.IntegerField(blank=True, null=True)
    fs_dna_num_right_gaps = models.IntegerField(blank=True, null=True)
    fs_dna_num_inner_gaps = models.IntegerField(blank=True, null=True)
    fs_dna_num_identical_cols = models.IntegerField(blank=True, null=True)
    fs_dna_identical_block = models.CharField(max_length=255, blank=True, null=True)
    fs_dna_trimmed_len = models.IntegerField(blank=True, null=True)
    num_species = models.IntegerField(blank=True, null=True)
    num_jobs = models.IntegerField(blank=True, null=True)
    num_core_fs = models.IntegerField(blank=True, null=True)
    gene_len_left = models.FloatField(blank=True, null=True)
    kaks_left = models.FloatField(blank=True, null=True)
    meme_width = models.IntegerField(blank=True, null=True)
    meme_evalue = models.FloatField(blank=True, null=True)
    gene_len_full = models.FloatField(blank=True, null=True)
    kaks_full = models.FloatField(blank=True, null=True)
    num_plus_fs = models.IntegerField(blank=True, null=True)
    num_minus_fs = models.IntegerField(blank=True, null=True)
    subst_per_site = models.FloatField(blank=True, null=True)
    num_down_rbs = models.IntegerField(blank=True, null=True)
    down_rbs_score = models.FloatField(blank=True, null=True)
    meme_period_1 = models.FloatField(blank=True, null=True)
    meme_period_2 = models.FloatField(blank=True, null=True)
    meme_period_3 = models.FloatField(blank=True, null=True)
    meme_profile_sum = models.FloatField(blank=True, null=True)
    main_signal_id = models.IntegerField(blank=True, null=True)
    main_signal_fract = models.FloatField(blank=True, null=True)
    stop_stop_len = models.FloatField(blank=True, null=True)
    ss50_word_exp = models.FloatField(blank=True, null=True)
    ss50_word_var = models.FloatField(blank=True, null=True)
    ss50_identity = models.FloatField(blank=True, null=True)
    best_word_score = models.FloatField(blank=True, null=True)
    best_word_num_fs = models.IntegerField(blank=True, null=True)
    best_word = models.CharField(max_length=255, blank=True, null=True)
    best_p_word = models.CharField(max_length=255, blank=True, null=True)
    best_p_word_num_fs = models.IntegerField(blank=True, null=True)
    trf_num_fs = models.IntegerField(blank=True, null=True)
    trf_unit_num_fs = models.IntegerField(blank=True, null=True)
    trf_unit_seq = models.CharField(max_length=255, blank=True, null=True)
    stop1_start2_num_fs = models.IntegerField(blank=True, null=True)
    poly_at_num_fs = models.IntegerField(blank=True, null=True)
    user = models.ForeignKey('Users', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'cofs'


class ExtDbs(models.Model):
    db_id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=255)
    descr = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'ext_dbs'


class FeatureTags(models.Model):
    f = models.ForeignKey('Features', models.DO_NOTHING)
    tag = models.CharField(max_length=255, blank=True, null=True)
    value = models.CharField(max_length=10000, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'feature_tags'


class Features(models.Model):
    f_id = models.IntegerField(primary_key=True)
    job = models.ForeignKey('Jobs', models.DO_NOTHING)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()
    primary_tag = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'features'


class FsDbs(models.Model):
    fs_id = models.IntegerField()
    db_id = models.IntegerField()
    ext_id = models.CharField(max_length=255)

    class Meta:
        managed = False
        db_table = 'fs_dbs'


class GmGenes(models.Model):
    gene_id = models.IntegerField(primary_key=True)
    job_id = models.IntegerField(blank=True, null=True)
    gene = models.IntegerField(db_column='Gene', blank=True, null=True)  # Field name made lowercase.
    strand = models.CharField(db_column='Strand', max_length=1, blank=True, null=True)  # Field name made lowercase.
    gene_left = models.IntegerField(db_column='Gene_left', blank=True, null=True)  # Field name made lowercase.
    gene_right = models.IntegerField(db_column='Gene_right', blank=True, null=True)  # Field name made lowercase.
    len = models.IntegerField(blank=True, null=True)
    class_field = models.IntegerField(db_column='Class', blank=True, null=True)  # Field name made lowercase. Field renamed because it was a Python reserved word.
    spacer = models.IntegerField(db_column='Spacer', blank=True, null=True)  # Field name made lowercase.
    rbs_score = models.FloatField(db_column='RBS_score', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'gm_genes'


class GtFs(models.Model):
    fs_id = models.IntegerField(primary_key=True)
    job = models.ForeignKey('Jobs', models.DO_NOTHING, blank=True, null=True)
    fs_coord = models.IntegerField(db_column='FS_coord', blank=True, null=True)  # Field name made lowercase.
    strand = models.CharField(max_length=1, blank=True, null=True)
    gene_left = models.IntegerField(db_column='Gene_Left', blank=True, null=True)  # Field name made lowercase.
    gene_right = models.IntegerField(db_column='Gene_Right', blank=True, null=True)  # Field name made lowercase.
    gene_nc_len = models.IntegerField(db_column='Gene_NC_Len', blank=True, null=True)  # Field name made lowercase.
    fragment_left = models.IntegerField(db_column='Fragment_Left', blank=True, null=True)  # Field name made lowercase.
    fragment_right = models.IntegerField(db_column='Fragment_Right', blank=True, null=True)  # Field name made lowercase.
    fs_path_score = models.IntegerField(db_column='FS_Path_Score', blank=True, null=True)  # Field name made lowercase.
    wo_fs_path_score = models.IntegerField(db_column='Wo_FS_Path_score', blank=True, null=True)  # Field name made lowercase.
    fs_score = models.FloatField(db_column='FS_Score', blank=True, null=True)  # Field name made lowercase.
    lg_fs_score = models.FloatField(db_column='LG_FS_Score', blank=True, null=True)  # Field name made lowercase.
    filter = models.CharField(max_length=255, blank=True, null=True)
    gene_seq = models.TextField(blank=True, null=True)
    prot_seq = models.TextField(blank=True, null=True)
    prot_fs_coord = models.IntegerField(blank=True, null=True)
    type = models.CharField(max_length=2, blank=True, null=True)
    down_dist = models.IntegerField(blank=True, null=True)
    down_rbs = models.FloatField(blank=True, null=True)
    init_gene_seq = models.TextField(blank=True, null=True)
    stop_stop_gene_seq = models.TextField(blank=True, null=True)
    ss50_frame_seq = models.CharField(max_length=255, blank=True, null=True)
    ss50_frame_signal = models.CharField(max_length=255, blank=True, null=True)
    ss50_frame_weight = models.FloatField(blank=True, null=True)
    gene_gc = models.FloatField(blank=True, null=True)
    pseudogene_f = models.ForeignKey(Features, models.DO_NOTHING, blank=True, null=True)
    gbk_programmed = models.CharField(max_length=255, blank=True, null=True)
    down_struct_num = models.IntegerField(blank=True, null=True)
    down_struct_energy = models.FloatField(blank=True, null=True)
    ss_start = models.IntegerField(blank=True, null=True)
    relative_fs_coord = models.FloatField(blank=True, null=True)
    gbk_variable_prot = models.CharField(max_length=255, blank=True, null=True)
    trf_repeat_seq = models.CharField(max_length=255, blank=True, null=True)
    trf_unit_seq = models.CharField(max_length=255, blank=True, null=True)
    trf_copy_number = models.FloatField(blank=True, null=True)
    trf_repeat_gc = models.FloatField(blank=True, null=True)
    stop1_start2_dist = models.IntegerField(blank=True, null=True)
    poly_at_seq = models.CharField(max_length=255, blank=True, null=True)
    cai_orf1 = models.FloatField(blank=True, null=True)
    user = models.ForeignKey('Users', models.DO_NOTHING, blank=True, null=True)
    name = models.CharField(max_length=255, blank=True, null=True)
    descr = models.CharField(max_length=255, blank=True, null=True)
    ss_frame_seq = models.TextField(blank=True, null=True)
    exon_junct_dist = models.IntegerField(blank=True, null=True)
    num_exons = models.IntegerField(blank=True, null=True)
    is_as_fs = models.IntegerField(blank=True, null=True)
    origin = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gt_fs'


class GtfsBlastp(models.Model):
    fs_id = models.IntegerField()
    hit_id = models.CharField(max_length=255)
    identity = models.FloatField(blank=True, null=True)
    ali_len = models.IntegerField(blank=True, null=True)
    mismatches = models.IntegerField(blank=True, null=True)
    gap_openings = models.IntegerField(blank=True, null=True)
    q_start = models.IntegerField(blank=True, null=True)
    q_end = models.IntegerField(blank=True, null=True)
    h_start = models.IntegerField(blank=True, null=True)
    h_end = models.IntegerField(blank=True, null=True)
    evalue = models.CharField(max_length=255, blank=True, null=True)
    score = models.FloatField(blank=True, null=True)
    is_fs_hit = models.IntegerField(blank=True, null=True)
    fs_to_hit_edge_dist = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gtfs_blastp'


class GtfsBlastpEuk(models.Model):
    fs_id = models.IntegerField()
    hit_id = models.CharField(max_length=255)
    identity = models.FloatField(blank=True, null=True)
    ali_len = models.IntegerField(blank=True, null=True)
    mismatches = models.IntegerField(blank=True, null=True)
    gap_openings = models.IntegerField(blank=True, null=True)
    q_start = models.IntegerField(blank=True, null=True)
    q_end = models.IntegerField(blank=True, null=True)
    h_start = models.IntegerField(blank=True, null=True)
    h_end = models.IntegerField(blank=True, null=True)
    evalue = models.CharField(max_length=255, blank=True, null=True)
    score = models.FloatField(blank=True, null=True)
    is_fs_hit = models.IntegerField(blank=True, null=True)
    fs_to_hit_edge_dist = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gtfs_blastp_euk'


class GtfsExons(models.Model):
    qname = models.CharField(db_column='qName', max_length=255, blank=True, null=True)  # Field name made lowercase.
    qsize = models.IntegerField(db_column='qSize', blank=True, null=True)  # Field name made lowercase.
    qstart = models.IntegerField(db_column='qStart', blank=True, null=True)  # Field name made lowercase.
    qend = models.IntegerField(db_column='qEnd', blank=True, null=True)  # Field name made lowercase.
    tname = models.CharField(db_column='tName', max_length=255, blank=True, null=True)  # Field name made lowercase.
    tstart = models.IntegerField(db_column='tStart', blank=True, null=True)  # Field name made lowercase.
    tend = models.IntegerField(db_column='tEnd', blank=True, null=True)  # Field name made lowercase.
    blockcount = models.IntegerField(db_column='blockCount', blank=True, null=True)  # Field name made lowercase.
    blocksizes = models.CharField(db_column='blockSizes', max_length=255, blank=True, null=True)  # Field name made lowercase.
    qstarts = models.CharField(db_column='qStarts', max_length=255, blank=True, null=True)  # Field name made lowercase.
    tstarts = models.CharField(db_column='tStarts', max_length=255, blank=True, null=True)  # Field name made lowercase.
    mismatches = models.IntegerField(db_column='misMatches', blank=True, null=True)  # Field name made lowercase.
    src = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gtfs_exons'


class GtfsPfam(models.Model):
    fs = models.ForeignKey(GtFs, models.DO_NOTHING)
    pfam_name = models.CharField(max_length=255)
    pfam_acc = models.CharField(max_length=255)
    pfam_type = models.CharField(max_length=255, blank=True, null=True)
    pfam_class = models.CharField(max_length=255, blank=True, null=True)
    start = models.IntegerField(blank=True, null=True)
    end = models.IntegerField(blank=True, null=True)
    ali_start = models.IntegerField(blank=True, null=True)
    ali_end = models.IntegerField(blank=True, null=True)
    hmm_start = models.IntegerField(blank=True, null=True)
    hmm_end = models.IntegerField(blank=True, null=True)
    bitscore = models.FloatField(blank=True, null=True)
    evalue = models.CharField(max_length=255, blank=True, null=True)
    significant = models.IntegerField(blank=True, null=True)
    evidence = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gtfs_pfam'


class GtfsRpf(models.Model):
    fs = models.ForeignKey(GtFs, models.DO_NOTHING, blank=True, null=True)
    rpf_id = models.CharField(max_length=255, blank=True, null=True)
    coord = models.IntegerField(blank=True, null=True)
    mism = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gtfs_rpf'


class GtfsSignals(models.Model):
    fs_id = models.IntegerField()
    signal_id = models.IntegerField()
    motif = models.CharField(max_length=64)

    class Meta:
        managed = False
        db_table = 'gtfs_signals'
        unique_together = (('fs_id', 'signal_id'),)


class GtfsWords(models.Model):
    fs_id = models.IntegerField()
    word = models.CharField(max_length=255)
    coord = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gtfs_words'


class JobDbs(models.Model):
    job_id = models.IntegerField()
    db_id = models.IntegerField()
    ext_id = models.CharField(max_length=255)

    class Meta:
        managed = False
        db_table = 'job_dbs'


class Jobs(models.Model):
    job_id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=255, blank=True, null=True)
    descr = models.TextField(blank=True, null=True)
    path = models.CharField(max_length=255, blank=True, null=True)
    c_date = models.DateTimeField(blank=True, null=True)
    type = models.CharField(max_length=255, blank=True, null=True)
    status = models.CharField(max_length=255, blank=True, null=True)
    pbs_id = models.CharField(max_length=255, blank=True, null=True)
    user = models.ForeignKey('Users', models.DO_NOTHING, blank=True, null=True)
    email = models.CharField(max_length=255, blank=True, null=True)
    seq_gc = models.FloatField(blank=True, null=True)
    seq_len = models.IntegerField(blank=True, null=True)
    num_fs = models.IntegerField(blank=True, null=True)
    num_b_fs = models.IntegerField(blank=True, null=True)
    num_p_fs = models.IntegerField(blank=True, null=True)
    num_bp_fs = models.IntegerField(blank=True, null=True)
    num_r_fs = models.IntegerField(blank=True, null=True)
    num_bpr_fs = models.IntegerField(blank=True, null=True)
    kingdom = models.CharField(max_length=255, blank=True, null=True)
    species = models.CharField(max_length=255, blank=True, null=True)
    taxa = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'jobs'


class Labels(models.Model):
    label_id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=255)
    descr = models.CharField(max_length=255, blank=True, null=True)
    logo_type = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'labels'


class PashaHits(models.Model):
    q_cluster = models.CharField(max_length=255)
    q_frameshift = models.CharField(max_length=255)
    q_start = models.IntegerField()
    q_end = models.IntegerField()
    h_fs_id = models.IntegerField()
    h_start = models.IntegerField()
    h_end = models.IntegerField()
    h_fs_ali_coord = models.IntegerField(blank=True, null=True)
    evalue = models.CharField(max_length=255, blank=True, null=True)
    identity = models.FloatField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'pasha_hits'


class PfamWords(models.Model):
    pfam_name = models.CharField(max_length=255)
    word = models.CharField(max_length=255)
    obs_num = models.FloatField(blank=True, null=True)
    exp_num = models.FloatField(blank=True, null=True)
    weight = models.FloatField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'pfam_words'


class RecodeGtfs(models.Model):
    fs_id = models.IntegerField()
    recode_id_product = models.CharField(max_length=255, blank=True, null=True)
    fs_start = models.IntegerField(blank=True, null=True)
    fs_end = models.IntegerField(blank=True, null=True)
    recode_start = models.IntegerField(blank=True, null=True)
    recode_end = models.IntegerField(blank=True, null=True)
    evalue = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'recode_gtfs'
        unique_together = (('recode_id_product', 'fs_id'),)


class RecodeProducts(models.Model):
    id = models.IntegerField(primary_key=True)
    id_product = models.CharField(max_length=64)
    id_molecule = models.IntegerField()
    id_recode = models.CharField(max_length=64)
    name = models.CharField(max_length=64)
    sequence = models.TextField(blank=True, null=True)
    modification = models.TextField(blank=True, null=True)
    coordinates = models.TextField(blank=True, null=True)
    description = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'recode_products'


class RecodeRecoding(models.Model):
    id = models.IntegerField(primary_key=True)
    id_recode = models.CharField(max_length=64)
    id_product = models.CharField(max_length=64)
    event = models.IntegerField()
    experimental = models.CharField(max_length=1)
    position = models.BigIntegerField(blank=True, null=True)
    codon = models.CharField(max_length=64, blank=True, null=True)
    upstream = models.CharField(max_length=64, blank=True, null=True)
    esite = models.CharField(max_length=64, blank=True, null=True)
    asite = models.CharField(max_length=64, blank=True, null=True)
    psite = models.CharField(max_length=64, blank=True, null=True)
    downstream = models.CharField(max_length=64, blank=True, null=True)
    model = models.TextField(blank=True, null=True)
    description = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'recode_recoding'

"""
class Sessions(models.Model):
    id = models.CharField(unique=True, max_length=32)
    a_session = models.TextField()

    class Meta:
        managed = False
        db_table = 'sessions'
"""

class Signals(models.Model):
    signal_id = models.IntegerField(primary_key=True)
    pattern_re = models.CharField(max_length=64)
    min_len = models.IntegerField()
    mechanism = models.CharField(max_length=64, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'signals'


class Users(models.Model):
    id = models.IntegerField(primary_key=True)
    c_date = models.DateTimeField()
    name = models.CharField(unique=True, max_length=255)
    descr = models.TextField(blank=True, null=True)
    pass_field = models.CharField(db_column='pass', max_length=255, blank=True, null=True)  # Field renamed because it was a Python reserved word.
    email = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'users'


class ValidFsV(models.Model):
    type = models.CharField(max_length=2)
    job_id = models.IntegerField(blank=True, null=True)
    fs_id = models.IntegerField(unique=True)
    fs_coord = models.IntegerField(blank=True, null=True)
    prot_fs_coord = models.IntegerField(blank=True, null=True)
    strand = models.CharField(max_length=1, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'valid_fs_v'


class ValidationBlastp(models.Model):
    fs_id = models.IntegerField(blank=True, null=True)
    fs_coord = models.IntegerField(db_column='FS_coord', blank=True, null=True)  # Field name made lowercase.
    fs_prot_len = models.IntegerField(db_column='FS_prot_len', blank=True, null=True)  # Field name made lowercase.
    fs_coord_in_protein = models.IntegerField(db_column='FS_coord_in_protein', blank=True, null=True)  # Field name made lowercase.
    num_hits = models.IntegerField(db_column='Num_hits', blank=True, null=True)  # Field name made lowercase.
    fs_hit = models.CharField(db_column='FS_hit', max_length=255, blank=True, null=True)  # Field name made lowercase.
    fs_hit_start = models.IntegerField(db_column='FS_hit_start', blank=True, null=True)  # Field name made lowercase.
    fs_hit_end = models.IntegerField(db_column='FS_hit_end', blank=True, null=True)  # Field name made lowercase.
    fs_hit_evalue = models.FloatField(db_column='FS_hit_evalue', blank=True, null=True)  # Field name made lowercase.
    min_dist_from_fs_to_hit_edge = models.IntegerField(db_column='Min_dist_from_fs_to_hit_edge')  # Field name made lowercase.
    notes = models.TextField(db_column='Notes', blank=True, null=True)  # Field name made lowercase.
    rid = models.CharField(db_column='RID', max_length=255, blank=True, null=True)  # Field name made lowercase.
    fs_prot_seq = models.TextField(db_column='FS_prot_seq', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'validation_blastp'


class ValidationPfam(models.Model):
    fs_id = models.IntegerField(blank=True, null=True)
    fs_coord = models.IntegerField(db_column='FS_coord', blank=True, null=True)  # Field name made lowercase.
    fs_prot_len = models.IntegerField(db_column='FS_prot_len', blank=True, null=True)  # Field name made lowercase.
    fs_coord_in_protein = models.IntegerField(db_column='FS_coord_in_protein', blank=True, null=True)  # Field name made lowercase.
    total_num_domains = models.IntegerField(db_column='Total_num_domains', blank=True, null=True)  # Field name made lowercase.
    num_fs_domains = models.IntegerField(db_column='Num_fs_domains', blank=True, null=True)  # Field name made lowercase.
    fs_domain = models.CharField(db_column='FS_domain', max_length=255, blank=True, null=True)  # Field name made lowercase.
    fs_domain_start = models.IntegerField(db_column='FS_domain_start', blank=True, null=True)  # Field name made lowercase.
    fs_domain_end = models.IntegerField(db_column='FS_domain_end', blank=True, null=True)  # Field name made lowercase.
    fs_domain_ali_start = models.IntegerField(db_column='FS_domain_ali_start', blank=True, null=True)  # Field name made lowercase.
    fs_domain_ali_end = models.IntegerField(db_column='FS_domain_ali_end', blank=True, null=True)  # Field name made lowercase.
    fs_domain_hmm_start = models.IntegerField(db_column='FS_domain_hmm_start', blank=True, null=True)  # Field name made lowercase.
    fs_domain_hmm_end = models.IntegerField(db_column='FS_domain_hmm_end', blank=True, null=True)  # Field name made lowercase.
    fs_domain_bitscore = models.IntegerField(db_column='FS_domain_bitscore', blank=True, null=True)  # Field name made lowercase.
    fs_domain_evalue = models.FloatField(db_column='FS_domain_evalue', blank=True, null=True)  # Field name made lowercase.
    min_dist_fs_to_domain_edge = models.IntegerField(db_column='Min_dist_fs_to_domain_edge')  # Field name made lowercase.
    notes = models.TextField(db_column='Notes', blank=True, null=True)  # Field name made lowercase.
    job_id = models.CharField(db_column='Job_id', max_length=255, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'validation_pfam'

"""
class ZzzTempNn(models.Model):
    qid = models.IntegerField(db_column='QID')  # Field name made lowercase.
    id = models.IntegerField(db_column='ID', blank=True, null=True)  # Field name made lowercase.
    val = models.IntegerField(db_column='VAL', blank=True, null=True)  # Field name made lowercase.
    descr = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'zzz_temp_nn'


class ZzzTempNs(models.Model):
    qid = models.IntegerField(db_column='QID', blank=True, null=True)  # Field name made lowercase.
    id = models.IntegerField(db_column='ID', blank=True, null=True)  # Field name made lowercase.
    val = models.CharField(db_column='VAL', max_length=255, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'zzz_temp_ns'
"""
