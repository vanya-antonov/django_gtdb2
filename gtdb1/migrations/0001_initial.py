# Generated by Django 2.1.3 on 2018-12-07 15:17

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='CofGtfs',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('is_core', models.IntegerField()),
            ],
            options={
                'db_table': 'cof_gtfs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='CofHits',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('q_fs_id', models.IntegerField()),
                ('q_start', models.IntegerField()),
                ('q_end', models.IntegerField()),
                ('q_fs_ali_coord', models.IntegerField()),
                ('h_fs_id', models.IntegerField()),
                ('h_start', models.IntegerField()),
                ('h_end', models.IntegerField()),
                ('h_fs_ali_coord', models.IntegerField()),
                ('ali_len', models.IntegerField()),
                ('score', models.IntegerField()),
                ('evalue', models.CharField(max_length=255)),
                ('identity', models.FloatField()),
            ],
            options={
                'db_table': 'cof_hits',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='CofLabels',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('label_id', models.IntegerField()),
            ],
            options={
                'db_table': 'cof_labels',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Cofs',
            fields=[
                ('cof_id', models.IntegerField(primary_key=True, serialize=False)),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('c_date', models.DateTimeField()),
                ('num_b_fs', models.IntegerField(blank=True, null=True)),
                ('num_p_fs', models.IntegerField(blank=True, null=True)),
                ('num_r_fs', models.IntegerField(blank=True, null=True)),
                ('num_bp_fs', models.IntegerField(blank=True, null=True)),
                ('dir', models.CharField(blank=True, max_length=255, null=True)),
                ('has_paml', models.IntegerField(blank=True, null=True)),
                ('num_fs', models.IntegerField(blank=True, null=True)),
                ('num_genus', models.IntegerField(blank=True, null=True)),
                ('kingdoms', models.CharField(blank=True, max_length=32, null=True)),
                ('fs_dna_num_left_gaps', models.IntegerField(blank=True, null=True)),
                ('fs_dna_num_right_gaps', models.IntegerField(blank=True, null=True)),
                ('fs_dna_num_inner_gaps', models.IntegerField(blank=True, null=True)),
                ('fs_dna_num_identical_cols', models.IntegerField(blank=True, null=True)),
                ('fs_dna_identical_block', models.CharField(blank=True, max_length=255, null=True)),
                ('fs_dna_trimmed_len', models.IntegerField(blank=True, null=True)),
                ('num_species', models.IntegerField(blank=True, null=True)),
                ('num_jobs', models.IntegerField(blank=True, null=True)),
                ('num_core_fs', models.IntegerField(blank=True, null=True)),
                ('gene_len_left', models.FloatField(blank=True, null=True)),
                ('kaks_left', models.FloatField(blank=True, null=True)),
                ('meme_width', models.IntegerField(blank=True, null=True)),
                ('meme_evalue', models.FloatField(blank=True, null=True)),
                ('gene_len_full', models.FloatField(blank=True, null=True)),
                ('kaks_full', models.FloatField(blank=True, null=True)),
                ('num_plus_fs', models.IntegerField(blank=True, null=True)),
                ('num_minus_fs', models.IntegerField(blank=True, null=True)),
                ('subst_per_site', models.FloatField(blank=True, null=True)),
                ('num_down_rbs', models.IntegerField(blank=True, null=True)),
                ('down_rbs_score', models.FloatField(blank=True, null=True)),
                ('meme_period_1', models.FloatField(blank=True, null=True)),
                ('meme_period_2', models.FloatField(blank=True, null=True)),
                ('meme_period_3', models.FloatField(blank=True, null=True)),
                ('meme_profile_sum', models.FloatField(blank=True, null=True)),
                ('main_signal_id', models.IntegerField(blank=True, null=True)),
                ('main_signal_fract', models.FloatField(blank=True, null=True)),
                ('stop_stop_len', models.FloatField(blank=True, null=True)),
                ('ss50_word_exp', models.FloatField(blank=True, null=True)),
                ('ss50_word_var', models.FloatField(blank=True, null=True)),
                ('ss50_identity', models.FloatField(blank=True, null=True)),
                ('best_word_score', models.FloatField(blank=True, null=True)),
                ('best_word_num_fs', models.IntegerField(blank=True, null=True)),
                ('best_word', models.CharField(blank=True, max_length=255, null=True)),
                ('best_p_word', models.CharField(blank=True, max_length=255, null=True)),
                ('best_p_word_num_fs', models.IntegerField(blank=True, null=True)),
                ('trf_num_fs', models.IntegerField(blank=True, null=True)),
                ('trf_unit_num_fs', models.IntegerField(blank=True, null=True)),
                ('trf_unit_seq', models.CharField(blank=True, max_length=255, null=True)),
                ('stop1_start2_num_fs', models.IntegerField(blank=True, null=True)),
                ('poly_at_num_fs', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'cofs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='CofWords',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('cof_id', models.CharField(max_length=255)),
                ('word', models.CharField(max_length=255)),
                ('obs_num', models.FloatField(blank=True, null=True)),
                ('weight', models.FloatField(blank=True, null=True)),
            ],
            options={
                'db_table': 'cof_words',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='ExtDbs',
            fields=[
                ('db_id', models.IntegerField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=255)),
                ('descr', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'ext_dbs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Features',
            fields=[
                ('f_id', models.IntegerField(primary_key=True, serialize=False)),
                ('start', models.IntegerField()),
                ('end', models.IntegerField()),
                ('strand', models.IntegerField()),
                ('primary_tag', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'features',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='FeatureTags',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tag', models.CharField(blank=True, max_length=255, null=True)),
                ('value', models.CharField(blank=True, max_length=10000, null=True)),
            ],
            options={
                'db_table': 'feature_tags',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='FsDbs',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fs_id', models.IntegerField()),
                ('db_id', models.IntegerField()),
                ('ext_id', models.CharField(max_length=255)),
            ],
            options={
                'db_table': 'fs_dbs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GmGenes',
            fields=[
                ('gene_id', models.IntegerField(primary_key=True, serialize=False)),
                ('job_id', models.IntegerField(blank=True, null=True)),
                ('gene', models.IntegerField(blank=True, db_column='Gene', null=True)),
                ('strand', models.CharField(blank=True, db_column='Strand', max_length=1, null=True)),
                ('gene_left', models.IntegerField(blank=True, db_column='Gene_left', null=True)),
                ('gene_right', models.IntegerField(blank=True, db_column='Gene_right', null=True)),
                ('len', models.IntegerField(blank=True, null=True)),
                ('class_field', models.IntegerField(blank=True, db_column='Class', null=True)),
                ('spacer', models.IntegerField(blank=True, db_column='Spacer', null=True)),
                ('rbs_score', models.FloatField(blank=True, db_column='RBS_score', null=True)),
            ],
            options={
                'db_table': 'gm_genes',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GtFs',
            fields=[
                ('fs_id', models.IntegerField(primary_key=True, serialize=False)),
                ('fs_coord', models.IntegerField(blank=True, db_column='FS_coord', null=True)),
                ('strand', models.CharField(blank=True, max_length=1, null=True)),
                ('gene_left', models.IntegerField(blank=True, db_column='Gene_Left', null=True)),
                ('gene_right', models.IntegerField(blank=True, db_column='Gene_Right', null=True)),
                ('gene_nc_len', models.IntegerField(blank=True, db_column='Gene_NC_Len', null=True)),
                ('fragment_left', models.IntegerField(blank=True, db_column='Fragment_Left', null=True)),
                ('fragment_right', models.IntegerField(blank=True, db_column='Fragment_Right', null=True)),
                ('fs_path_score', models.IntegerField(blank=True, db_column='FS_Path_Score', null=True)),
                ('wo_fs_path_score', models.IntegerField(blank=True, db_column='Wo_FS_Path_score', null=True)),
                ('fs_score', models.FloatField(blank=True, db_column='FS_Score', null=True)),
                ('lg_fs_score', models.FloatField(blank=True, db_column='LG_FS_Score', null=True)),
                ('filter', models.CharField(blank=True, max_length=255, null=True)),
                ('gene_seq', models.TextField(blank=True, null=True)),
                ('prot_seq', models.TextField(blank=True, null=True)),
                ('prot_fs_coord', models.IntegerField(blank=True, null=True)),
                ('type', models.CharField(blank=True, max_length=2, null=True)),
                ('down_dist', models.IntegerField(blank=True, null=True)),
                ('down_rbs', models.FloatField(blank=True, null=True)),
                ('init_gene_seq', models.TextField(blank=True, null=True)),
                ('stop_stop_gene_seq', models.TextField(blank=True, null=True)),
                ('ss50_frame_seq', models.CharField(blank=True, max_length=255, null=True)),
                ('ss50_frame_signal', models.CharField(blank=True, max_length=255, null=True)),
                ('ss50_frame_weight', models.FloatField(blank=True, null=True)),
                ('gene_gc', models.FloatField(blank=True, null=True)),
                ('gbk_programmed', models.CharField(blank=True, max_length=255, null=True)),
                ('down_struct_num', models.IntegerField(blank=True, null=True)),
                ('down_struct_energy', models.FloatField(blank=True, null=True)),
                ('ss_start', models.IntegerField(blank=True, null=True)),
                ('relative_fs_coord', models.FloatField(blank=True, null=True)),
                ('gbk_variable_prot', models.CharField(blank=True, max_length=255, null=True)),
                ('trf_repeat_seq', models.CharField(blank=True, max_length=255, null=True)),
                ('trf_unit_seq', models.CharField(blank=True, max_length=255, null=True)),
                ('trf_copy_number', models.FloatField(blank=True, null=True)),
                ('trf_repeat_gc', models.FloatField(blank=True, null=True)),
                ('stop1_start2_dist', models.IntegerField(blank=True, null=True)),
                ('poly_at_seq', models.CharField(blank=True, max_length=255, null=True)),
                ('cai_orf1', models.FloatField(blank=True, null=True)),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('descr', models.CharField(blank=True, max_length=255, null=True)),
                ('ss_frame_seq', models.TextField(blank=True, null=True)),
                ('exon_junct_dist', models.IntegerField(blank=True, null=True)),
                ('num_exons', models.IntegerField(blank=True, null=True)),
                ('is_as_fs', models.IntegerField(blank=True, null=True)),
                ('origin', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'gt_fs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GtfsBlastp',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fs_id', models.IntegerField()),
                ('hit_id', models.CharField(max_length=255)),
                ('identity', models.FloatField(blank=True, null=True)),
                ('ali_len', models.IntegerField(blank=True, null=True)),
                ('mismatches', models.IntegerField(blank=True, null=True)),
                ('gap_openings', models.IntegerField(blank=True, null=True)),
                ('q_start', models.IntegerField(blank=True, null=True)),
                ('q_end', models.IntegerField(blank=True, null=True)),
                ('h_start', models.IntegerField(blank=True, null=True)),
                ('h_end', models.IntegerField(blank=True, null=True)),
                ('evalue', models.CharField(blank=True, max_length=255, null=True)),
                ('score', models.FloatField(blank=True, null=True)),
                ('is_fs_hit', models.IntegerField(blank=True, null=True)),
                ('fs_to_hit_edge_dist', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'gtfs_blastp',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GtfsBlastpEuk',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fs_id', models.IntegerField()),
                ('hit_id', models.CharField(max_length=255)),
                ('identity', models.FloatField(blank=True, null=True)),
                ('ali_len', models.IntegerField(blank=True, null=True)),
                ('mismatches', models.IntegerField(blank=True, null=True)),
                ('gap_openings', models.IntegerField(blank=True, null=True)),
                ('q_start', models.IntegerField(blank=True, null=True)),
                ('q_end', models.IntegerField(blank=True, null=True)),
                ('h_start', models.IntegerField(blank=True, null=True)),
                ('h_end', models.IntegerField(blank=True, null=True)),
                ('evalue', models.CharField(blank=True, max_length=255, null=True)),
                ('score', models.FloatField(blank=True, null=True)),
                ('is_fs_hit', models.IntegerField(blank=True, null=True)),
                ('fs_to_hit_edge_dist', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'gtfs_blastp_euk',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GtfsExons',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('qname', models.CharField(blank=True, db_column='qName', max_length=255, null=True)),
                ('qsize', models.IntegerField(blank=True, db_column='qSize', null=True)),
                ('qstart', models.IntegerField(blank=True, db_column='qStart', null=True)),
                ('qend', models.IntegerField(blank=True, db_column='qEnd', null=True)),
                ('tname', models.CharField(blank=True, db_column='tName', max_length=255, null=True)),
                ('tstart', models.IntegerField(blank=True, db_column='tStart', null=True)),
                ('tend', models.IntegerField(blank=True, db_column='tEnd', null=True)),
                ('blockcount', models.IntegerField(blank=True, db_column='blockCount', null=True)),
                ('blocksizes', models.CharField(blank=True, db_column='blockSizes', max_length=255, null=True)),
                ('qstarts', models.CharField(blank=True, db_column='qStarts', max_length=255, null=True)),
                ('tstarts', models.CharField(blank=True, db_column='tStarts', max_length=255, null=True)),
                ('mismatches', models.IntegerField(blank=True, db_column='misMatches', null=True)),
                ('src', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'gtfs_exons',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GtfsPfam',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pfam_name', models.CharField(max_length=255)),
                ('pfam_acc', models.CharField(max_length=255)),
                ('pfam_type', models.CharField(blank=True, max_length=255, null=True)),
                ('pfam_class', models.CharField(blank=True, max_length=255, null=True)),
                ('start', models.IntegerField(blank=True, null=True)),
                ('end', models.IntegerField(blank=True, null=True)),
                ('ali_start', models.IntegerField(blank=True, null=True)),
                ('ali_end', models.IntegerField(blank=True, null=True)),
                ('hmm_start', models.IntegerField(blank=True, null=True)),
                ('hmm_end', models.IntegerField(blank=True, null=True)),
                ('bitscore', models.FloatField(blank=True, null=True)),
                ('evalue', models.CharField(blank=True, max_length=255, null=True)),
                ('significant', models.IntegerField(blank=True, null=True)),
                ('evidence', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'gtfs_pfam',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GtfsRpf',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('rpf_id', models.CharField(blank=True, max_length=255, null=True)),
                ('coord', models.IntegerField(blank=True, null=True)),
                ('mism', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'gtfs_rpf',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GtfsSignals',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fs_id', models.IntegerField()),
                ('signal_id', models.IntegerField()),
                ('motif', models.CharField(max_length=64)),
            ],
            options={
                'db_table': 'gtfs_signals',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='GtfsWords',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fs_id', models.IntegerField()),
                ('word', models.CharField(max_length=255)),
                ('coord', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'gtfs_words',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='JobDbs',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('job_id', models.IntegerField()),
                ('db_id', models.IntegerField()),
                ('ext_id', models.CharField(max_length=255)),
            ],
            options={
                'db_table': 'job_dbs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Jobs',
            fields=[
                ('job_id', models.IntegerField(primary_key=True, serialize=False)),
                ('name', models.CharField(blank=True, max_length=255, null=True)),
                ('descr', models.TextField(blank=True, null=True)),
                ('path', models.CharField(blank=True, max_length=255, null=True)),
                ('c_date', models.DateTimeField(blank=True, null=True)),
                ('type', models.CharField(blank=True, max_length=255, null=True)),
                ('status', models.CharField(blank=True, max_length=255, null=True)),
                ('pbs_id', models.CharField(blank=True, max_length=255, null=True)),
                ('email', models.CharField(blank=True, max_length=255, null=True)),
                ('seq_gc', models.FloatField(blank=True, null=True)),
                ('seq_len', models.IntegerField(blank=True, null=True)),
                ('num_fs', models.IntegerField(blank=True, null=True)),
                ('num_b_fs', models.IntegerField(blank=True, null=True)),
                ('num_p_fs', models.IntegerField(blank=True, null=True)),
                ('num_bp_fs', models.IntegerField(blank=True, null=True)),
                ('num_r_fs', models.IntegerField(blank=True, null=True)),
                ('num_bpr_fs', models.IntegerField(blank=True, null=True)),
                ('kingdom', models.CharField(blank=True, max_length=255, null=True)),
                ('species', models.CharField(blank=True, max_length=255, null=True)),
                ('taxa', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'jobs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Labels',
            fields=[
                ('label_id', models.IntegerField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=255)),
                ('descr', models.CharField(blank=True, max_length=255, null=True)),
                ('logo_type', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'labels',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='PashaHits',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('q_cluster', models.CharField(max_length=255)),
                ('q_frameshift', models.CharField(max_length=255)),
                ('q_start', models.IntegerField()),
                ('q_end', models.IntegerField()),
                ('h_fs_id', models.IntegerField()),
                ('h_start', models.IntegerField()),
                ('h_end', models.IntegerField()),
                ('h_fs_ali_coord', models.IntegerField(blank=True, null=True)),
                ('evalue', models.CharField(blank=True, max_length=255, null=True)),
                ('identity', models.FloatField(blank=True, null=True)),
            ],
            options={
                'db_table': 'pasha_hits',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='PfamWords',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pfam_name', models.CharField(max_length=255)),
                ('word', models.CharField(max_length=255)),
                ('obs_num', models.FloatField(blank=True, null=True)),
                ('exp_num', models.FloatField(blank=True, null=True)),
                ('weight', models.FloatField(blank=True, null=True)),
            ],
            options={
                'db_table': 'pfam_words',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='RecodeGtfs',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fs_id', models.IntegerField()),
                ('recode_id_product', models.CharField(blank=True, max_length=255, null=True)),
                ('fs_start', models.IntegerField(blank=True, null=True)),
                ('fs_end', models.IntegerField(blank=True, null=True)),
                ('recode_start', models.IntegerField(blank=True, null=True)),
                ('recode_end', models.IntegerField(blank=True, null=True)),
                ('evalue', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'recode_gtfs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='RecodeProducts',
            fields=[
                ('id', models.IntegerField(primary_key=True, serialize=False)),
                ('id_product', models.CharField(max_length=64)),
                ('id_molecule', models.IntegerField()),
                ('id_recode', models.CharField(max_length=64)),
                ('name', models.CharField(max_length=64)),
                ('sequence', models.TextField(blank=True, null=True)),
                ('modification', models.TextField(blank=True, null=True)),
                ('coordinates', models.TextField(blank=True, null=True)),
                ('description', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'recode_products',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='RecodeRecoding',
            fields=[
                ('id', models.IntegerField(primary_key=True, serialize=False)),
                ('id_recode', models.CharField(max_length=64)),
                ('id_product', models.CharField(max_length=64)),
                ('event', models.IntegerField()),
                ('experimental', models.CharField(max_length=1)),
                ('position', models.BigIntegerField(blank=True, null=True)),
                ('codon', models.CharField(blank=True, max_length=64, null=True)),
                ('upstream', models.CharField(blank=True, max_length=64, null=True)),
                ('esite', models.CharField(blank=True, max_length=64, null=True)),
                ('asite', models.CharField(blank=True, max_length=64, null=True)),
                ('psite', models.CharField(blank=True, max_length=64, null=True)),
                ('downstream', models.CharField(blank=True, max_length=64, null=True)),
                ('model', models.TextField(blank=True, null=True)),
                ('description', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'recode_recoding',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Signals',
            fields=[
                ('signal_id', models.IntegerField(primary_key=True, serialize=False)),
                ('pattern_re', models.CharField(max_length=64)),
                ('min_len', models.IntegerField()),
                ('mechanism', models.CharField(blank=True, max_length=64, null=True)),
            ],
            options={
                'db_table': 'signals',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Users',
            fields=[
                ('id', models.IntegerField(primary_key=True, serialize=False)),
                ('c_date', models.DateTimeField()),
                ('name', models.CharField(max_length=255, unique=True)),
                ('descr', models.TextField(blank=True, null=True)),
                ('pass_field', models.CharField(blank=True, db_column='pass', max_length=255, null=True)),
                ('email', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'users',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='ValidationBlastp',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fs_id', models.IntegerField(blank=True, null=True)),
                ('fs_coord', models.IntegerField(blank=True, db_column='FS_coord', null=True)),
                ('fs_prot_len', models.IntegerField(blank=True, db_column='FS_prot_len', null=True)),
                ('fs_coord_in_protein', models.IntegerField(blank=True, db_column='FS_coord_in_protein', null=True)),
                ('num_hits', models.IntegerField(blank=True, db_column='Num_hits', null=True)),
                ('fs_hit', models.CharField(blank=True, db_column='FS_hit', max_length=255, null=True)),
                ('fs_hit_start', models.IntegerField(blank=True, db_column='FS_hit_start', null=True)),
                ('fs_hit_end', models.IntegerField(blank=True, db_column='FS_hit_end', null=True)),
                ('fs_hit_evalue', models.FloatField(blank=True, db_column='FS_hit_evalue', null=True)),
                ('min_dist_from_fs_to_hit_edge', models.IntegerField(db_column='Min_dist_from_fs_to_hit_edge')),
                ('notes', models.TextField(blank=True, db_column='Notes', null=True)),
                ('rid', models.CharField(blank=True, db_column='RID', max_length=255, null=True)),
                ('fs_prot_seq', models.TextField(blank=True, db_column='FS_prot_seq', null=True)),
            ],
            options={
                'db_table': 'validation_blastp',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='ValidationPfam',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fs_id', models.IntegerField(blank=True, null=True)),
                ('fs_coord', models.IntegerField(blank=True, db_column='FS_coord', null=True)),
                ('fs_prot_len', models.IntegerField(blank=True, db_column='FS_prot_len', null=True)),
                ('fs_coord_in_protein', models.IntegerField(blank=True, db_column='FS_coord_in_protein', null=True)),
                ('total_num_domains', models.IntegerField(blank=True, db_column='Total_num_domains', null=True)),
                ('num_fs_domains', models.IntegerField(blank=True, db_column='Num_fs_domains', null=True)),
                ('fs_domain', models.CharField(blank=True, db_column='FS_domain', max_length=255, null=True)),
                ('fs_domain_start', models.IntegerField(blank=True, db_column='FS_domain_start', null=True)),
                ('fs_domain_end', models.IntegerField(blank=True, db_column='FS_domain_end', null=True)),
                ('fs_domain_ali_start', models.IntegerField(blank=True, db_column='FS_domain_ali_start', null=True)),
                ('fs_domain_ali_end', models.IntegerField(blank=True, db_column='FS_domain_ali_end', null=True)),
                ('fs_domain_hmm_start', models.IntegerField(blank=True, db_column='FS_domain_hmm_start', null=True)),
                ('fs_domain_hmm_end', models.IntegerField(blank=True, db_column='FS_domain_hmm_end', null=True)),
                ('fs_domain_bitscore', models.IntegerField(blank=True, db_column='FS_domain_bitscore', null=True)),
                ('fs_domain_evalue', models.FloatField(blank=True, db_column='FS_domain_evalue', null=True)),
                ('min_dist_fs_to_domain_edge', models.IntegerField(db_column='Min_dist_fs_to_domain_edge')),
                ('notes', models.TextField(blank=True, db_column='Notes', null=True)),
                ('job_id', models.CharField(blank=True, db_column='Job_id', max_length=255, null=True)),
            ],
            options={
                'db_table': 'validation_pfam',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='ValidFsV',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('type', models.CharField(max_length=2)),
                ('job_id', models.IntegerField(blank=True, null=True)),
                ('fs_id', models.IntegerField(unique=True)),
                ('fs_coord', models.IntegerField(blank=True, null=True)),
                ('prot_fs_coord', models.IntegerField(blank=True, null=True)),
                ('strand', models.CharField(blank=True, max_length=1, null=True)),
            ],
            options={
                'db_table': 'valid_fs_v',
                'managed': False,
            },
        ),
    ]
