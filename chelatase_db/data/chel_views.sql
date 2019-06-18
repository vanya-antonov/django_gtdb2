
drop table if exists chel_feats_v;
create table chel_feats_v as
select f.id, f.seq_id, s.org_id,
    f.type, f.name, f.descr, f.start, f.end, f.strand, f.parent_id,
    (select t.num from feat_params t where t.parent_id=f.id and t.name='translation') AS prot_len,
    (select t.value from feat_params t where t.parent_id=f.id and t.name='gene') AS gene,
    (select t.value from feat_params t where t.parent_id=f.id and t.name='chel_gene') AS chel_gene,
    (select t.value from feat_params t where t.parent_id=f.id and t.name='chel_subunit') AS chel_subunit,
    (select t.num   from feat_params t where t.parent_id=f.id and t.name='chel_evalue') AS chel_evalue,
    (select fs.len from feat_fshifts ff, fshifts fs where ff.feat_id=f.id and fs.id=ff.fshift_id limit 1) AS fs_len,
    (select fs.len from feats t, feat_fshifts ff, fshifts fs where t.parent_id=f.id and ff.feat_id=t.id and fs.id=ff.fshift_id limit 1) AS child_fs_len,
    (select count(1) from feats t where t.parent_id=f.id) AS num_kids,
    o.name AS org_name, o.phylum, o.kingdom, o.genus,
    (select t.data from feat_params t where t.parent_id=f.id and t.name='translation') AS translation,
    (select t.data from feat_params t where t.parent_id=f.id and t.name='seq_nt') AS seq_nt
from orgs o, seqs s, feats f
where s.org_id=o.id and f.seq_id=s.id;
CREATE INDEX chelfeatsv_id_i ON chel_feats_v(id);
CREATE INDEX chelfeatsv_seqid_i ON chel_feats_v(seq_id);
CREATE INDEX chelfeatsv_orgid_i ON chel_feats_v(org_id);


drop table if exists chel_orgs_v;
create table chel_orgs_v as
select o.id, o.name, o.genus, o.phylum, o.kingdom,
    (select t.value from org_params t where t.parent_id=o.id and t.name='taxonomy' and t.num=2) AS tax2,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_subunit='L') AS num_L,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_subunit='M') AS num_M,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_subunit='S') AS num_S,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_subunit='M' and t.fs_len is NULL) AS num_M_zero,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_subunit='M' and t.fs_len =-1) AS num_M_minus,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_subunit='M' and t.fs_len = 1) AS num_M_plus,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='cobN') AS num_cobN,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='chlH') AS num_chlH,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='bchH') AS num_bchH,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='cobT') AS num_cobT,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='chlD') AS num_chlD,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='bchD') AS num_bchD,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='cobS') AS num_cobS,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='chlI') AS num_chlI,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene='bchI') AS num_bchI,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='cobN') AS evalue_cobN,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='chlH') AS evalue_chlH,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='bchH') AS evalue_bchH,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='cobT') AS evalue_cobT,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='chlD') AS evalue_chlD,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='bchD') AS evalue_bchD,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='cobS') AS evalue_cobS,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='chlI') AS evalue_chlI,
    (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene='bchI') AS evalue_bchI
from orgs o
where exists (select 1 from chel_feats_v t where t.org_id=o.id and t.chel_subunit='M');
CREATE UNIQUE INDEX chelorgsv_id_i ON chel_orgs_v(id);

