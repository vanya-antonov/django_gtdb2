## search_TTA_genes.pl
Search TTA codons in sequences of *file.gbk*.
Creates a file (table) with statistics: **statistic_TTA-vs-FS_genes.tsv** or stdout output.

### USAGE

./search_TTA_genes.pl *file.gbk* [OPTIONS]

Here *file.gbk* - input GenBank file only

Example:
```bash
./search_TTA_genes.pl NC_010572.1.gbk -auto -output=stdout

./search_TTA_genes.pl NC_010572.1.gbk -auto
# or same
./search_TTA_genes.pl NC_010572.1.gbk -save_tta_fna NC_010572.1.TTA.fna -save_tta_faa NC_010572.1.TTA.faa -output NC_010572.1.TTA.tsv
```

### OUTPUT TABLE FORMAT:

| column | name                  | description                                        | example |
|:------:|:----------------------|:---------------------------------------------------|:--------|
| 1      | TAXON                 | NCBI taxon ID of organism                          | 227882  |
| 2      | ORG\_ID               | internal organism ID                               | 167     |
| 3      | ORG_NAME              | organism name                                      | Streptomyces avermitilis MA-4680 = NBRC 14893 |
| 4      | NUM_GENES             | total number of CDS genes for organism             | 8025    |
| 5      | NUM_TTA_GENES         | total number of genes with TTA codon (TTA-genes)   | 274     |
| 6      | NUM_FS_GENES          | total number of FS-genes. Empty for *--skip_fs* mode | 1188    |
| 7      | NUM_FS_and_TTA_GENES  | intersection of FS-genes with TTA-genes            | 29      |
| 8      | NUM_COFS              | total number of clusters that include FS-genes     | 293     |
| 9      | NUM_FS_GENES_in_COFS  | total number of FS-genes in clusters               | 300     |
| 10     | NUM_TTA_GENES_in_COFS | total number of TTA-genes that are similar to FS-genes from ALL clusters | 124 |
| 11     | NUM_FS_and_TTA_GENES_in_COFS | intersection of FS-genes with TTA-genes in clusters | 11 |
| 12     | ACC_IDs               | Accession ID(;s) of locus/genomic sequence(s)      | NC_003155.5;NC_004719.1 |
| 13     | COF_IDs               | List of **Cluster_ID=number genes** with TTA codon | 1004907=3;1004909=1;... |
| 14     | FS_IDs                | List of **FS-TTA-gene_ID=gid1,gid2,...**           | 74297=NC_003155.5:p25699.4695.3,NC_003155.5:m28699.1078.3;... |
| 15     | WOFS\_IDs             | List of **Cluster_ID=gid1,gid2,...** without FS-genes | 1004907=NC_003155.5:m5172533.817.1,NC_003155.5:p9004239.817.1;... |
| 16     | GENE\_IDs             | List of all CDS genes as **gid:locus_tag**         | NC_003155.5:m1002287.1724.0:SAVERM_RS04685;... |

:star: **FS-TTA-gene_ID** is internal gene ID for FS-gene with TTA-codon, e.g. 74297

:star: **Cluster_ID** is internal cluster ID, e.g. 1004907

:star: **gid** is **acc_id:\<strand\>start.length.gtag**

:star: **acc_id** is Accession ID of sequence, e.g. NC_003155.5

:star: **strand** is **p** (positive) for **(+)**, or **n** (negative) for **(-)** strand.

:star: **start** is 1-based start location of the CDS-gene (*sloc*).

:star: **length** is length of the CDS-gene. This is mostly a positive integer. Sometimes it can be negative (*length=eloc-start*).

:star: **gtag** is 1 for TTA-gene, 2 - FS-gene, 3 - (TTA+FS)-gene, 0 - ordinary gene.

:star: **locus_tag** is name of LOCUS for the gene, e.g. SAVERM_RS04685

:exclamation: Fields ( ORG_ID, NUM_FS_GENES, NUM_FS_and_TTA_GENES, NUM_COFS, NUM_FS_GENES_in_COFS,
NUM_TTA_GENES_in_COFS, NUM_FS_and_TTA_GENES_in_COFS, COF_IDs, FS_IDs, WOFS_IDs ) are (empty | 0) for *--skip_fs* option.

#### NOTES (without *--skip_fs* option):
1. a configuration DB file *django_gtdb2/django/mysite/local_settings.json* (or specified by *--cfg_db*) is required.
2. BLAST program is required.


## add_TTA_genes_in_FEATS.pl
Fills and/or updates **feats**, and **feat_params** tables of **GTDB2** database (default)
with statistics from *statistic_TTA-vs-FS_genes.tsv* file.

### USAGE

./add_TTA_genes_in_FEATS.pl *[statistic.tsv]* [OPTIONS]

Here *statistic.tsv* - By default, *statistic_TTA-vs-FS_genes.tsv*.

Example:
```bash
./add_TTA_genes_in_FEATS.pl
# or
./add_TTA_genes_in_FEATS.pl  statistic_TTA-vs-FS_genes.tsv
```

## add_statistic_TTA-FS_in_ORG_PARAMS.pl
Fills and/or updates **org_params** table of **GTDB2** database (default) with statistics
from the *statistic.tsv* file

### USAGE

./add_statistic_TTA-FS_in_ORG_PARAMS.pl *statistic.tsv* [OPTIONS]

Here *statistic.tsv* - By default, for Bacteria is *statistic_TTA-vs-FS_genes.tsv*,
for Viruses|Phages is *statistic_Streptomyces_phages_TTA-vs-FS_genes.tsv*.

Example:
```bash
./add_statistic_TTA-FS_in_ORG_PARAMS.pl  statistic_Streptomyces_phages_TTA-vs-FS_genes.tsv -taxonomy Phage
```

## get_hosts_from_GenBank.pl

