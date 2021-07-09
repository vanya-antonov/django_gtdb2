#!/usr/bin/perl -w

use strict;
use warnings;

use File::Spec;
use Getopt::Long;
use DBI;

# Custom libraries
use MyLibGT::DBapp qw( :all );

my $VERSION = '1.04';

###
# Default Options
my $INFILE = 'statistic_TTA-vs-FS_genes.tsv';
my $DJANGO = '/home/alessandro/src/D/django_gtdb2/django';
my $CFG_DB = "$DJANGO/mysite/local_settings.json";
my $SESSION = 'default'; # 'gtdb2_cof'

###
# Parse input data
GetOptions(
	'infile=s'   => \$INFILE,
	'cfg_db=s'   => \$CFG_DB,
	'session=s'  => \$SESSION,
) or die &usage();

my $EMSG = "\x1b[31mERROR\x1b[0m";

# Read DataBase configurations
my( $dbh, $gtdb_dir ) = read_configDB( $SESSION, $CFG_DB );

my $save_db = SaveDB;

my $i;
open INFILE, $INFILE or die "Can't found $INFILE: $!\n";
while(<INFILE>){
	s/^\s+|\s+$//g;
	next if /^$/ || /^#/;
	next unless /^\d/;
=comment
# TAXONOMY=Bacteria
TAXON	ORG_ID	ORG_NAME	NUM_GENES	NUM_TTA_GENES	NUM_FS_GENES	NUM_FS_and_TTA_GENES	NUM_COFS	NUM_FS_GENES_in_COFS	NUM_TTA_GENES_in_COFS	NUM_FS_and_TTA_GENES_in_COFS	ACC_IDs	COF_IDs	FS_IDs	WOFS_IDs	GENE_IDs
1116232	1170	Streptomyces acidiscabies 84-104	9900	335	1717	20	230	230	147	7	NZ_AHBF01000244.1;NZ_AHBF01000243.1;...	1000540=1;1000572=1;...	11985=NZ_AHBF01000003.1:p97179.1633.3;12041=NZ_AHBF01000006.1:m116982.1976.3;...	1000540=NZ_AHBF01000158.1:m14454.629.1;1000572=NZ_AHBF01000060.1:m39197.1184.1;...	NZ_AHBF01000244.1:p44.992.0:ON12_RS48620;NZ_AHBF01000243.1:p118.413.0:ON12_RS52005;...
=cut
	my( $TAXON, $ORG_ID, $ORG_NAME, $NUM_GENES, $NUM_TTA_GENES, $NUM_FS_GENES, $NUM_FS_and_TTA_GENES,
		$NUM_COFS, $NUM_FS_GENES_in_COFS, $NUM_TTA_GENES_in_COFS, $NUM_FS_and_TTA_GENES_in_COFS,
		$ACC_IDs, $COF_IDs, $FS_IDs, $WOFS_IDs, $GENE_IDs ) = split /\t/;

	print ++$i, ") Prepare $ORG_ID - $ORG_NAME\n";

	my %ag;	# for all genes
	for my $g ( split ';', $GENE_IDs ){	# NZ_AHBF01000001.1:p95458.461.0:ON12_RS00405;NZ_AHBF01000001.1:p99943.1112.0:ON12_RS00420;...
		my( $gid, $locus_tag ) = $g=~/^(.*?):([^:]+)$/;

		my( $acc, $strand, $sloc, $eloc, $gtag, $f_TTA ) = parse_gid( $gid );
		next unless $f_TTA;	# discard non-TTA genes/locuses

		$ag{ $gid } = $locus_tag;
	}

	next unless $save_db;

	# Add TTA with FS-genes
	&create_TTA_records( $ORG_ID, $FS_IDs, 'FS_ID', $dbh, $DJANGO, \%ag );

	# Add TTA without FS-genes
	&create_TTA_records( $ORG_ID, $WOFS_IDs, 'CLUSTER', $dbh, $DJANGO, \%ag );

	# Add the rest
	if( scalar( keys %ag ) ){
		my $rest = '0=' . join(',', keys %ag );
		&create_TTA_records( $ORG_ID, $rest, '', $dbh, $DJANGO, \%ag );
	}

}
close INFILE;

exit;


sub create_TTA_records {
	my( $ORG_ID, $data_ids, $data_msg, $dbh, $DJANGO, $ag ) = @_;

	for my $d ( split ';', $data_ids ){	# 13441=NZ_AHBF01000125.1:p110537.1148.3,NZ_AHBF01000125.1:p109080.1901.3;13538=NZ_AHBF01000147.1:p11590.1100.3;...
		my( $msg_id, $gene_ids ) = split '=', $d; # 13441=NZ_AHBF01000125.1:p110537.1148.3,NZ_AHBF01000125.1:p109080.1901.3

		for my $gid ( split ',', $gene_ids ){ # NZ_AHBF01000125.1:p110537.1148.3,NZ_AHBF01000125.1:p109080.1901.3

			next unless exists $ag->{ $gid };

			my $locus = $ag->{ $gid };
			delete $ag->{ $gid };	# so as not to process again later

			my( $acc, $strand, $sloc, $eloc, $gtag, $f_TTA ) = parse_gid( $gid );	# in --> NZ_AHBF01000125.1:p110537.1148.3
			next unless $f_TTA;	# non-TTA-gene

			# Создание для локуса(gene) записей в таблицах seqs, feats, feat_params
			`python3 $DJANGO/manage.py  run_method_on_object  get_or_create_feat_from_locus_tag  Org $ORG_ID  $locus  CDS`;

			my $TTAs = &search_TTA_codons( $dbh, $locus, $strand, $sloc, $eloc );
			next unless $TTAs;	# not found TTA

			&save_DB_info( $dbh, $locus, 'tta__start_coord_tta', $TTAs, ( $data_msg ? qq{data="$data_msg=$msg_id"} : qq{data=NULL} ) );
		}
	}

}


sub search_TTA_codons {
	my( $dbh, $locus, $strand, $sloc, $eloc ) = @_;

	my( $parent_id ) = $dbh->selectrow_array( qq{ SELECT id FROM feats WHERE name="$locus" LIMIT 1 } );
	unless( $parent_id ){
		warn "ATTENTION: '$locus' is not in the GTDB2!	\x1b[31mDISCARDED\x1b[0m.";
		return;
	}

	# Get Nucleotide sequence (CDS)
	# e.g. 'attattttcgct...'
	my $fna_seq = lc( $dbh->selectrow_array( qq{ SELECT data FROM feat_params WHERE name='seq_nt' AND parent_id=$parent_id LIMIT 1 } ) );

	# Search TTA codons
	my @TTAs; # points of TTA_codons
	pos $fna_seq = 0;
	while( $fna_seq=~/tta/g ){
		my $s = $-[0]; # start point of TTA

		next if $s % 3; # Take only ORF TTA coordinate

		# Adjust TTA-point to absolute position
		push @TTAs, ($strand > 0) ? $sloc + $s - 1 : $eloc - $s - 1;
	}

	return \@TTAs;
}


sub save_DB_info {
	my( $dbh, $locus, $tag, $TTA_COORDs, @qq ) = @_;

	my( $parent_id ) = $dbh->selectrow_array( qq{ SELECT id FROM feats WHERE name="$locus" LIMIT 1 } );
	unless( $parent_id ){
		warn "ATTENTION: '$locus' is not in the GTDB2!	\x1b[31mDISCARDED\x1b[0m.";
		return;
	}

	# Collect ALL id of locus
	my @ids;
	my $sth_id = $dbh->prepare( qq{ SELECT id FROM feat_params WHERE name="$tag" AND parent_id=$parent_id ORDER BY id } );
	$sth_id->execute;
	while( my( $id ) = $sth_id->fetchrow_array ) {
		push @ids, $id;
	}
	my $done = $sth_id->finish();

	for my $tta_coord ( sort {$a <=> $b} @$TTA_COORDs ){

		my $qu = join ', ', qq{name="$tag"}, qq{parent_id=$parent_id}, qq{num=$tta_coord}, @qq;

		if( ~~@ids ){
			my $id = shift @ids;
			$dbh->do( qq{ UPDATE feat_params SET $qu WHERE id=$id } ) or warn $dbh->errstr;
		}else{
			$dbh->do( qq{ INSERT feat_params SET $qu } ) or warn $dbh->errstr;
		}
	}

	for( @ids ){	# Сlean extra
		$dbh->do( qq{ DELETE FROM feat_params WHERE id=$_ } ) or warn $dbh->errstr;
	}

}


sub usage {
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = "\x1b[32m" . File::Spec->splitpath($0) . "\x1b[0m";
	return"$msg
$script version $VERSION
    Fills or updates 'feats', and 'feat_params' tables of 'GTDB2' database (default)
    with statistics from the <statistic.tsv> file

USAGE:
    $script [statistic.tsv] [OPTIONS]

EXAMPLE:
    ./$script statistics_NUM_FS_and_TTA_GENES.results.tsv

HERE:
    <statistic.tsv>  --  By default, 'statistics_NUM_FS_and_TTA_GENES.results.tsv',

OPTIONS:
    --cfg_db         --  DB configuration file. By default, 'django_gtdb2/django/mysite/local_settings.json'
    --session        --  session name/key in the DB configuration file. By default, 'default', i.e. 'gtdb2'
";
}

