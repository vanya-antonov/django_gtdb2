#!/usr/bin/perl -w

use strict;
use warnings;
use 5.018;

use File::Spec;
use Getopt::Long;
use DBI;

# Custom libraries
use MyLibGT::DBapp qw( :all );

my $VERSION = '0.01';

###
# Default Options
my $INFILE  = 'Conservative_TTA_codons_in_clusters.STATISTICS.csv';
my $CFG_DB = '/home/alessandro/src/D/django_gtdb2/django/mysite/local_settings.json';
my $SESSION = 'default'; # 'gtdb2_cof'
my $HELP;

###
# Parse input data
GetOptions(
	'infile=s'   => \$INFILE,
	'cfg_db=s'   => \$CFG_DB,
	'session=s'  => \$SESSION,
	'help'       => \$HELP,
) or &usage();

&usage() if $HELP;

my $infile = $ARGV[0] || $INFILE || &usage('<file> not specified for processing!');

# Read DataBase configurations
my( $dbh, $gtdb_dir ) = read_configDB( $SESSION, $CFG_DB );

my $save_db = &SaveDB;	# 1 -- Save into 'cof_params' table of DB

my $i;

# Читаем файл и записываем в DB
open( INFILE, $infile ) or &usage("\x1b[31mERROR\x1b[0m: Can't find $infile. $!");
while(<INFILE>){
	chomp;
	next if /^$/ || /^#/;	# Drop comments

=comment
# cof_id	cluster_size	#_TTA_genes	#_FS-TTA_genes	average_distance
# 1-----	2-----------	3----------	4-------------	5---------------
1004907	131	64	6	172.5
1004909	12	5	0	
=cut

	my( $cof_id, $NUM_ORGS, $NUM_TTA_GENES, $NUM_FS_and_TTA_GENES, $average_distance ) = split /\t/;
	next if $cof_id !~/^\d+/ or ! $NUM_ORGS;

	print ++$i, ") Preparing cluster: $cof_id\n";

	# Check CLASTER in DB
	unless( $dbh->selectrow_array( qq{ SELECT id FROM cof_params WHERE parent_id=$cof_id LIMIT 1 } ) )
	{
		--$i;
		warn"	absent in the database:	SKIPPED\n";
		next;
	}

	$_ ||= 0 for $NUM_TTA_GENES, $NUM_FS_and_TTA_GENES;
;

	my %tags =(
		'num_orgs'              => [ qq{num=$NUM_ORGS} ],
		'tta__num_tta_genes'    => [ qq{num=$NUM_TTA_GENES} ],
		'tta__num_tta_fs_genes' => [ qq{num=$NUM_FS_and_TTA_GENES} ],
	);

	if( $NUM_FS_and_TTA_GENES )
	{
		$average_distance = sprintf "%.1f", $average_distance;
		push @{ $tags{'tta__num_tta_fs_genes'} }, qq{data="avg.Distance=$average_distance"};
	}
	else
	{
		push @{ $tags{'tta__num_tta_fs_genes'} }, qq{data=NULL};
	}

	# Upload into DB (COF_PARAMS table) statistics of FS-,TTA-genes
	for( keys %tags ){
		&save_DB_info( $dbh, $cof_id, $_, $tags{$_} ) if $save_db;
	}

}
close INFILE;
$dbh->disconnect();

die &usage("\x1b[31mERROR\x1b[0m: wrong format or empty input.") unless $i;

exit;


sub save_DB_info {
	my( $dbh, $cof_id, $tag, $qq ) = @_;

	my $qu = join ', ', qq{name="$tag"}, qq{parent_id=$cof_id}, @{ $qq };

	my( $id ) = $dbh->selectrow_array( qq{ SELECT id FROM cof_params WHERE name="$tag" AND parent_id=$cof_id LIMIT 1 } );
	if( $id ){
		$dbh->do( qq{ UPDATE cof_params SET $qu WHERE id = $id } ) or warn $dbh->errstr;
	}else{
		$dbh->do( qq{ INSERT cof_params SET $qu } ) or warn $dbh->errstr;
	}
}

sub usage
{
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = "\x1b[32m" . File::Spec->splitpath($0) . "\x1b[0m";
	my $text = "$msg
$script version $VERSION
    Fills and/or updates 'cof_params' table of 'GTDB2' database (default)
    with statistics from <statistic.tsv> file

USAGE:
    $script [statistic.tsv] [OPTIONS]

EXAMPLE:
    $script  Conservative_TTA_codons_in_clusters.STATISTICS.csv
or
    $script -infile Conservative_TTA_codons_in_clusters.STATISTICS.csv

HERE:
    [statistic.tsv]         -- file obtained by 'search_conserv_TTA_in_cluster_align.pl' script.
                               By default, 'Conservative_TTA_codons_in_clusters.STATISTICS.csv'
OPTIONS:
    --infile                --  Input <statistic.tsv>
    --cfg_db                --  DB configuration file. By default, 'django_gtdb2/django/mysite/local_settings.json'
    --session               --  session name/key in the DB configuration file. By default, 'gtdb2'
    --help

NOTES:
  a configuration DB file 'django_gtdb2/django/mysite/local_settings.json' (or specified by --cfg_db) is required.

";

	die $text;
}

