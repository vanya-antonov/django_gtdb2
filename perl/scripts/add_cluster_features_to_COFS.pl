#!/usr/bin/perl -w

use strict;
use warnings;

use File::Spec;
use Getopt::Long;
use DBI;

# Custom libraries
use MyLibGT::DBapp qw( :all );

my $VERSION = '0.01';

###
# Default Options
my $INFILE  = 'Streptomyces_cluster_features.txt';
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

my $infile = $ARGV[0] || $INFILE || &usage('<cluster_features.file> not specified for processing!');

# Read DataBase configurations
my( $dbh, $gtdb_dir ) = read_configDB( $SESSION, $CFG_DB );

my $save_db = &SaveDB;	# 1 -- Save into 'cofs' table of DB

my @data;
my $name_max;
my $cur_id;

open FTSV, $infile or &usage();
while(<FTSV>)
{
	chomp;
	next if /^$/ || /^#/;	# Drop comments

=comment
# COF_ID	NUM_FEATURES	FEATURE	FS_and_TTA_GENES	ACC_IDs
1006623	5	NA	3597;4385;63652;65439;68745	NZ_CP013129.1;NC_010572.1;NZ_MVFC01000021.1;NZ_CP010407.1;NZ_AJSZ01000509.1
1005426	7	mfs transporter	NZ_AORZ01000011.1_m33953.1256.1;NZ_CP016279.1_p6921417.1268.1;NZ_GG657754.1_m8862580.1259.1;NZ_JNWJ01000118.1_p7603.1253.1;NZ_JNYR01000007.1_m350268.1238.1;NZ_LIQS01000531.1_p4988.1280.1;NZ_MUME01000338.1_m4787.1283.1	NZ_AORZ01000011.1;NZ_CP016279.1;NZ_GG657754.1;NZ_JNWJ01000118.1;NZ_JNYR01000007.1;NZ_LIQS01000531.1;NZ_MUME01000338.1
=cut

	my( $cof_id, $num_features, $feature, $FS_and_TTA_genes, $ACC_IDs ) = split /\t/;

	next if $feature && $feature eq 'NA';

	if( ! defined($cur_id) or $cof_id =~/_+END_+/ or $cur_id ne $cof_id )
	{
		# Search gene in DB
		if( defined($cur_id) && $dbh->selectrow_array( qq{SELECT COUNT(*) AS n FROM cofs WHERE id = $cur_id LIMIT 1} ) )
		{
			next unless $save_db; # Не сохранять

			my $name_clust = $dbh->quote( $name_max );

			my $descr_clust = 'NULL';
			if( @data )
			{
				$descr_clust = join ';', @data;
				for( $descr_clust )
				{
					s/^(.{252}).*/$1\.\.\./ if length($_) > 255;
					$_ = $dbh->quote($_);
				}
			}

			$dbh->do( qq{UPDATE cofs SET c_date=NOW(), name = $name_clust, descr = $descr_clust WHERE id = $cur_id LIMIT 1} );
		}

		@data = ();
		next if $cof_id =~/_+END_+/;

		$name_max = "$feature ($num_features)";
		$cur_id = $cof_id;
	}
	else
	{
		push @data, "$feature ($num_features)";
	}


}
close FTSV;
$dbh->disconnect();

exit;


sub usage
{
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = "\x1b[32m" . File::Spec->splitpath($0) . "\x1b[0m";
	my $text = "$msg
$script version $VERSION
    Fills or updates 'cofs' table of 'GTDB2' database (default)
    with cluster features from <cluster_features.file>

USAGE:
    $script [cluster_features.file]  [OPTIONS]

EXAMPLE:
    $script  Streptomyces_cluster_features.txt
or
    $script -infile Streptomyces_cluster_features.txt

HERE:
    [cluster_features.file] -- file obtained by 'get_cluster_features.pl' script.
                               By default, 'Streptomyces_cluster_features.txt'
OPTIONS:
    --infile                --  Input <cluster_features.file>
    --cfg_db                --  DB configuration file. By default, 'django_gtdb2/django/mysite/local_settings.json'
    --session               --  session name/key in the DB configuration file. By default, 'gtdb2'
    --help

NOTES:
  a configuration DB file 'django_gtdb2/django/mysite/local_settings.json' (or specified by --cfg_db) is required.

";

	die $text;
}

