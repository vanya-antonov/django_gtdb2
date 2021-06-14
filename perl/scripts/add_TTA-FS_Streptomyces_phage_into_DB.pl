#!/usr/bin/perl -w

use strict;
use warnings;

use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use DBI;
use JSON;

my $VERSION = '1.03';

###
# Default Options
my $CFG_DB = '/home/alessandro/src/D/django_gtdb2/django/mysite/local_settings.json';
my $SESSION = 'default'; # 'gtdb2_cof'

###
# Parse input data
GetOptions(
	'cfg_db=s'       => \$CFG_DB,
	'session=s'      => \$SESSION,
) or die &usage();

my $infile = $ARGV[0] || 'statistic_Streptomyces_phages_TTA-FS_genes.tsv';

# Read DataBase configurations
my( $dbh, $gtdb_dir ) = &read_configDB( $SESSION, $CFG_DB );

my $save_db = &SaveDB;

my $i = 0;

# Читаем файл и записываем в DB
open( INFILE, $infile ) or die &usage("\x1b[31mERROR\x1b[0m: Can't find $infile. $!");
while(<INFILE>){
	s/^\s+|\s+$//g;
	next if /^$/ || /^#/;

	if(/^PHAGE_ID/){
		next if $_ eq "PHAGE_ID	PHAGE_NAME	PHAGE_GENOME_LEN	HOST_ID(gtdb2)	HOST_NAME	PHAGE_NUM_GENES	PHAGE_NUM_TTA_GENES	PHAGE_NUM_FS_GENES	PHAGE_NUM_FS_and_TTA_GENES	FS_and_TTA_LIST(;)";
		&usage("\x1b[31mERROR\x1b[0m: wrong format of input. $!");
	}

=comment
PHAGE_ID	PHAGE_NAME	PHAGE_GENOME_LEN	HOST_ID(gtdb2)	HOST_NAME	PHAGE_NUM_GENES	PHAGE_NUM_TTA_GENES	PHAGE_NUM_FS_GENES	PHAGE_NUM_FS_and_TTA_GENES	FS_and_TTA_LIST(;)
5895	Streptomyces phage VWB	49220	3718	Streptomyces venezuelae	61	0	5	0	
5921	Streptomyces phage Joe	48941	972	Streptomyces coelicolor A3(2)	81	2	4	2	39050-39574:-:39435:-1;39050-39574:-:39342:+1
=cut

	my( $PHAGE_ID, $PHAGE_NAME, $PHAGE_GENOME_LEN, $HOST_ID_gtdb2, $HOST_NAME, $PHAGE_NUM_GENES,
$PHAGE_NUM_TTA_GENES, $PHAGE_NUM_FS_GENES, $PHAGE_NUM_FS_and_TTA_GENES, $FS_and_TTA_LIST ) = split /\t/;

	next if $PHAGE_ID !~/^\d/;

	print ++$i, ") Preparing phage: $PHAGE_NAME\n";

	unless( $PHAGE_NUM_GENES ){
		warn "	PHAGE_NUM_GENES = 0:	SKIPPED\n";
		next;
	}

	# Get phage in DB
	unless( $dbh->selectrow_array( qq{ SELECT id FROM org_params WHERE name='taxonomy' AND value='Viruses' AND parent_id = $PHAGE_ID LIMIT 1 } ) ){
		warn"	absent in the database:	SKIPPED\n";
		next;
	}

	for( $PHAGE_NUM_TTA_GENES, $PHAGE_NUM_FS_GENES, $PHAGE_NUM_FS_and_TTA_GENES ){
		$_ ||= 0;
	}

	my $qu_FS_and_TTA_LIST = 'data=' . ($FS_and_TTA_LIST ? qq{"$FS_and_TTA_LIST"} : 'NULL');

	next unless $save_db;

	my %tags =(
		'tta__virus_host_name'  => [ qq{value="$HOST_NAME"}, qq{num=$HOST_ID_gtdb2} ],
		'num_annotated_genes'   => [ qq{num=$PHAGE_NUM_GENES} ],
		'tta__num_tta_genes'    => [ qq{num=$PHAGE_NUM_TTA_GENES} ],
		'num_fshifts'           => [ qq{num=$PHAGE_NUM_FS_GENES} ],
		'tta__num_tta_fs_genes' => [ qq{num=$PHAGE_NUM_FS_and_TTA_GENES}, $qu_FS_and_TTA_LIST ],
	);
	# Upload into DB (ORG_PARAMS table) statistics TTA-genes
	for( keys %tags ){
		&save_DB_info( $dbh, $PHAGE_ID, $_, $tags{$_} );
	}

}
close INFILE;
$dbh->disconnect();

exit;


sub save_DB_info {
	my( $dbh, $PHAGE_ID, $tag, $qq ) = @_;

	my $qu = join ', ', qq{name="$tag"}, qq{parent_id=$PHAGE_ID}, @{ $qq };

	my( $id ) = $dbh->selectrow_array( qq{ SELECT id FROM org_params WHERE name="$tag" AND parent_id = $PHAGE_ID LIMIT 1 } );
	if( $id ){
		$dbh->do( qq{ UPDATE org_params SET $qu WHERE id = $id } ) or warn $dbh->errstr;
	}else{
		$dbh->do( qq{ INSERT org_params SET $qu } ) or warn $dbh->errstr;
	}
}


sub read_configDB {
	my( $session, $cfgs ) = @_;

	# Read DataBase configuration
	my $json_set = do {
		open( my $json_fh, $cfgs )
			or &usage("\x1b[31mERROR\x1b[0m: Can't open $cfgs file. $!");

		local $/;
		<$json_fh>;
	};

	my $ref = decode_json( $json_set );
	&usage("\x1b[31mERROR\x1b[0m: No DB settings exist. $!")
		if !exists( $ref->{'DATABASES'} ) or !exists( $ref->{'DATABASES'}{ $session } );

	my $db_name  = $ref->{'DATABASES'}{ $session }{'NAME'}   || 'gtdb2'; # 'gtdb2_cof'
	my $user     = $ref->{'DATABASES'}{ $session }{'USER'};
	my $password = $ref->{'DATABASES'}{ $session }{'PASSWORD'};
	my $host     = $ref->{'DATABASES'}{ $session }{'HOST'}   || 'localhost';
	my $engine   = $ref->{'DATABASES'}{ $session }{'DBI'} || 'DBI:mysql:database';
	my $port     = $ref->{'DATABASES'}{ $session }{'PORT'}   || 3306;

	my $dbh = DBI->connect("$engine=$db_name;host=$host;port=$port", $user, $password,
						{ RaiseError => 0, PrintError => 1, AutoCommit => 1} );

	my $gtdb_dir = $ref->{'DATABASES'}{ $session }{'GTDB_DIR'}; # "/home/gtdb/data/gtdb2/"

	return( $dbh, $gtdb_dir );
}


sub SaveDB {
	my $save_db = 1;	# 1=save into DB

	while( 1 ){
		last unless $save_db;

		print"Do You want save/update info DB (y or n)? ";
		$_ = <STDIN>;
		last if /^y(?:es)?/i;
		$save_db = 0 if /^no?/i;
	}

	$save_db;
}


sub usage {
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = "\x1b[32m" . File::Spec->splitpath($0) . "\x1b[0m";
	return"$msg
$script version $VERSION

USAGE:
    $script [statistic.tsv] [OPTIONS]

EXAMPLE:
    $script statistic_Streptomyces_phages_TTA-FS_genes.tsv

HERE:
    <statistic.tsv>  --  By default, 'statistic_Streptomyces_phages_TTA-FS_genes.tsv'

OPTIONS:
    --cfg_db         --  DB configuration file. By default, 'django_gtdb2/django/mysite/local_settings.json'
    --session        --  session name/key in the DB configuration file. By default, 'default', i.e. 'gtdb2'
";
}
