#!/usr/bin/perl -w

use strict;
use warnings;

use File::Spec;
use Getopt::Long;
use DBI;
use JSON;

my $VERSION = '1.02';

###
# Default Options
my $INFILE = 'statistics_NUM_FS_and_TTA_GENES.results.tsv';
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
my( $dbh, $gtdb_dir ) = &read_configDB( $SESSION, $CFG_DB );

my $save_db = &SaveDB;

# FS_ID     	inner ID for FS-gene
# TTA_FS_gene	FS-gene ID, в формате SEQ_ID:START-END
# SLOC      	start location
# ELOC      	end location
# STRAND    	strand
# SEQ_ID    	ID последовательности
# ORG_ID    	inner ID организма
# ORG_NAME  	название организма
# FSHIFT_COORD	абсолютная координата СРС по нт последовательности
# SOURCE    	source of sequence for TTA-codon search (gtdb2 | GBK)
# TTA_COORD 	абсолютная(ые) координата кодона ТТА по нт последовательности
# DIST_FS_TTA	среднее расстояние между СРС и ТТА кодоном (в нт последовательности)
# LINK_TTA1 	ссылка_1 на NCBI на кодон ТТА
# LINK_TTA2 	ссылка_2 на NCBI на кодон ТТА (если есть 2 TTA-codons)
#
open INFILE, $INFILE or die "Can't found $INFILE: $!\n";
while(<INFILE>){
	s/^\s+|\s+$//g;
	next if /^$/ || /^#/;
	next unless /^\d/;
=comment
# FS_ID	TTA_FS_gene	SLOC	ELOC	STRAND	SEQ_ID	ORG_ID	ORG_NAME	FSHIFT_COORD	SOURCE	TTA_COORD	DIST_FS_TTA	LINK_TTA1	LINK_TTA2
# 1----	2----------	3---	4---	5-----	6-----	7-----	8-------	9-----------	10----	11-------	12---------	13-------	14-------
74297	SAVERM_RS00545=NC_003155.5:25699-30394	29513	30391	-1	NC_003155.5	167	Streptomyces avermitilis MA-4680 = NBRC 14893	30205	GBK	25702	4502	=HYPERLINK("https://www.ncbi.nlm.nih.gov/nuccore/NC_003155.5?report=graph&from=25699&to=25705", "TTA_codon")
=cut

	my( $fs_id, $TTA_FS_gene, $SLOC, $ELOC, $strand, $SEQ_ID, $ORG_ID, $ORG_NAME, $FSHIFT_COORD, $SOURCE, $TTA_COORD, $DIST_FS_TTA ) = split /\t/;

	my( $locus_tag, $location ) = split /=/, $TTA_FS_gene;	# SAVERM_RS00545=NC_003155.5:25699-30394

	next unless $save_db;

	# Создание для гена/локуса записей в таблицах seqs, feats, feat_params
	`python3 $DJANGO/manage.py  run_method_on_object  get_or_create_feat_from_locus_tag  Org $ORG_ID $locus_tag CDS`;

	$TTA_COORD =~s/;\d+$//; # Оставляем только 1й TTA (before)
	&save_DB_info( $dbh, $locus_tag, 'tta__start_coord_tta', qq{num=$TTA_COORD}, qq{data="FS_ID=$fs_id;SOURCE=$SOURCE"} );

}
close INFILE;

exit;


sub save_DB_info {
	my( $dbh, $locus_tag, $tag, @qq ) = @_;

	my( $parent_id ) = $dbh->selectrow_array( qq{ SELECT id FROM feats WHERE name="$locus_tag" LIMIT 1 } );

	my $qu = join ', ', qq{name="$tag"}, qq{parent_id=$parent_id}, @qq;

	my( $id ) = $dbh->selectrow_array( qq{ SELECT id FROM feat_params WHERE name="$tag" AND parent_id=$parent_id LIMIT 1 } );
	if( $id ){
		$dbh->do( qq{ UPDATE feat_params SET $qu WHERE id=$id } ) or warn $dbh->errstr;
	}else{
		$dbh->do( qq{ INSERT feat_params SET $qu } ) or warn $dbh->errstr;
	}

}


sub SaveDB {
	my $save_db = 1;	# 1=save into DB

	while( 1 ){
		last unless $save_db;

		print"Do You want save/update info DB (\x1b[31;1m y\x1b[0m or \x1b[32;1m n\x1b[0m )? ";
		$_ = <STDIN>;
		last if /^y(?:es)?/i;
		$save_db = 0 if /^no?/i;
	}

	$save_db;
}


sub read_configDB {
	my( $session, $cfgs ) = @_;

	# Read DataBase configuration
	my $json_set = do {
		open( my $json_fh, $cfgs )
			or die "\x1b[31mERROR\x1b[0m: Can't open $cfgs file: $!";

		local $/;
		<$json_fh>;
	};

	my $ref = decode_json( $json_set );
	die "\x1b[31mERROR\x1b[0m: No DB settings exist: $!"
		if !exists( $ref->{'DATABASES'} ) or !exists( $ref->{'DATABASES'}{ $session } );

	my $db_name  = $ref->{'DATABASES'}{ $session }{'NAME'}     || 'gtdb2'; # 'gtdb2_cof'
	my $user     = $ref->{'DATABASES'}{ $session }{'USER'};
	my $password = $ref->{'DATABASES'}{ $session }{'PASSWORD'};
	my $host     = $ref->{'DATABASES'}{ $session }{'HOST'}     || 'localhost';
	my $engine   = $ref->{'DATABASES'}{ $session }{'DBI'}      || 'DBI:mysql:database';
	my $port     = $ref->{'DATABASES'}{ $session }{'PORT'}     || 3306;

	my $dbh = DBI->connect("$engine=$db_name;host=$host;port=$port", $user, $password,
						{ RaiseError => 0, PrintError => 1, AutoCommit => 1} );

	my $gtdb_dir = $ref->{'DATABASES'}{ $session }{'GTDB_DIR'}; # "/home/gtdb/data/gtdb2"
	$gtdb_dir =~s/\/+$//;

	return( $dbh, $gtdb_dir );
}


sub usage {
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = "\x1b[32m" . File::Spec->splitpath($0) . "\x1b[0m";
	return"$msg
$script version $VERSION
    Fills or updates 'feats', 'feat_params' tables of 'GTDB2' database (default) with statistics
    from the <statistic.tsv> file

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

