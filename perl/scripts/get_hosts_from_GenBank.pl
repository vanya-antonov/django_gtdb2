#!/usr/bin/perl -w

use strict;
use warnings;
no strict 'refs';

use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use DBI;
use JSON;

my $VERSION = '1.04';

###
# Default Options
my $OUTPUT;
my $AUTO;
my $HEADER;
my $ECHO;
my $CFG_DB = '/home/alessandro/src/D/django_gtdb2/django/mysite/local_settings.json';
my $SESSION = 'default'; # 'gtdb2_cof'

###
# Parse input data
GetOptions(
	'cfg_db=s'   => \$CFG_DB,
	'session=s'  => \$SESSION,
	'output=s'   => \$OUTPUT,
	'auto'       => \$AUTO,
	'header'     => \$HEADER,
	'echo'       => \$ECHO,
) or die &usage();

my $infile = $ARGV[0] || die &usage(); # 'NZ_CP023977.1.gbk';

# Check GenBank format
open GBF, $infile or die &usage();
$_ = <GBF>;
close GBF;

die "\x1b[31mERROR\x1b[0m: Invalid GenBank file format" unless /^LOCUS\s+/;

# Read DataBase configurations
my( $dbh, $gtdb_dir ) = &read_configDB( $SESSION, $CFG_DB );

if( defined $AUTO ){
	(my $fl = $infile) =~s/\.[^\.]+$//;

	$OUTPUT ||= "$fl.tsv";
}

$OUTPUT = undef if $OUTPUT && $OUTPUT =~/^stdout$/i;

###
my $START_TIME = time;

my $seqio_obj = Bio::SeqIO->new(
	-file   => $infile,
	-format => 'genbank'
	);

my %all;

while( my $seq_obj = $seqio_obj->next_seq ){

	my $acc_id = $seq_obj->accession;
	my $ver = $seq_obj->seq_version;
	$acc_id .= ".$ver" if defined $ver;

	my $species_string = $seq_obj->species->node_name; # ORGANISM
	for( $species_string ){
		s/^\s+|\s+$//g;
		s/\s+/ /g;
	}

	my %dt;
	my( $organism, $taxon, @hosts );

	for my $feat_obj ($seq_obj->get_SeqFeatures) {
		my $primary_tag = $feat_obj->primary_tag;

		if( $primary_tag eq 'source'){

			# from a line like: /organism="Streptomyces griseus subsp. griseus NBRC 13350"
			if( $feat_obj->has_tag('organism') ){
				for( $feat_obj->get_tag_values('organism') ){
					s/^\s+|\s+$//g;
					s/\s+/ /g;
					$organism = $_;
					last;
				}
			}

			if( $feat_obj->has_tag('db_xref') ){
				for( $feat_obj->get_tag_values('db_xref') ){
					if( /^taxon:(\d+)/ ){ # from a line like: /db_xref="taxon:455632"
						$taxon = $1;
						$taxon =~s/\D+//g;
						last;
					}
				}
			}

			# from a line like: /host="Streptomyces venezuelae ATCC 10712"
			if( $feat_obj->has_tag('host') ){
				for( $feat_obj->get_tag_values('host') ){
					s/^\s+|\s+$//g;
					s/\s+/ /g;
					push @hosts, $_;
				}
			}

			# from a line like: /lab_host="Streptomyces venezuelae ATCC 10712"
			if( $feat_obj->has_tag('lab_host') ){
				for( $feat_obj->get_tag_values('lab_host') ){
					s/^\s+|\s+$//g;
					s/\s+/ /g;
					push @hosts, $_;
				}
			}

			next;
		}

	}

# 1
	$taxon ||= $species_string;
# 2
	$all{ $taxon }{'ORG_NAME'} = $organism;
# 3, 4
	for my $host ( @hosts ){	# 'Streptomyces venezuelae ATCC 10712'
		my $h = $host;
		for( $h ){
			s/\s*Streptomyces\*//;
			s/\s+\d+//g;

			s/^\s+|\s+$//g;
			s/\s+/ /g;
		}

		my @hh = split /\s+/, $h;
		while( ~~@hh ){
			my $msk = join ' ', @hh;

			my( $host_id, $name ) = $dbh->selectrow_array( qq{ SELECT id, name FROM orgs WHERE name RLIKE "$msk" LIMIT 1 } );
			if( $host_id ){
				$all{ $taxon }{'IDs'}{ $host_id }{'gbk'} = $host;
				$all{ $taxon }{'IDs'}{ $host_id }{'db'} = $name;
				last;
			}
			pop @hh;
		}

	}

	if( ! exists( $all{ $taxon }{'IDs'} ) and ~~@hosts ){
		$all{ $taxon }{'IDs'}{-1}{'gbk'} = join ';', @hosts;
	}
}

print "\n# Elapsed time: ".(time - $START_TIME)." sec\n" if $ECHO;

exit unless scalar( keys %all );

my $head_out = join "\t", qw{TAXON ORG_NAME HOST_IDs HOST_NAMES_db HOST_NAMES_gbk};

# Save collection
if( $OUTPUT ){
	open OFILE, ">$OUTPUT";
	print OFILE "$head_out\n" if $HEADER;
}else{
	print "$head_out\n" if $HEADER;
}

for my $taxon ( sort keys %all ){
	my( @HOST_IDs, @HOST_NAMES_db, @HOST_NAMES_gbk );
	for( sort{ $a <=> $b } keys %{ $all{ $taxon }{'IDs'} } ){
		push @HOST_IDs, $_;
		push @HOST_NAMES_db, $all{ $taxon }{'IDs'}{$_}{'db'};
		push @HOST_NAMES_gbk, $all{ $taxon }{'IDs'}{$_}{'gbk'};
	}

	my $out = join("\t", $taxon, $all{ $taxon }{'ORG_NAME'}, join(';', @HOST_IDs), join(';', @HOST_NAMES_db), join(';', @HOST_NAMES_gbk) );

	if( $OUTPUT ){
		print OFILE "$out\n";
	}else{
		print "$out\n";
	}
}

exit;


# Read DataBase configuration
sub read_configDB {
	my( $session, $cfgs ) = @_;

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
	my $engine   = $ref->{'DATABASES'}{ $session }{'DBI'}    || 'DBI:mysql:database';
	my $port     = $ref->{'DATABASES'}{ $session }{'PORT'}   || 3306;

	my $dbh = DBI->connect("$engine=$db_name;host=$host;port=$port", $user, $password,
						{ RaiseError => 0, PrintError => 1, AutoCommit => 1} );

	my $gtdb_dir = $ref->{'DATABASES'}{ $session }{'GTDB_DIR'}; # "/home/gtdb/data/gtdb2/"
	$gtdb_dir =~s/\/+$//;

	return( $dbh, $gtdb_dir );
}


sub usage
{
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = "\x1b[32m" . File::Spec->splitpath($0) . "\x1b[0m";
	return"$msg
$script version $VERSION

USAGE:
    $script <file.gbk> [OPTIONS]

EXAMPLE:
    $script NC_010572.1.gbk -auto -output=stdout

    $script NC_010572.1.gbk -auto
or same
    $script NC_010572.1.gbk -output NC_010572.1.tsv

HERE:
    <file.gbk>   -- input GenBank file only

OPTIONS:
    --output  <file.tsv|stdout>  --  output table. Default STDOUT output
    --header                     --  output header of table
    --echo                       --  output echo messages
    --auto                       --  autocomplete options: --output
    --cfg_db                     --  DB configuration file. By default, 'django_gtdb2/django/mysite/local_settings.json'
    --session                    --  session name/key in the DB configuration file. By default, 'gtdb2'

OUTPUT TABLE FORMAT:
  1.TAXON                        --  NCBI taxon ID of organism from GenBank
  2.ORG_NAME                     --  Organism (virus) name from GenBank
  3.HOST_IDs                     --  Host identifier(;s) for HOST_NAMES_db from GTDB
  4.HOST_NAMEs_db                --  Host name(;s) existing in GTDB
  5.HOST_NAMES_gbk               --  Host name(;s) from GenBank
";
}

