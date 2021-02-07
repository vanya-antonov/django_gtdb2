#!/usr/bin/perl --

use strict;
use warnings;

# $Id$

###
# Ivan Antonov (antonov1986@gmail.com)
#
# Наполняет / Изменяет таблицы:
#	'cof_gtfs' (необходима, из-за хранения промежуточных кластеров)
#	'cofs'
# Использует таблицы:
#	'cof_hits'
#	'fshifts'

$|++; # Turn off buffering

# http://www.shlomifish.org/lecture/LM-Solve/slides/exotic-bugs/recursion_limit.html
$DB::deep = 10000;
no warnings "recursion";

use Data::Dumper;
use File::Spec;
use Cwd 'abs_path';

use MyLib::classes::GeneTackDB;
use MyLib::BaseUtil qw(ah2a);

die &usage() if @ARGV != 2;

###
# CONSTANTS
my $LABEL_ID = undef;
my $USER_ID  = undef;


###
my $START_TIME = time;

&run(
	evalue_thr => $ARGV[0],
	subset     => $ARGV[1],
	label      => $LABEL_ID,
	user_id    => $USER_ID,
);

warn "\nElapsed time: ".(time - $START_TIME)." sec\n";
exit;

###
# SUBROUTINES
sub run
{
	my %opts = @_;

	my $gtdb = MyLib::classes::GeneTackDB->new( MyLib::BaseUtil->new() );

	my $all_fs_ids = [];
	my %fs_only    = ();
	if( $opts{subset} eq '__ALL__' )
	{
		$all_fs_ids = ah2a('QID', $gtdb->{bu}->exec_SQL_ar('SELECT DISTINCT q_fs_id AS qid FROM cof_hits') );

#		warn "Removing COFs....";
#		$gtdb->delete_cof($_->{COF_ID}) foreach @{$gtdb->{bu}->exec_SQL_ar('select cof_id from cofs')};
	}
	elsif( -e $opts{subset} ) # for file FS_ONLY.txt
	{
		@$all_fs_ids = grep { length } split /[\n\r]+/, `cat $opts{subset}`;
	}
	else		# for COF_ID
	{
		$all_fs_ids = $gtdb->get_all_fs_for_cof( $opts{subset} );
#		$all_fs_ids = ah2a('ID', $gtdb->{bu}->exec_SQL_ar('SELECT DISTINCT id FROM fshifts WHERE cof_id=?', $opts{subset} ) );

		$gtdb->delete_cof( $opts{subset} );

		%fs_only = map { $_ => 1 } @$all_fs_ids;
	}

	&find_all_cofs( $gtdb, $opts{evalue_thr}, $all_fs_ids,
		label   => $opts{label},
		user_id => $opts{user_id},
		fs_only => keys(%fs_only) ? \%fs_only : undef,
	);

	$gtdb->{bu}->commit;
}


sub find_all_cofs
{
	my( $gtdb, $evalue_thr, $all_fs_ids, %opts ) = @_;
	my $bu = $gtdb->{bu};

	while(my $q_id = shift @$all_fs_ids )
	{
		print "[".@$all_fs_ids."] Processing query $q_id...\n";

		# If this FS is already clustered
		warn "\tQuery $q_id is already in COF\n" and next if $gtdb->get_cofs_for_fs($q_id)->[0];

		my $q_TYPE = $gtdb->{bu}->exec_SQL_ar('SELECT len FROM fshifts WHERE id=?', $q_id )->[0]{LEN};

		warn "FS $q_id does not have type!!!" and next unless $q_TYPE;

		my $cluster = [];
		my $visited_ids = {};
		&add_query_to_cluster($gtdb, $q_id, $cluster, $evalue_thr, $q_TYPE, $visited_ids, %opts);

		if( @$cluster > 1) # cluster consists more 1 FS
		{
			my $cof = $gtdb->create_new_cof( $cluster, user_id => $opts{user_id} );

			print "\tCOF $cof created with ".@$cluster." elements\n";
		}
	}
}


sub add_query_to_cluster
{
	my( $gtdb, $q_id, $cluster, $evalue, $type, $visited_ids, %opts ) = @_;

	$visited_ids = {} unless $visited_ids;

	return if exists $visited_ids->{$q_id};

	$visited_ids->{$q_id} = undef;

	return if defined( $opts{fs_only} ) and ! exists( $opts{fs_only}{$q_id} );

	# if it is already in COFS
	if( my $cof_id = $gtdb->get_cofs_for_fs($q_id)->[0] )
	{
		# "$q_id is already in a COF -- mearge this COF into a new one
		my $cof_fs = $gtdb->get_all_fs_for_cof( $cof_id );

		warn "\tMearging with COF $cof_id (".@$cof_fs." elements)\n";

		push @$cluster, @$cof_fs;
		$visited_ids->{$_} = undef foreach @$cof_fs;

		$gtdb->delete_cof( $cof_id );

		return;
	}

	push @$cluster, $q_id;

	# Try to add hits of its hits  --  recursion!!!
	my $hits = $gtdb->{bu}->exec_SQL_ar(q[
		SELECT DISTINCT ch.h_fs_id AS FS_ID, ch.evalue, fs.len AS TYPE
		FROM cof_hits AS ch, fshifts AS fs
		WHERE ch.q_fs_id=?
		AND fs.id = ch.h_fs_id
		AND fs.len=?
	], $q_id, $type);

	for( @$hits ){
		next if $_->{EVALUE} > $evalue;

		&add_query_to_cluster( $gtdb, $_->{FS_ID}, $cluster, $evalue, $type, $visited_ids, %opts );
	}
}


sub usage
{
	my($msg) = @_;
	$msg = $msg ? "$msg\n" : '';
	my $script = File::Spec->splitpath($0);
	return"$msg
DESCRIPTION:
    The script finds new COFs from all the fs (__ALL__) or it can split existing big COFs (if COF_ID specified)

    FS_ONLY.txt  --  recluster only these frameshfits (NO HEADER!)

USAGE:
    $script  <EVALUE_THR>   <FS_ONLY.txt | COF_ID | __ALL__>

EXAMPLE:
    $script  1e-50  __ALL__

";
}

