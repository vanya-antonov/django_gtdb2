#!/usr/bin/perl --

use strict;
use warnings;

# $Id$

###
# Ivan Antonov (antonov1986@gmail.com)
#

$| = 1; # Turn off buffering

use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Cwd 'abs_path';

use Bio::Seq;
use Bio::SearchIO;

use MyLib::classes::GeneTackDB;
use MyLib::BaseUtil qw(ah2a min);

###
# CONSTANTS
my $MIN_FS_DIST = 10;   # FS_DIST -- minimal distance from FS to BLAST alignment border

###
# Parse input data
GetOptions() || die usage();

die usage() if @ARGV!=1;

###
my $START_TIME = time;

run(
	blast_fn   => $ARGV[0],
);

warn "\nElapsed time: ".(time-$START_TIME)." sec\n";
###

###
# SUBROUTINES
sub run
{
	my %opts = @_;

	my $bu   = MyLib::BaseUtil->new();
	my $gtdb = MyLib::classes::GeneTackDB->new( $bu );

#	$bu->exec_SQL_nr('ALTER TABLE cof_hits DROP FOREIGN KEY cofhits_qfsid_fk;');
#	$bu->exec_SQL_nr('ALTER TABLE cof_hits DROP FOREIGN KEY cofhits_hfsid_fk;');
#	$bu->exec_SQL_nr('ALTER TABLE cof_hits DROP INDEX cofhits_qfsid_hfsid_i;');
#	$bu->exec_SQL_nr('delete from cof_hits');

	hits2db($gtdb, $opts{blast_fn});
	$bu->commit;

	# Check for cofhits_qfsid_fk index existence (MyLib::classes::GeneTackDB)
	my $chk = $bu->exec_SQL_ar( q[SELECT COUNT(CONSTRAINT_NAME) AS N FROM information_schema.KEY_COLUMN_USAGE
			WHERE CONSTRAINT_SCHEMA='gtdb2_cof' AND TABLE_NAME='cof_hits' AND COLUMN_NAME=?], 'q_fs_id')->[0]{N};

	$bu->exec_SQL_nr(q[ALTER TABLE cof_hits ADD CONSTRAINT FOREIGN KEY cofhits_qfsid_fk (q_fs_id) REFERENCES fshifts(id);])
if !$chk or $chk < 2;

	# Check for Index existence
	$chk = $bu->exec_SQL_ar( q[SELECT COUNT(CONSTRAINT_NAME) AS N FROM information_schema.KEY_COLUMN_USAGE
			WHERE CONSTRAINT_SCHEMA='gtdb2_cof' AND TABLE_NAME='cof_hits' AND COLUMN_NAME=?], 'h_fs_id')->[0]{N};

	$bu->exec_SQL_nr(q[ALTER TABLE cof_hits ADD CONSTRAINT FOREIGN KEY cofhits_hfsid_fk (h_fs_id) REFERENCES fshifts(id);])
if !$chk or $chk < 2;

	# Check for Index existence
	unless( $bu->exec_SQL_ar( q[SELECT COUNT(1) AS N FROM information_schema.STATISTICS
WHERE table_schema=DATABASE() AND table_name='cof_hits' AND index_name=?], 'cofhits_qfsid_hfsid_i')->[0]{N} ){

	    $bu->exec_SQL_nr('CREATE UNIQUE INDEX cofhits_qfsid_hfsid_i ON cof_hits(q_fs_id, h_fs_id);');
	}

	$bu->commit;
}

sub hits2db
{
	my($gtdb, $blast_fn) = @_;

	my $in = new Bio::SearchIO(-format => 'blast', -file => $blast_fn);
	while( my $res = $in->next_result )
	{
		my( $q_id )  = $res->query_accession =~ /(\d+)/;

#		my $q_info = $gtdb->get_ALL_info_about_fs($q_id); # MyLib::classes::GeneTackDB

		my $q_prot_fs_coord = $gtdb->{bu}->exec_SQL_ar(
			q[SELECT num FROM fshift_params WHERE name='seq_prot_n' AND parent_id=?], $q_id)->[0]{NUM};

		next unless $q_prot_fs_coord;
#		die "Unknown FS_ID '$q_id'" if !$q_prot_fs_coord;

		warn "Processing hits for '$q_id'...\n";
		foreach my $tp ( @{ MyLib::classes::GeneTackDB::get_all_TP_hit_hsp($res, $q_prot_fs_coord, $MIN_FS_DIST) } )
		{
			my $hsp   = $tp->{hsp};
			my $hit   = $tp->{hit};
			my($h_id) = $hit->accession =~ /(\d+)/;
			next if $h_id == $q_id;

			# my $h_info = $gtdb->get_ALL_info_about_fs($h_id);

			my $h_prot_fs_coord = $gtdb->{bu}->exec_SQL_ar( q[
			SELECT num FROM fshift_params WHERE name='seq_prot_n' AND parent_id=?], $h_id )->[0]{NUM};


			next unless $h_prot_fs_coord;
#			die "Unknown FS_ID '$h_id'" if ! $h_prot_fs_coord;

			# check if hit_FS is far enough from ali border
			my $h_fs_dist = min( $h_prot_fs_coord - $hsp->start('hit'), $tp->{hsp}->end('hit') - $h_prot_fs_coord );
			next if $h_fs_dist < $MIN_FS_DIST;

			my $q_fs_ali_coord = prot_coord_to_ali_coord( $q_prot_fs_coord - $hsp->start('query') + 1, $hsp->query_string );
			my $h_fs_ali_coord = prot_coord_to_ali_coord( $h_prot_fs_coord - $hsp->start('hit') + 1,   $hsp->hit_string );

			# if hit FS is not int he ali OR hit_FS and q_FS are too far apart
			next if !$h_fs_ali_coord || abs( $q_fs_ali_coord - $h_fs_ali_coord ) > 50;

			my $q_start = $tp->{hsp}->start('query');
			my $q_end   = $tp->{hsp}->end('query');
			my $h_start = $tp->{hsp}->start('hit');
			my $h_end   = $tp->{hsp}->end('hit');
			my $ali_len = $tp->{hsp}->length('total');
			my $evalue  = $tp->{hsp}->evalue;
			my $identity= $tp->{hsp}->percent_identity;
			my $score   = $tp->{hsp}->score;

			$gtdb->{bu}->exec_SQL_nr( q[
				INSERT INTO cof_hits
				(q_fs_id, q_start, q_end, q_fs_ali_coord, h_fs_id, h_start, h_end, h_fs_ali_coord, ali_len, score, evalue, identity)
				VALUES(?,?,?,?,?,?,?,?,?,?,?,?)
			], $q_id, $q_start, $q_end, $q_fs_ali_coord, $h_id, $h_start, $h_end, $h_fs_ali_coord, $ali_len, $score, $evalue, $identity
			);
		}
	}
}


sub prot_coord_to_ali_coord
{
	my($prot_coord, $ali_str) = @_;

	my $ali_len      = length($ali_str);
	my $ali_coord    = $prot_coord;
	my $coord_prefix = substr($ali_str, 0, $ali_coord);
	while( $coord_prefix =~ s/-//g )
	{
		my $extra_len = $prot_coord-length($coord_prefix);

		# FS coord is not covered at all -- happens in the case of hits
		return undef if $ali_coord+$extra_len > $ali_len;

		$coord_prefix = $coord_prefix.substr($ali_str, $ali_coord, $extra_len);
		$ali_coord += $extra_len;
	}

	return $ali_coord;
}

sub usage
{
	my($msg) = @_;
	$msg = $msg ? $msg."\n" : '';
	my $script = File::Spec->splitpath($0);
	return"$msg
DESCRIPTION:

USAGE:
    $script   <HITS.blast>

OPTIONS:

";
}

