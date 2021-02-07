package MyLib::classes::GeneTackDB;

use strict;
use warnings;

# $Id$

###
# Ivan Antonov (antonov1986@gmail.com)
#
###
#
# $self
#   |
#   |--->{tbl}  --  MyLib::classes::Table object
#   |
#   |--->{all_fs}   --   reference to hash of hashes
#   |       |
#   |       |--->{ $fs_coord }
#   |                  |
#   |                  |--->{FS_coord}
#   |                  |
#   |                 ...
#   |             Other keys according to headers in 
#   |             the .fsmarkgm file
#   |                 ...
#   |
#   |--->{CACHE}
#   |       |
#

use Data::Dumper;
use Carp qw(confess);
use Bio::SeqIO;
use MyLib::BaseUtil qw(min ah2a);


##############
# CONSTRUCTOR
sub new
{
	my $class = shift;
	my($bu) = @_;

	my $self = bless {
		bu  => $bu,
		dir => '/home/alessandro/gtdb',
	}, $class;

	return $self;
}

#################
# PUBLIC METHODS

sub create_new_cof
{
	my $self = shift;
	my( $fs_arr, %opts ) = @_;

	my $cof_id  = $self->{bu}->get_random_DB_ID('cofs', 'cof_id');
	my $cof_dir = "$self->{dir}/cofs/$cof_id";
	`mkdir -p $cof_dir`;

	$self->{bu}->exec_SQL_nr(q[
		INSERT INTO cofs(cof_id, dir, user_id, c_date) VALUES(?,?,?,NOW())
	], $cof_id, $cof_dir, $opts{user_id} );
	foreach my $fs_id ( @$fs_arr )
	{
		$self->{bu}->exec_SQL_nr('INSERT INTO cof_gtfs(cof_id, fs_id, is_core) VALUES(?,?,?)', $cof_id, $fs_id, 1);
	}

	# calculate some values
	$self->get_ALL_info_about_cof( $cof_id );

	return $cof_id;
}


sub get_cofs_for_fs
{
	my $self = shift;
	my($fs_id) = @_;
	my $res = ah2a('COF_ID', $self->{bu}->exec_SQL_ar('select cof_id from cof_gtfs where fs_id=?',$fs_id));
	return @$res > 0 ? $res : [];
}


sub get_info_about_cof
{
	my $self = shift;
	my($cof_id) = @_;

	my $info = $self->{bu}->exec_SQL_ar(q[
		SELECT cof_id, name, kingdoms, num_fs, num_core_fs, num_b_fs, num_p_fs, num_r_fs, num_bp_fs, num_genus, num_species, num_jobs,
		       has_paml, dir, fs_dna_num_left_gaps, fs_dna_num_right_gaps, fs_dna_num_inner_gaps, fs_dna_trimmed_len,
		       fs_dna_num_identical_cols, fs_dna_identical_block, c_date, gene_len_left, stop_stop_len, kaks_left, gene_len_full, kaks_full,
			   meme_evalue, meme_width, num_plus_fs, num_minus_fs, subst_per_site, num_down_rbs, down_rbs_score,
			   meme_period_1, meme_period_2, meme_period_3, meme_profile_sum, meme_profile_sum / meme_width AS meme_profile_aver,
		       ss50_word_exp, ss50_word_var, ss50_identity,
			   best_word_score, best_word_num_fs, best_word, best_p_word, best_p_word_num_fs
		FROM cofs WHERE cof_id=?
	], $cof_id)->[0];
	$info->{FS_DNA_IDENTICAL_BLOCK_LEN} = length($info->{FS_DNA_IDENTICAL_BLOCK}) if $info->{FS_DNA_IDENTICAL_BLOCK};

	if( !defined $info->{NUM_FS} )
	{
		$info->{NUM_FS} = scalar( @{$self->get_all_fs_for_cof($cof_id)} );
		$self->{bu}->exec_SQL_nr('UPDATE cofs SET num_fs=? WHERE cof_id=?', $info->{NUM_FS}, $cof_id);
	}
	if( !defined $info->{NUM_CORE_FS} )
	{
		$info->{NUM_CORE_FS} = scalar( @{$self->get_core_fs_for_cof($cof_id)} );
		$self->{bu}->exec_SQL_nr('UPDATE cofs SET num_core_fs=? WHERE cof_id=?', $info->{NUM_CORE_FS}, $cof_id);
	}
	if( !defined $info->{NUM_PLUS_FS} )
	{
		$info->{NUM_PLUS_FS} = $self->{bu}->exec_SQL_ar(q[
			SELECT COUNT(DISTINCT id) AS c FROM fshifts
			WHERE len = 1 AND cof_id=?
		], $cof_id )->[0]{C};
		$self->{bu}->exec_SQL_nr('UPDATE cofs SET num_plus_fs=? WHERE cof_id=?', $info->{NUM_PLUS_FS}, $cof_id);
	}
	if( !defined $info->{NUM_MINUS_FS} )
	{
		$info->{NUM_MINUS_FS} = $self->{bu}->exec_SQL_ar(q[
			SELECT COUNT(DISTINCT id) AS c FROM fshifts
			WHERE len = -1 AND cof_id=?
		], $cof_id )->[0]{C};
		$self->{bu}->exec_SQL_nr('UPDATE cofs SET num_minus_fs=? WHERE cof_id=?', $info->{NUM_MINUS_FS}, $cof_id);
	}

	return $info;
}


###
# The method doesn't return {FS_IDS} anymore, because it takes too much time to get them.
# Use get_all_fs_for_cof($cof_id) method to get them
#
sub get_ALL_info_about_cof
{
	my $self = shift;
	my($cof_id) = @_;

	my $info = $self->get_info_about_cof( $cof_id );

	return $info;
}

sub get_all_fs_for_cof
{
	my $self = shift;
	my( $cof_id ) = @_;
	return ah2a('FS_ID', $self->{bu}->exec_SQL_ar('SELECT DISTINCT fs_id FROM cof_gtfs WHERE cof_id=?', $cof_id) );
}


sub get_core_fs_for_cof
{
	my $self = shift;
	my($cof_id) = @_;
	return ah2a('FS_ID', $self->{bu}->exec_SQL_ar('select fs_id from cof_gtfs where cof_id=? and is_core=1', $cof_id ));
}


sub delete_cof
{
	my $self = shift;
	my($cof_id) = @_;

	my $cof = $self->get_ALL_info_about_cof( $cof_id );
	system("rm -rvf $cof->{DIR}");

	$self->{bu}->exec_SQL_nr('DELETE FROM cof_gtfs WHERE cof_id=?', $cof_id );
#	$self->{bu}->exec_SQL_nr('UPDATE fshifts SET cof_id=NULL WHERE cof_id=?', $cof_id );	# ???
	$self->{bu}->exec_SQL_nr('DELETE FROM cofs WHERE cof_id=?', $cof_id );
}


sub get_all_TP_hit_hsp
{
	my($result, $fs_coord, $tp_shoulder) = @_;

	my $all_tp = [];
	while(my $hit = $result->next_hit)
	{
		my($tp_hsp, $best_dist) = (undef, 0);
		while(my $hsp = $hit->next_hsp)
		{
			if($hsp->start('query') < $fs_coord - $tp_shoulder && $fs_coord + $tp_shoulder < $hsp->end('query'))
			{
				my $dist = min($fs_coord - $hsp->start('query'), $hsp->end('query') - $fs_coord);
				if( $dist > $best_dist )
				{
					$tp_hsp    = $hsp;
					$best_dist = $dist;
				}
			}
		}

		if( $tp_hsp )
		{
			push @$all_tp, {
				hit      => $hit,
				hsp      => $tp_hsp,
				min_dist => $best_dist,
			};
		}
	}

	return $all_tp;
}

1;
