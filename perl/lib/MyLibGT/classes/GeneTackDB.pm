package MyLibGT::classes::GeneTackDB;

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
#   |--->{tbl}  --  MyLibGT::classes::Table object
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
use MyLibGT::BaseUtil qw(min ah2a);


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

	my $cof_id = $self->{bu}->exec_SQL_ar('SELECT IFNULL(MAX(id)+1,1000000) AS next_id FROM cofs')->[0]{NEXT_ID};
#	my $cof_dir = "$self->{dir}/cofs/$cof_id";
#	`mkdir -p $cof_dir`;

	$self->{bu}->exec_SQL_nr('INSERT INTO cofs(id, user_id, c_date) VALUES(?,?,NOW())', $cof_id, $opts{user_id} );

	foreach my $fs_id ( @$fs_arr )
	{
		$self->{bu}->exec_SQL_nr('UPDATE fshifts SET cof_id=? WHERE id=?', $cof_id, $fs_id );
	}

	my $num_fs = ~~@$fs_arr;

	if( $self->{bu}->exec_SQL_ar('SELECT COUNT(1) AS N FROM cof_params WHERE name="num_fs" AND parent_id=?', $cof_id)->[0]{N} ){
		$self->{bu}->exec_SQL_nr('UPDATE cof_params SET num=? WHERE name="num_fs" AND parent_id=?', $num_fs, $cof_id );
	}else{
		$self->{bu}->exec_SQL_nr('INSERT INTO cof_params(name, num, parent_id) VALUES("num_fs",?,?)', $num_fs, $cof_id );
	}

#	if( $self->{bu}->exec_SQL_ar('SELECT COUNT(1) AS N FROM cof_params WHERE name="dir" AND parent_id=?', $cof_id)->[0]{N} ){
#		$self->{bu}->exec_SQL_nr('UPDATE cof_params SET value=? WHERE name="dir" AND parent_id=?', $cof_dir, $cof_id );
#	}else{
#		$self->{bu}->exec_SQL_nr('INSERT INTO cof_params(name, value, parent_id) VALUES("dir",?,?)', $cof_dir, $cof_id );
#	}

	my $num_plus_fs = $self->{bu}->exec_SQL_ar('SELECT COUNT(DISTINCT id) AS N FROM fshifts WHERE len>0 AND cof_id=?', $cof_id)->[0]{N} || 0;

	if( $self->{bu}->exec_SQL_ar('SELECT COUNT(1) AS N FROM cof_params WHERE name="num_plus_fs" AND parent_id=?', $cof_id)->[0]{N} ){
		$self->{bu}->exec_SQL_nr('UPDATE cof_params SET num=? WHERE name="num_plus_fs" AND parent_id=?', $num_plus_fs, $cof_id );
	}else{
		$self->{bu}->exec_SQL_nr('INSERT INTO cof_params(name, num, parent_id) VALUES("num_plus_fs",?,?)', $num_plus_fs, $cof_id );
	}

#	my $num_minus_fs = $self->{bu}->exec_SQL_ar('SELECT COUNT(DISTINCT id) AS N FROM fshifts WHERE len<0 AND cof_id=?', $cof_id)->[0]{N} || 0;
	my $num_minus_fs = $num_fs - $num_plus_fs;

	if( $self->{bu}->exec_SQL_ar('SELECT COUNT(1) AS N FROM cof_params WHERE name="num_minus_fs" AND parent_id=?', $cof_id)->[0]{N} ){
		$self->{bu}->exec_SQL_nr('UPDATE cof_params SET num=? WHERE name="num_minus_fs" AND parent_id=?', $num_minus_fs, $cof_id );
	}else{
		$self->{bu}->exec_SQL_nr('INSERT INTO cof_params(name, num, parent_id) VALUES("num_minus_fs",?,?)', $num_minus_fs, $cof_id );
	}

	return $cof_id;
}


sub get_cofs_for_fs
{
	my $self = shift;
	my($fs_id) = @_;
	my $res = ah2a('COF_ID', $self->{bu}->exec_SQL_ar('SELECT cof_id FROM fshifts WHERE id=?', $fs_id ));
	return @$res > 0 ? $res : [];
}


sub delete_cof
{
	my $self = shift;
	my( $cof_id ) = @_;
	return unless $cof_id;

#	my $cof_dir = $self->{bu}->exec_SQL_ar('SELECT value AS dir FROM cof_params WHERE name="dir" AND parent_id=? LIMIT 1', $cof_id )->[0]{DIR};
#	system("rm -rvf $cof_dir") if -d $cof_dir;

	$self->{bu}->exec_SQL_nr('UPDATE fshifts SET cof_id=NULL WHERE cof_id=?', $cof_id );
	$self->{bu}->exec_SQL_nr('DELETE FROM cof_params WHERE name="dir" AND parent_id=?', $cof_id );
	$self->{bu}->exec_SQL_nr('DELETE FROM cof_params WHERE name="num_fs" AND parent_id=?', $cof_id );
	$self->{bu}->exec_SQL_nr('DELETE FROM cof_params WHERE name="num_plus_fs" AND parent_id=?', $cof_id );
	$self->{bu}->exec_SQL_nr('DELETE FROM cof_params WHERE name="num_minus_fs" AND parent_id=?', $cof_id );

	# Should be the last!!!
	$self->{bu}->exec_SQL_nr('DELETE FROM cofs WHERE id=?', $cof_id );
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
