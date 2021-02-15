package MyLib::BaseUtil;

use strict;
use warnings;

# $Id$

##################################################################
# Antonov Ivan                  ver.1.02
#
##################################################################
#
# $self
#   |
#   |--->{mode}             -  'cgi' or 'script'
#   |
#   |--->{dbh}
#   |
#   |--->{data_path}        -  the value from MyLib::Local module
#   |--->{tmp_dir}          -  the value from MyLib::Local module
#   |--->{mylib_path}       -  the value from MyLib::Local module
#   |--->{exe_path}         -  the value from MyLib::Local module
#   |--->{script_path}      -  the value from MyLib::Local module
#   |--->{templates_path}   -  the value from MyLib::Local module
#   |
#   |--->{dir_info}         -  hashref. Information about all folders
#   |        |
#   |        |
#   |        |--->{$dir_key}
#   |               |
#   |               |--->{dir_name}
#   |               |--->{descr}
#   |               |--->{file_fmt}
#   |               |--->{file_ext}
#   |               |--->{ver}
#   |
#   |--->{time_debug}           -  is used in 'exec_SQL_ar_dbg' and 'exec_SQL_nr_dbg' methods ???
#   |--->{time_debug_limit}     -  is used in 'exec_SQL_ar_dbg' and 'exec_SQL_nr_dbg' methods ???
#   |
#   |--->{_temp_var}
#   |--->{_unlink_at_destr}    -  reference to array with filenames
#   |--->{_stdout_fn}
#   |--->{_stderr_fn}
#   |
#   |--->{CACHE}
#   |       |
#   |       |
#

use Data::Dumper;
use Carp qw(cluck confess);
use Time::HiRes qw(gettimeofday tv_interval);
use DBI;

use MyLib::Local;

use Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA       = qw(Exporter);
@EXPORT_OK = ( qw(
	min
	aoh2array
	ah2a
));

#############
# CONSTANTS

my $MODES = {
	'script' => 1,
	'cgi'    => 1,
};

##############
# CONSTRUCTOR

###
# 'debug_sql'  --  exec exec_SQL_nr_dbg and exec_SQL_ar_dbg
sub new
{
	my $class = shift;
	my(%opts) = @_;

	my $mode = $opts{mode} || 'script';
	confess("Unknown mode '$mode'") unless $MODES->{$mode};

	my $self = bless {
		warn_debg        => $opts{warn_debg} || 0,
		mode             => $mode,
		tmp_dir          => ( $opts{tmp_dir} || MyLib::Local::TMP_DIR() ),
		db_name          => $MyLib::Local::DB_ACCOUNTS{default}{db_name},
		db_account       => $opts{db_account},
		debug_sql        => $opts{debug_sql},
		dbi_decode       => $opts{dbi_decode},
		time_debug       => exists $opts{time_debug_limit},
		time_debug_limit => $opts{time_debug_limit},
		_unlink_at_destr => [],
		_temp_var        => 0,
	}, $class;

	$self->{_stdout_fn} = $self->get_tmp_filename(prefix => '_stdout.');
	$self->{_stderr_fn} = $self->get_tmp_filename(prefix => '_stderr.');

	return $self;
}

##################
# PRIVATE METHODS
sub DESTROY
{
	my $self = shift;
	if ( $self->{dbh} )
	{
		$self->{dbh}->rollback unless $self->{dbh}{AutoCommit};

		# With this line it constantly warns "DBD::mysql::db do failed: MySQL server has gone away..."
#		$self->{dbh}->disconnect;
	}

	# Remove temporary files
	foreach my $fn ( @{$self->{_unlink_at_destr}} )
	{
		-d $fn ? rmdir($fn) : unlink($fn);
	}
}

##################################################################
# PUBLIC METHODS
sub dbh
{
	my $self = shift;
	return $self->{dbh} if $self->{dbh};

	$self->{dbh} = DBI->connect( MyLib::Local::DBH_CONNECT($self->{db_account}), MyLib::Local::DBH_ATTR($self->{db_account}) ) ||
		confess "Can't connect to database '$self->{db_account}': $DBI::errstr";

	if( $self->{dbi_decode} && $self->{dbi_decode} eq 'utf8' )
	{
		$self->exec_SQL_nr("SET character_set_database='utf8';");
		$self->exec_SQL_nr("SET character_set_server='utf8';");
		$self->exec_SQL_nr("SET character_set_client='utf8'");
		$self->exec_SQL_nr("SET character_set_connection='utf8'");
		$self->exec_SQL_nr("SET character_set_results='utf8'");
	}

	return $self->{dbh};
}

sub commit
{
	my $self = shift;
	$self->dbh->commit() or confess "Can't commit: ".$DBI::errstr;
}

sub exec_SQL_nr
{
	my $self = shift;
	my $sql  = shift;

	eval
	{
		my $sth = $self->dbh->prepare_cached($sql);
		confess("Can't prepare [$sql]: $DBI::errstr") unless $sth;

		my $res = $sth->execute( @_ );
		confess("Can't execute [$sth->{Statement}]: $DBI::errstr") unless $res;

		$sth->finish;
		1;
	}
	or do
	{
		confess "\n\n--- SQL ERROR ---\n$@\n---\n$sql\n---\n".Data::Dumper->Dump([\@_], ['PARAMS'])."---\n";
	};
}

sub _exec_SQL_ar
{
	my $self = shift;
	my $sql  = shift;

	my $sth = $self->dbh->prepare_cached( $sql );
	confess "Can't prepare [$sql]: $DBI::errstr" unless $sth;

	my $res = $sth->execute( @_ );
	confess "Can't execute [$sth->{Statement}]: $DBI::errstr" unless $res;

	my @array = ();
	while(my $hash = $sth->fetchrow_hashref)
	{
		if( $self->{dbi_decode} )
		{
			$hash->{$_} = decode($self->{dbi_decode}, $hash->{$_}) foreach keys %$hash;
		}
		push(@array, $hash);
	}

	$sth->finish;

	return \@array;
}

sub exec_SQL_ar
{
	my $self = shift;
	return $self->{debug_sql} ? $self->exec_SQL_ar_dbg(@_) : $self->_exec_SQL_ar(@_);
}

sub exec_SQL_ar_dbg
{
	my $self = shift;
	my $sql  = shift;

	my @array = ();
	eval
	{
		# start time
		my $t0 = [gettimeofday];

		@array = @{$self->_exec_SQL_ar($sql, @_)};

		# elapsed time
		if( $self->{time_debug})
		{
			my $elapsed = tv_interval($t0);
			$sql =~ s/[\n\r\t]/ /gs;
			$sql =~ s/^\s*(\S.*\S)\s*$/$1/;
			print "$elapsed\t$sql\n" if $elapsed > $self->{time_debug_limit};
		}

		1;
	}
	or do
	{
		confess "\n\n--- SQL ERROR ---\n$@\n---\n$sql\n---\n".Data::Dumper->Dump([\@_], ['PARAMS'])."---\n";
	};

	return \@array;
}

sub get_tmp_dir
{
	my $self = shift;
	return $self->{tmp_dir};
}

sub get_tmp_filename
{
	my $self = shift;
	my(%opts) = @_;
	my $del  = exists $opts{del}    ? $opts{del}    : 1;
	my $dir  = exists $opts{dir}    ? $opts{dir}    : $self->get_tmp_dir;
	my $ext  = exists $opts{ext}    ? ".$opts{ext}" : '';
	my $pref = exists $opts{prefix} ? $opts{prefix} : '';

	my $fn;
	while(1)
	{
		$fn = $self->{tmp_dir}.$pref.int(rand(1000000)).$ext;
		last unless -e $fn;
	}

	$self->add_file_to_delete_at_destroy($fn) if $del;

	return $fn;
}

sub add_file_to_delete_at_destroy
{
	my $self = shift;
	my(@files) = @_;
	push( @{$self->{_unlink_at_destr}}, @files );
}

sub ah2a{ return aoh2array(@_); }
sub aoh2array
{
	my($key, $aoh) = @_;
	return [ map { $_->{$key} } @$aoh ];
}

#########################
# Mathematical functions

sub min
{
	return undef unless @_;

	my $min = shift;
	return $min unless @_;

	for( @_ )
	{
		$min = $_ if $_ < $min;
	}
	return $min;
}


1;

