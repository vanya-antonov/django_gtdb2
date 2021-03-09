package MyLibGT::Local;

# $Id$

use Carp 'confess';

###########################
# set up $ENV{HOME} below #
###########################

our %DFLT_DBH_ATTR = (
	RaiseError       => 1,
	AutoCommit       => 0,
	FetchHashKeyName => 'NAME_uc',
	LongReadLen      => 1024 * 10,
	LongTruncOk      => 1,
);


our %DB_ACCOUNTS = (
	ivdb => {
		dbi_str  => "dbi:mysql:ivanyaco_toilet:ivanya.com",
		db_user  => "ivanyaco_toilet",
		db_pass  => "123",
		dbh_attr => \%DFLT_DBH_ATTR,
	},
	svdb => {
		dbi_str  => "dbi:mysql:slavanya_deals:localhost",
		db_user  => "slavanya_deals",
		db_pass  => "123",
		dbh_attr => \%DFLT_DBH_ATTR,
	},
	gtdb2_cof => {
		dbi_str  => "dbi:mysql:gtdb2_cof:localhost",
		db_name  => "gtdb2_cof",
		db_user  => "genetack",
		db_pass  => "secret",
		dbh_attr => \%DFLT_DBH_ATTR,
	},
);
$DB_ACCOUNTS{default} = $DB_ACCOUNTS{gtdb2_cof};

BEGIN
{
	$ENV{HOME} ||= '/home/alessandro/';
}

sub TMP_DIR { return '/tmp/'; }

sub DBH_CONNECT
{
	my( $account_name ) = @_;
	$account_name = 'default' if !$account_name;
	my $db = $DB_ACCOUNTS{ $account_name };
	confess("DB Account '$account_name' is not defined!!!") if !defined $db;
	return($db->{dbi_str}, $db->{db_user}, $db->{db_pass});
}

sub DBH_ATTR
{
	my($account_name) = @_;
	$account_name = 'default' if !$account_name;
	my $db = $DB_ACCOUNTS{ $account_name };
	confess("DB Account '$account_name' is not defined!!!") if !defined $db;
	return $db->{dbh_attr};
}

1;

