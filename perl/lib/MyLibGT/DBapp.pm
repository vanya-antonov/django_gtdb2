package MyLibGT::DBapp;

use 5.008008;
use strict;
use warnings;

use Carp;
use DBI;
use JSON;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw(
		parse_gid
		SaveDB
		read_configDB
	) ],
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw( );

our $VERSION = '0.01';

use base qw(Exporter);


sub parse_gid {
	my( $gid ) = @_;

	my( $acc, $location ) = split ':', $gid;	# NC_003155.5:p9004239.817.1
	my( $strand, $sloc, $sz, $gtag ) = $location=~/(p|m)(\d+)\.([\d\-]+)\.(\d+)$/i;
	my $eloc = $sloc + $sz;
	$strand = ($strand=~/p/i) ? +1 : -1;
	my $f_TTA = ($gtag & 0b01 ) ? 1 : 0; # TTA-gene

	return( $acc, $strand, $sloc, $eloc, $gtag, $f_TTA );
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


sub read_configDB
{
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

	my $db_name  = $ref->{'DATABASES'}{ $session }{'NAME'}   || 'gtdb2'; # 'gtdb2_cof'
	my $user     = $ref->{'DATABASES'}{ $session }{'USER'};
	my $password = $ref->{'DATABASES'}{ $session }{'PASSWORD'};
	my $host     = $ref->{'DATABASES'}{ $session }{'HOST'}   || 'localhost';
	my $engine   = $ref->{'DATABASES'}{ $session }{'DBI'}    || 'DBI:mysql:database';
	my $port     = $ref->{'DATABASES'}{ $session }{'PORT'}   || 3306;

	my $dbh = DBI->connect("$engine=$db_name;host=$host;port=$port", $user, $password,
						{ RaiseError => 0, PrintError => 1, AutoCommit => 1} );

	my $gtdb_dir = $ref->{'DATABASES'}{ $session }{'GTDB_DIR'}; # "/home/gtdb/data/gtdb2"
	$gtdb_dir =~s/\/+$//;

	return( $dbh, $gtdb_dir );
}


1;

__END__

=head1 NAME

MyLibGT::DBapp - Some applications to GTDB2 database.

=head1 SYNOPSIS

  use MyLibGT::DBapp

=head1 AUTHOR

Alessandro Gorohovski, E<lt>an.gorohovski@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2020-2021 by A. N. Gorohovski

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.

=cut

