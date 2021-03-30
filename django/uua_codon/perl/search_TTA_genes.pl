#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use DBI;

my $VERSION = '1.03';

###
# Default Options
my $OUTPUT;
my $SAVE_FNA_SEQS;
my $SAVE_FAA_SEQS;
my $AUTO;

###
# Parse input data
GetOptions(
	'output=s'        => \$OUTPUT,
	'save_fna_seqs=s' => \$SAVE_FNA_SEQS,
	'save_faa_seqs=s' => \$SAVE_FAA_SEQS,
	'auto'            => \$AUTO,
) or die &usage();

my $infile = $ARGV[0] || die &usage(); # 'NZ_CP023977.1.gbk';

# Check GenBank format
open GBF, $infile or die &usage();
my $chk = 0;
while( <GBF> ){
	++$chk if /^(?:LOCUS|DEFINITION|ACCESSION)\s+/;
	last if $. >= 3;
}
close GBF;

die "\x1b[31mERROR\x1b[0m: Invalid GenBank file format" if $chk < 3;

# Read DataBase configuration
my $db_cfg = 'db.cfg';
die "\x1b[31mERROR\x1b[0m: Can't open $db_cfg file" unless -s $db_cfg;
our $dbh; do "./$db_cfg";

if( defined $AUTO ){
	(my $fl = $infile) =~s/\.[^\.]+$//;
	$OUTPUT ||= "$fl.tsv";
	$SAVE_FNA_SEQS ||= "$fl.fna";
	$SAVE_FAA_SEQS ||= "$fl.faa";
}

$OUTPUT = undef if $OUTPUT && $OUTPUT =~/^stdout$/i;

###
my $START_TIME = time;

$SAVE_FNA_SEQS && open FFNA, ">$SAVE_FNA_SEQS";
$SAVE_FAA_SEQS && open FFAA, ">$SAVE_FAA_SEQS";

my $head_out = join "\t", 
	qw{ ACC_ID_SEQ
	ORG_NAME
	NUM_TTA_GENES
	NUM_FS_GENES
	NUM_FS_AND_TTA_GENES
	NUM_COFs
	NUM_FS_GENES_IN_COFS
	NUM_TTA_GENES_IN_COFS
	NUM_FS_AND_TTA_GENES_IN_COFS
	};

if( $OUTPUT ){
	open OFILE, ">$OUTPUT";
	print OFILE "$head_out\n";
}else{
	print "$head_out\n";
}



my $seqio_obj = Bio::SeqIO->new(
	-file   => $infile,
	-format => 'genbank'
	);

while( my $seq_obj = $seqio_obj->next_seq ){

	my $acc_id = $seq_obj->accession;
	my $species_string = $seq_obj->species->node_name; # ORGANISM

	my %dt;
	my $organism;
	for my $feat_obj ($seq_obj->get_SeqFeatures) {
		my $primary_tag = $feat_obj->primary_tag;

#	if( $primary_tag eq 'source'){
#		$organism = join '', $feat_obj->get_tag_values('organism') if $feat_obj->has_tag('organism');
#		next;
#	}

		next if $primary_tag ne 'CDS'; #  && ( $primary_tag ne 'gene');

		my @pp; # protein_id for FNA
		if( $feat_obj->has_tag('protein_id') ){
			push @pp, $feat_obj->get_tag_values('protein_id'); # e.g. 'WP_100112605.1', from a line like /protein_id="WP_100112605.1"
		}

		my @gg; # gene name(s) for FNA
		if( $feat_obj->has_tag('gene') ){
			push @gg, $feat_obj->get_tag_values('gene');	# e.g. 'NDP', from a line like '/gene="NDP"'

		}elsif( $feat_obj->has_tag('locus_tag') ){
			push @gg, $feat_obj->get_tag_values('locus_tag');
		}

		my $head = shift @pp;
		if( $head ){
			next if exists $dt{ $head };

		}elsif( ~~@gg ){
			$head = $gg[0];
			next if exists( $dt{ $head } ) and $primary_tag eq 'gene';

		}else{
			next;
		}

		# Nucleotide sequence (CDS)
		my $fna_seq = lc $feat_obj->spliced_seq->seq; # e.g. 'ATTATTTTCGCT...'

		# Search TTA codons
		my @TTAs;
		pos $fna_seq = 0;
		while( $fna_seq=~/tta/g ){
			my $s = $-[0]; # start point of TTA
#			my $e = $+[0]; # end point

			next if $s % 3; # Take only ORF TTA coordinate

			push @TTAs, $s;
		}
		next unless @TTAs;

		$dt{ $head }{'TTAs'} = \@TTAs;

		# Capitalize TTA-codons in sequence (CDS)
		substr( $fna_seq, $_, 3 ) = 'TTA' for @TTAs;
		$dt{ $head }{'fna'} = $fna_seq;

		$dt{ $head }{'genes'} = join ',', @gg if @gg;

		# Save protein sequence
		$dt{ $head }{'faa'} = uc(join '', $feat_obj->get_tag_values('translation')) if $feat_obj->has_tag('translation');

		# Location
		$dt{ $head }{'start'} = $feat_obj->location->start;
		$dt{ $head }{'end'}   = $feat_obj->location->end;

	}

	my( $num_fs_genes ) = $dbh->selectrow_array( qq{ SELECT COUNT(*) AS num_fs FROM fshifts WHERE seq_id RLIKE "$acc_id" GROUP BY seq_id } );

	my( $num_cofs ) = $dbh->selectrow_array( qq{ SELECT COUNT(DISTINCT cof_id) AS num_cofs FROM fshifts
WHERE cof_id IS NOT NULL AND seq_id RLIKE "$acc_id" GROUP BY seq_id } );

	my( $num_fs_genes_in_cofs ) = $dbh->selectrow_array( qq{ SELECT COUNT(cof_id) AS num_fs_cofs FROM fshifts
WHERE cof_id IS NOT NULL AND seq_id RLIKE "$acc_id" GROUP BY seq_id } );


#ACC_ID_SEQ, ORG_NAME, NUM_TTA_GENES, NUM_FS_GENES, NUM_FS_AND_TTA_GENES, NUM_COFs, NUM_FS_GENES_IN_COFS, NUM_TTA_GENES_IN_COFS, NUM_FS_AND_TTA_GENES_IN_COFS
#1           2         3              4             .                     6         7                     .                      .                          
	my $out = join("\t", $acc_id, $species_string, scalar(keys %dt), $num_fs_genes, '', 
			$num_cofs, $num_fs_genes_in_cofs  );

	if( $OUTPUT ){
		print OFILE "$out\n";
	}else{
		print "$out\n";
	}

	if( $SAVE_FNA_SEQS ){
		for( sort{ $dt{$a}{'start'} <=> $dt{$b}{'start'} } keys %dt ){
			print FFNA ">$_ gene=$dt{$_}{'genes'} location=$dt{$_}{'start'}..$dt{$_}{'end'}\n$dt{$_}{'fna'}\n";
		}
	}

	if( $SAVE_FAA_SEQS ){
		for( sort{ $dt{$a}{'start'} <=> $dt{$b}{'start'} } keys %dt ){
			next unless exists $dt{$_}{'faa'};
			print FFAA ">$_ gene=$dt{$_}{'genes'} location=$dt{$_}{'start'}..$dt{$_}{'end'}\n$dt{$_}{'faa'}\n";
		}
	}

}

print "\n# Elapsed time: ".(time - $START_TIME)." sec\n";
exit;


sub usage
{
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = File::Spec->splitpath($0);
	return"$msg
$0 version $VERSION

USAGE:
    $0 <file.gbk> [OPTIONS]

EXAMPLE:
    $0 NZ_CP023977.1.gbk -auto -output=stdout

    $0 NZ_CP023977.1.gbk -auto
OR
    $0 NZ_CP023977.1.gbk -save_fna_seq NZ_CP023977.1.fna -save_faa_seq NZ_CP023977.1.faa -output NZ_CP023977.1.tsv

HERE:
    <file.gbk>   -- input GenBank file only

OPTIONS:
    --output  <file.tsv|stdout>  --  output table. Default STDOUT output
    --save_fna_seqs  <file.fna>  --  save nt sequences of CDS(s) with TTA
    --save_faa_seqs  <file.faa>  --  save sequences of all protein(s) with TTA (Leu)
    --auto                       --  autocomplete options: --output, --save_fna_seqs, --save_faa_seqs

NOTES:
    A configuration DB file 'db.cfg' is required

";
}
