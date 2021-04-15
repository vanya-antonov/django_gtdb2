#!/usr/bin/perl -w

use strict;
use warnings;
no strict 'refs';

use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use DBI;

my $VERSION = '1.04';

# BLAST options
my $evalue = '1e-10';
my $num_threads = 14;

# DB for aligment by BLASTp
my $prot_db = '/home/alessandro/BLASTp/db/Streptomyces.seq_prot.faa';

###
# Default Options
my $OUTPUT;
my $SAVE_FNA_SEQS;
my $SAVE_FAA_SEQS;
my $AUTO;
my $WOFS;

###
# Parse input data
GetOptions(
	'output=s'        => \$OUTPUT,
	'save_fna_seqs=s' => \$SAVE_FNA_SEQS,
	'save_faa_seqs=s' => \$SAVE_FAA_SEQS,
	'auto'            => \$AUTO,
	'wofs'            => \$WOFS,
) or die &usage();

my $infile = $ARGV[0] || die &usage(); # 'NZ_CP023977.1.gbk';

# Check GenBank format
open GBF, $infile or die &usage();
$_ = <GBF>;
close GBF;

die "\x1b[31mERROR\x1b[0m: Invalid GenBank file format" unless /^LOCUS\s+/;

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

my( $WOFS_TSV, $WOFS_FNA, $WOFS_FAA );

if( defined $WOFS ){
	(my $fl = $infile) =~s/\.[^\.]+$//;

	$WOFS_TSV = "$fl.without_fs.tsv";
	$WOFS_FNA = "$fl.without_fs.fna";
	$WOFS_FAA = "$fl.without_fs.faa";

	open WOFSTSV, ">$WOFS_TSV";
	open WOFSFNA, ">$WOFS_FNA";
	open WOFSFAA, ">$WOFS_FAA";
}

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
	my( $organism, $taxon, $num_fs_and_TTA_genes, $num_TTA_genes_in_cofs, $num_fs_and_TTA_genes_in_cofs );

	for my $feat_obj ($seq_obj->get_SeqFeatures) {
		my $primary_tag = $feat_obj->primary_tag;

		if( $primary_tag eq 'source'){

			# from a line like /organism="Streptomyces griseus subsp. griseus NBRC 13350"
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
					if( /^taxon:(\d+)/ ){ # from a line like /db_xref="taxon:455632"
						$taxon = $1;
						$taxon =~s/\D+//g;
						last;
					}
				}
			}

			next;
		}

		next if $primary_tag ne 'CDS'; #  && ( $primary_tag ne 'gene');
		my $type_seq = $primary_tag;

		# Gene Location
		my $start = $feat_obj->location->start;
		my $end   = $feat_obj->location->end;
		my $len   = $end - $start;

		# Sequence Identifier
		my $gene_id = "$acc_id:$start.$len";

		next if exists $dt{ $gene_id };

		# Get Nucleotide sequence (CDS)
		my $fna_seq = lc $feat_obj->spliced_seq->seq; # e.g. 'ATTATTTTCGCT...'

		# Search TTA codons
		my @TTAs;
		pos $fna_seq = 0;
		while( $fna_seq=~/tta/g ){
			my $s = $-[0]; # start point of TTA

			next if $s % 3; # Take only ORF TTA coordinate

			push @TTAs, $s;
		}
		next unless @TTAs; # gene/transcript without TTA-codons

		if( $feat_obj->has_tag('protein_id') ){
			# e.g. 'WP_100112605.1', from a line like /protein_id="WP_100112605.1"
			$dt{ $gene_id }{'proteins'} = join ',', $feat_obj->get_tag_values('protein_id'); # protein_id(s) for FNA/FAA
		}

		if( $feat_obj->has_tag('gene') ){
			# e.g. 'NDP', from a line like '/gene="NDP"'
			$dt{ $gene_id }{'genes'} = join ',', $feat_obj->get_tag_values('gene'); # gene name(s) for FNA/FAA

		}elsif( $feat_obj->has_tag('locus_tag') ){
			$dt{ $gene_id }{'genes'} = join ',', $feat_obj->get_tag_values('locus_tag');
		}

		# Save Gene Location
		$dt{ $gene_id }{'start'} = $start;
		$dt{ $gene_id }{'end'}   = $end;

		$dt{ $gene_id }{'TTAs'} = \@TTAs; # Save TTA-codons for the future...

		# Capitalize TTA-codons in sequence (CDS)
		substr( $fna_seq, $_, 3) = 'TTA' for @TTAs;
		$dt{ $gene_id }{'fna'} = $fna_seq;

		# Get protein sequence
		my $yes_cof;
		if( $feat_obj->has_tag('translation') ){
			$dt{ $gene_id }{'faa'} = uc(join '', $feat_obj->get_tag_values('translation'));

			# Aligment by BLASTp and search FS-ortholog and Clasters
			$yes_cof = &BLAST_algn( $dbh, $prot_db, $evalue, $num_threads, $dt{ $gene_id }{'faa'} );
			++$num_TTA_genes_in_cofs if $yes_cof;
		}

		# Search corresponding FS for TTA-genes: |start...{FS-coord}...end|
		my $num = $dbh->selectrow_array( qq{ SELECT COUNT(DISTINCT id) AS num FROM fshifts
WHERE seq_id LIKE "$acc_id" AND $start <= coord AND coord <= $end } );

		if( $num ){
			$num_fs_and_TTA_genes += $num;
			++$num_fs_and_TTA_genes_in_cofs if $yes_cof;

		}elsif( $WOFS ){
			my $strand = $feat_obj->location->strand;

			my $descr = '';
			if( $feat_obj->has_tag('product') ){
				$descr = join ';', $feat_obj->get_tag_values('product'); # from a line like: /product="VWA domain-containing protein"
			}

			print WOFSTSV join("\t", $gene_id, ($dt{ $gene_id }{'genes'}||''), ($dt{ $gene_id }{'proteins'}||''),
									$strand, $type_seq, $descr ), "\n";

			&save_fasta('WOFSFNA', 'fna', $gene_id, \%dt );
			&save_fasta('WOFSFAA', 'faa', $gene_id, \%dt );
		}

	}

	if( $SAVE_FNA_SEQS ){
		# Save nucleotide sequence
		&save_fasta('FFNA', 'fna', $_, \%dt ) for sort{ $dt{$a}{'start'} <=> $dt{$b}{'start'} } keys %dt;
	}

	if( $SAVE_FAA_SEQS ){
		# Save protein (AA) sequence
		&save_fasta('FFAA', 'faa', $_, \%dt ) for sort{ $dt{$a}{'start'} <=> $dt{$b}{'start'} } keys %dt;
	}

	my( $num_fs_genes ) = $dbh->selectrow_array( qq{ SELECT COUNT(*) AS num_fs FROM fshifts WHERE seq_id LIKE "$acc_id" GROUP BY seq_id } );

	my( $num_cofs ) = $dbh->selectrow_array( qq{ SELECT COUNT(DISTINCT cof_id) AS num_cofs FROM fshifts
WHERE cof_id IS NOT NULL AND seq_id LIKE "$acc_id" GROUP BY seq_id } );

	my( $num_fs_genes_in_cofs ) = $dbh->selectrow_array( qq{ SELECT COUNT(cof_id) AS num_fs_cofs FROM fshifts
WHERE cof_id IS NOT NULL AND seq_id LIKE "$acc_id" GROUP BY seq_id } );

	$taxon ||= $species_string;
	$all{ $taxon }{'ORG_NAME'} = $organism; # 1

	$all{ $taxon }{'NUM_TTA_GENES'} += scalar( keys %dt ); # 2
	$all{ $taxon }{'NUM_FS_GENES'} += $num_fs_genes || 0; # 3

	$all{ $taxon }{'NUM_FS_and_TTA_GENES'} += $num_fs_and_TTA_genes || 0; # 4

	$all{ $taxon }{'NUM_COFS'} += $num_cofs || 0; # 5
	$all{ $taxon }{'NUM_FS_GENES_in_COFS'} += $num_fs_genes_in_cofs || 0; # 6

	$all{ $taxon }{'NUM_TTA_GENES_in_COFS'} += $num_TTA_genes_in_cofs || 0; # 7
	$all{ $taxon }{'NUM_FS_and_TTA_GENES_in_COFS'} += $num_fs_and_TTA_genes_in_cofs || 0; # 8

	push @{ $all{ $taxon }{'ACC_IDs'} }, $acc_id; # 9

}

print "\n# Elapsed time: ".(time - $START_TIME)." sec\n";

unless( scalar( keys %all )){

	if( $WOFS ){
		close WOFSTSV;
		close WOFSFNA;
		close WOFSFAA;

		`rm $WOFS_TSV` if -z $WOFS_TSV;
		`rm $WOFS_FNA` if -z $WOFS_FNA;
		`rm $WOFS_FAA` if -z $WOFS_FAA;
	}

	if( $SAVE_FNA_SEQS ){
		close FFNA;
		`rm $SAVE_FNA_SEQS` if -z $SAVE_FNA_SEQS;
	}

	if( $SAVE_FAA_SEQS ){
		close FFAA;
		`rm $SAVE_FAA_SEQS` if -z $SAVE_FAA_SEQS;
	}

	exit;
}

my @hh = ( qw{
	ORG_NAME
	NUM_TTA_GENES
	NUM_FS_GENES
	NUM_FS_and_TTA_GENES
	NUM_COFS
	NUM_FS_GENES_in_COFS
	NUM_TTA_GENES_in_COFS
	NUM_FS_and_TTA_GENES_in_COFS
} );

my $head_out = join "\t", 'TAXON', @hh, 'ACC_IDs';

# Save collection of TTA-genes
if( $OUTPUT ){
	open OFILE, ">$OUTPUT";
	print OFILE "$head_out\n";
}else{
	print "$head_out\n";
}

for my $taxon ( sort keys %all ){
	my $out = join("\t", $taxon, @{ $all{ $taxon } }{ @hh }, join(';', @{ $all{ $taxon }{'ACC_IDs'} } ) );

	if( $OUTPUT ){
		print OFILE "$out\n";
	}else{
		print "$out\n";
	}

}

exit;


sub save_fasta {
	my( $fh, $stype, $head, $dt ) = @_;
	return unless exists $dt->{ $head }{ $stype };

	my @h;
	push @h, ">$head";
	push @h, "gene=$dt->{ $head }{'genes'}"       if $dt->{ $head }{'genes'};
	push @h, "protein=$dt->{ $head }{'proteins'}" if $dt->{ $head }{'proteins'};

	print $fh join(' ', @h), "\n$dt->{ $head }{ $stype }\n";
}


# Aligment by BLASTp: Search FS-ortholog and Clasters
sub BLAST_algn {
	my( $dbh, $prot_db, $evalue, $num_threads, $fasta ) = @_;
	my $num;

	# -max_hsps 1
	my $run = qq{echo $fasta | blastp -db $prot_db -outfmt "6 sacc" -num_threads $num_threads -evalue $evalue};

	open( HITSF, "$run 2>/dev/null |") or die "Can't run BLASTp: $!";
	while(<HITSF>){
		chomp;
		next if /^$/ or /^\D/;	# Skip empty line or barewords
=comment
NC_010572.1:66749.1145	6079	37.596	391	213	10	13	381	74	455	2.45e-64	217

my( $gene_id, $fs_id, $pident, $alen, $mismatches, $gaps, $qstart, $qend, $sstart, $send, $evalue, $bitscore ) = split /\t/;
=cut
		$num = $dbh->selectrow_array( qq{ SELECT COUNT(cof_id) FROM fshifts WHERE cof_id IS NOT NULL AND id=$_} ); # $fs_id

		last if $num > 0; # Take 1st best COF
	}
	close HITSF;

	return $num;
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
OR
    $script NC_010572.1.gbk -save_fna_seq NC_010572.1.fna -save_faa_seq NC_010572.1.faa -output NC_010572.1.tsv

HERE:
    <file.gbk>   -- input GenBank file only

OPTIONS:
    --output  <file.tsv|stdout>  --  output table. Default STDOUT output
    --save_fna_seqs  <file.fna>  --  save nt sequences of CDS(s) with TTA
    --save_faa_seqs  <file.faa>  --  save sequences of all protein(s) with TTA (Leu)
    --auto                       --  autocomplete options: --output, --save_fna_seqs, --save_faa_seqs
    --wofs                       --  save the collection of TTA-genes and their sequences without FrameShift

OUTPUT TABLE FORMAT:
  0.TAXON                        --  NCBI taxon id of organism
  1.ORG_NAME                     --  organism name
  2.NUM_TTA_GENES                --  total number of genes with TTA codon (TTA-genes)
  3.NUM_FS_GENES                 --  total number of FS-genes
  4.NUM_FS_and_TTA_GENES         --  intersection of FS-genes with TTA-genes
  5.NUM_COFS                     --  total number of clusters that include FS-genes
  6.NUM_FS_GENES_in_COFS         --  total number of FS-genes in clusters
  7.NUM_TTA_GENES_in_COFS        --  total number of TTA-genes that are 'similar' to FS-genes from ALL clusters
  8.NUM_FS_and_TTA_GENES_in_COFS --  
  9.ACC_IDs                      --  Accession ID(;s) of locus/genomic sequence(s)

NOTES:
  1. a configuration DB file 'db.cfg' is required.
  2. BLAST program is required.

";
}

