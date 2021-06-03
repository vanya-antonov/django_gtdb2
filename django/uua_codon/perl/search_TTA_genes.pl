#!/usr/bin/perl -w

use strict;
use warnings;
no strict 'refs';

use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use DBI;

my $VERSION = '1.08';

###
# Default Options
my $OUTPUT;
my $SAVE_FNA;
my $SAVE_FAA;
my $SAVE_TTA_FNA;
my $SAVE_TTA_FAA;
my $AUTO;
my $WOFS;
my $HEADER;
my $ECHO;
my $SKIP_FS;
my $EVALUE = 1e-10; # BLAST option
my $PROT_DB = '/home/alessandro/BLASTp/db/Streptomyces.seq_prot.faa'; # DB for aligment by BLASTp

# Number of CPUs (for BLAST option)
(my $CPUs = `cat /proc/cpuinfo | grep ^processor | wc -l`) =~s/\D+//g;
$CPUs += 0; $CPUs ||= 1;

my $THREADS = int(0.9*$CPUs + 1/2) || 1;

###
# Parse input data
GetOptions(
	'output=s'       => \$OUTPUT,
	'save_fna=s'     => \$SAVE_FNA,	# TODO
	'save_faa=s'     => \$SAVE_FAA,	# TODO
	'save_tta_fna=s' => \$SAVE_TTA_FNA,
	'save_tta_faa=s' => \$SAVE_TTA_FAA,
	'auto'           => \$AUTO,
	'wofs'           => \$WOFS,
	'header'         => \$HEADER,
	'echo'           => \$ECHO,
	'skip_fs'        => \$SKIP_FS,
	'evalue=f'       => \$EVALUE,
	'threads=i'      => \$THREADS,
	'prot_db=s'      => \$PROT_DB,
) or die &usage();

$THREADS = $CPUs if $THREADS > $CPUs;

my $infile = $ARGV[0] || die &usage(); # NC_010572.1.gbk

# Check GenBank format
open GBF, $infile or die &usage();
$_ = <GBF>;
close GBF;

die "\x1b[31mERROR\x1b[0m: Invalid GenBank file format" unless /^LOCUS\s+/;

# Read DataBase configuration
our $dbh;
unless( $SKIP_FS ){
	my $db_cfg = 'db.cfg';
	die "\x1b[31mERROR\x1b[0m: Can't open $db_cfg file" unless -s $db_cfg;
	do "./$db_cfg";
}

if( defined $AUTO ){
	(my $fl = $infile) =~s/\.[^\.]+$//;

	$OUTPUT ||= "$fl.TTA.tsv";
	$SAVE_TTA_FNA ||= "$fl.TTA.fna";
	$SAVE_TTA_FAA ||= "$fl.TTA.faa";
}

$OUTPUT = undef if $OUTPUT && $OUTPUT =~/^stdout$/i;

###
my $START_TIME = time;

$SAVE_FNA && open FFNA, ">$SAVE_FNA";
$SAVE_FAA && open FFAA, ">$SAVE_FAA";
$SAVE_TTA_FNA && open TTAFNA, ">$SAVE_TTA_FNA";
$SAVE_TTA_FAA && open TTAFAA, ">$SAVE_TTA_FAA";

my( $WOFS_TSV, $WOFS_FNA, $WOFS_FAA );

if( defined( $WOFS ) && ! defined( $SKIP_FS ) ){
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

	my %gg; # for all CDS-genes
	my %dt; # for TTA-genes only
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

		my $strand = $feat_obj->location->strand || 0;
		my $ts;
		if( $strand > 0 ){
			$strand = 1;
			$ts = 'p';
		}elsif( $strand < 0 ){
			$strand = -1;
			$ts = 'm';
		}else{
			$ts = 'n';
		}

		# Gene Identifier for FNA/FAA sequence
		my $gene_id = "$acc_id:$ts$start.$len";

		next if exists $gg{ $gene_id };
		$gg{ $gene_id } = undef;

		# DB search corresponding FS for gene (acc_id) and a specific strand: |start...{FS-coord}...end|
		my( $numFS, $fs_ids ) = $SKIP_FS ? (undef, undef) :
								$dbh->selectrow_array( qq{ SELECT COUNT(DISTINCT id) AS numFS, GROUP_CONCAT(DISTINCT id SEPARATOR ';') AS fs_ids
FROM fshifts WHERE seq_id LIKE "$acc_id" AND strand = $strand AND $start <= coord AND coord <= $end } );

		$gg{ $gene_id } = $fs_ids if $numFS;

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

### TTA-codon(s) exists

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
		my( $yes_cof, $cof_id );
		if( $feat_obj->has_tag('translation') ){
			$dt{ $gene_id }{'faa'} = uc(join '', $feat_obj->get_tag_values('translation'));

			# Aligment by BLASTp and search FS-ortholog and Clasters
			( $yes_cof, $cof_id ) = &BLAST_algn( $dbh, $PROT_DB, $EVALUE, $THREADS, $dt{ $gene_id }{'faa'} ) unless $SKIP_FS;
			if( $yes_cof ){
				++$num_TTA_genes_in_cofs;
				push @{ $dt{ $gene_id }{'cofs'} }, $cof_id;
			}
		}

		if( $numFS ){ # corresponding FS for TTA-genes
			$num_fs_and_TTA_genes += $numFS;
			++$num_fs_and_TTA_genes_in_cofs if $yes_cof;

			push @{ $dt{ $gene_id }{'fs_ids'} }, split(';', $fs_ids );

		}else{ # without FS
			push @{ $dt{ $gene_id }{'wofs'} }, $cof_id if $yes_cof;

			if( defined( $WOFS ) && ! defined( $SKIP_FS ) ){

				# from a line like: /product="VWA domain-containing protein"
				my $descr = $feat_obj->has_tag('product') ? join(';', $feat_obj->get_tag_values('product')) : '';

				print WOFSTSV join("\t", $gene_id, ($dt{ $gene_id }{'genes'}||''), ($dt{ $gene_id }{'proteins'}||''),
									$strand, $type_seq, $descr ), "\n";

				&save_fasta('WOFSFNA', 'fna', $gene_id, \%dt );
				&save_fasta('WOFSFAA', 'faa', $gene_id, \%dt );
			}
		}

	}

	if( $SAVE_TTA_FNA ){
		# Save nucleotide sequence
		&save_fasta('TTAFNA', 'fna', $_, \%dt ) for sort{ $dt{$a}{'start'} <=> $dt{$b}{'start'} } keys %dt;
	}

	if( $SAVE_TTA_FAA ){
		# Save protein (AA) sequence
		&save_fasta('TTAFAA', 'faa', $_, \%dt ) for sort{ $dt{$a}{'start'} <=> $dt{$b}{'start'} } keys %dt;
	}

	my( $num_fs_genes ) = $SKIP_FS ? undef :
					$dbh->selectrow_array( qq{ SELECT COUNT(*) AS num_fs FROM fshifts WHERE seq_id LIKE "$acc_id" GROUP BY seq_id } );

	my( $num_cofs ) = $SKIP_FS ? undef :
					$dbh->selectrow_array( qq{ SELECT COUNT(DISTINCT cof_id) AS num_cofs FROM fshifts
WHERE cof_id IS NOT NULL AND seq_id LIKE "$acc_id" GROUP BY seq_id } );

	my( $num_fs_genes_in_cofs ) = $SKIP_FS ? undef :
					$dbh->selectrow_array( qq{ SELECT COUNT(cof_id) AS num_fs_cofs FROM fshifts
WHERE cof_id IS NOT NULL AND seq_id LIKE "$acc_id" GROUP BY seq_id } );

# 1
	$taxon ||= $species_string;
# 2
	my( $org_id ) = $SKIP_FS ? undef : $dbh->selectrow_array( qq{ SELECT org_id FROM seqs WHERE id LIKE "$acc_id"} );
	$all{ $taxon }{'IDs'}{ $org_id } = undef if $org_id;
# 3
	$all{ $taxon }{'ORG_NAME'} = $organism;
# 4
	$all{ $taxon }{'NUM_GENES'} += scalar( keys %gg );
# 5
	$all{ $taxon }{'NUM_TTA_GENES'} += scalar( keys %dt );
# 6
	$all{ $taxon }{'NUM_FS_GENES'} += $num_fs_genes || 0;
# 7
	$all{ $taxon }{'NUM_FS_and_TTA_GENES'} += $num_fs_and_TTA_genes || 0;
# 8
	$all{ $taxon }{'NUM_COFS'} += $num_cofs || 0;
# 9
	$all{ $taxon }{'NUM_FS_GENES_in_COFS'} += $num_fs_genes_in_cofs || 0;
# 10
	$all{ $taxon }{'NUM_TTA_GENES_in_COFS'} += $num_TTA_genes_in_cofs || 0;
# 11
	$all{ $taxon }{'NUM_FS_and_TTA_GENES_in_COFS'} += $num_fs_and_TTA_genes_in_cofs || 0;
# 12
	push @{ $all{ $taxon }{'ACC_IDs'} }, $acc_id;

	for my $gene_id ( keys %dt ){

# 13: for List of Cluster_ID=number_gene(;s) with TTA codon, e.g. 1000515=3;1000517=1;...
		++$all{ $taxon }{'COF_IDs'}{$_} for @{ $dt{ $gene_id }{'cofs'} };
# 14: for List of Frameshift-TTA-gene ID(;s). Format: <fs_id>=<gene_id1,gene_id2,...>,
		push @{ $all{ $taxon }{'FS_IDs'}{$_} }, "$gene_id.3" for @{ $dt{ $gene_id }{'fs_ids'} };
# 15: for List of Cluster ID(;s) without FS-genes. Format: <cluster_id>=<gene_id1,gene_id2,...>,
		push @{ $all{ $taxon }{'WOFS_IDs'}{$_} }, "$gene_id.1" for @{ $dt{ $gene_id }{'wofs'} };

	}

# 16: for List of all CDS gene ID(;s). Format: <acc_id>:<strand><start>.<length>.<tag>
	for( sort keys %gg ){
		my $gtag = 0;
		if( exists $dt{$_} ){
			$gtag = 1; # tag_{TTA} = 1
			$gtag += 2 if exists $dt{$_}{'fs_ids'}; # tag_{TTA+FS} = 3

		}elsif( defined $gg{$_} ){ # tag_{FS} = 2;
			$gtag = 2;
		}
		push @{ $all{ $taxon }{'GENE_IDs'} }, "$_.$gtag";
	}

}

print "\n# Elapsed time: ".(time - $START_TIME)." sec\n" if $ECHO;

unless( scalar( keys %all )){

	if( defined( $WOFS ) && ! defined( $SKIP_FS ) ){
		close WOFSTSV;
		close WOFSFNA;
		close WOFSFAA;

		`rm $WOFS_TSV` if -z $WOFS_TSV;
		`rm $WOFS_FNA` if -z $WOFS_FNA;
		`rm $WOFS_FAA` if -z $WOFS_FAA;
	}

	if( $SAVE_FNA ){
		close FFNA;
		`rm $SAVE_FNA` if -z $SAVE_FNA;
	}

	if( $SAVE_FAA ){
		close FFAA;
		`rm $SAVE_FAA` if -z $SAVE_FAA;
	}

	if( $SAVE_TTA_FNA ){
		close TTAFNA;
		`rm $SAVE_TTA_FNA` if -z $SAVE_TTA_FNA;
	}

	if( $SAVE_TTA_FAA ){
		close TTAFAA;
		`rm $SAVE_TTA_FAA` if -z $SAVE_TTA_FAA;
	}

	exit;
}

my @hh = ( qw{
	ORG_ID
	ORG_NAME
	NUM_GENES
	NUM_TTA_GENES
	NUM_FS_GENES
	NUM_FS_and_TTA_GENES
	NUM_COFS
	NUM_FS_GENES_in_COFS
	NUM_TTA_GENES_in_COFS
	NUM_FS_and_TTA_GENES_in_COFS
} );

my $head_out = join "\t", 'TAXON', @hh, qw(ACC_IDs COF_IDs FS_IDs WOFS_IDs GENE_IDs);

# Save collection of TTA-genes
if( $OUTPUT ){
	open OFILE, ">$OUTPUT";
	print OFILE "$head_out\n" if $HEADER;
}else{
	print "$head_out\n" if $HEADER;
}

for my $taxon ( sort keys %all ){
	$all{ $taxon }{'ORG_ID'} = join ';', sort{ $a <=> $b } keys %{ $all{ $taxon }{'IDs'} };

	my $out = join("\t", $taxon, @{ $all{ $taxon } }{ @hh },
					join(';', @{ $all{ $taxon }{'ACC_IDs'} } ),
					join(';', map{"$_=$all{ $taxon }{'COF_IDs'}{$_}"} sort{ $a <=> $b } keys %{ $all{ $taxon }{'COF_IDs'} } ),
					join(';', map{"$_=". join(',', @{ $all{ $taxon }{'FS_IDs'}{$_} }) } sort{ $a <=> $b } keys %{ $all{ $taxon }{'FS_IDs'} } ),
					join(';', map{"$_=". join(',', @{ $all{ $taxon }{'WOFS_IDs'}{$_} }) } sort{ $a <=> $b } keys %{ $all{ $taxon }{'WOFS_IDs'} } ),
					join(';', @{ $all{ $taxon }{'GENE_IDs'} } ),
				);

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

	my $gtag = 1; # tag_{TTA} = 1
	$gtag += 2 if exists $dt->{$_}{'fs_ids'}; # tag_{TTA+FS} = 3

	my @h = (">$head.$gtag");	# gtag: 0=ordinary, 1=TTA, 2=FS, 3=TTA+FS
	push @h, "gene=$dt->{ $head }{'genes'}"       if $dt->{ $head }{'genes'};
	push @h, "protein=$dt->{ $head }{'proteins'}" if $dt->{ $head }{'proteins'};

	print $fh join(' ', @h), "\n$dt->{ $head }{ $stype }\n";
}


# Aligment by BLASTp: Search FS-ortholog and Clasters
sub BLAST_algn {
	my( $dbh, $prot_db, $evalue, $num_threads, $fasta ) = @_;
	my( $num, $cof_id );

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
		($num, $cof_id) = $dbh->selectrow_array( qq{ SELECT COUNT(cof_id) AS num, cof_id FROM fshifts WHERE cof_id IS NOT NULL AND id=$_ LIMIT 1} ); # $fs_id

		last if $num > 0; # Take 1st best COF
	}
	close HITSF;

	return( $num, $cof_id);
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
    $script NC_010572.1.gbk -save_tta_fna NC_010572.1.TTA.fna -save_tta_faa NC_010572.1.TTA.faa -output NC_010572.1.TTA.tsv

HERE:
    <file.gbk>   -- input GenBank file only

OPTIONS:
    --output  <file.tsv|stdout>  --  output table. Default STDOUT output
    --save_tta_fna  <file.fna>   --  save nt sequences of CDS(s) with TTA
    --save_tta_faa  <file.faa>   --  save sequences of all protein(s) with TTA (Leu)
    --save_fna  <file.fna>       --  save nt sequences of CDS(s). TODO
    --save_faa  <file.faa>       --  save sequences of all protein(s). TODO
    --header                     --  output header of table
    --echo                       --  output echo messages
    --auto                       --  autocomplete options: --output, --save_tta_fna, --save_tta_faa
    --wofs                       --  save the collection of TTA-genes and their sequences without FrameShift
    --evalue                     --  E-value BLAST option. Default 1e-10
    --threads                    --  num_threads BLAST option. Default (0.9*CPUs | 1)
    --prot_db                    --  protein DB for aligment by BLASTp. Default internal Streptomyces.seq_prot.faa DB
    --skip_fs                    --  skip FrameShift search

OUTPUT TABLE FORMAT:
  1.TAXON                        --  NCBI taxon ID of organism
  2.ORG_ID                       --  internal organism ID
  3.ORG_NAME                     --  organism name
  4.NUM_GENES                    --  total number of CDS genes for organism
  5.NUM_TTA_GENES                --  total number of genes with TTA codon (TTA-genes)
  6.NUM_FS_GENES                 --  total number of FS-genes. Empty for --skip_fs mode
  7.NUM_FS_and_TTA_GENES         --  intersection of FS-genes with TTA-genes
  8.NUM_COFS                     --  total number of clusters that include FS-genes
  9.NUM_FS_GENES_in_COFS         --  total number of FS-genes in clusters
 10.NUM_TTA_GENES_in_COFS        --  total number of TTA-genes that are 'similar' to FS-genes from ALL clusters
 11.NUM_FS_and_TTA_GENES_in_COFS --  intersection of FS-genes with TTA-genes in clusters
 12.ACC_IDs                      --  Accession ID(;s) of locus/genomic sequence(s), e.g. NC_003155.5;NC_004719.1;...
 13.COF_IDs                      --  List of Cluster_ID=number_gene(;s) with TTA codon, e.g. 1000515=3;1000517=1;...
 14.FS_IDs                       --  List of Frameshift-TTA-gene ID(;s). Format: <fs_id>=<gene_id1,gene_id2,...>,
                                       e.g. 74297=NC_003155.5:p25699.4695.3,NC_003155.5:m28699.1078.3;...
 15.WOFS_IDs                     --  List of Cluster ID(;s) without FS-genes. Format: <cluster_id>=<gene_id1,gene_id2,...>,
                                       e.g.: 1000568=NC_010572.1:m66749.1145.1,NC_010572.1:p8478036.1145.1;...
 16.GENE_IDs                     --  List of all CDS gene ID(;s). Format: <acc_id>:<strand><start>.<length>.<tag>,
                                       e.g.: NC_003155.5:m869.1085.0;NC_010572.1:m66749.1145.1;NC_010572.1:p8478036.1145.1;...

  Fields ( ORG_ID, NUM_FS_GENES, NUM_FS_and_TTA_GENES, NUM_COFS, NUM_FS_GENES_in_COFS,
         NUM_TTA_GENES_in_COFS, NUM_FS_and_TTA_GENES_in_COFS, COF_IDs, FS_IDs, WOFS_IDs ) are (empty | 0) for --skip_fs option

NOTES (without --skip_fs option):
  1. a configuration DB file 'db.cfg' is required.
  2. BLAST program is required.

";
}
