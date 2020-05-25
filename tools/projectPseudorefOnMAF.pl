#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
################################################################

#First pass script to project pseudoreference-based genotypes
# from one reference onto another via pairwise MAF of 1:1
# alignments between references produced by LAST

my $SCRIPTNAME = "projectPseudorefOnMAF.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

projectPseudorefOnMAF.pl - Project pseudoref SNPs on pairwise 1:1 MAF

=head1 SYNOPSIS

projectPseudorefOnMAF.pl [options] <Species 1 prefix> <Species 2 prefix>

 Options:
  --help,-h,-?          Print this help documentation
  --input_MAF,-i        Path to input MAF 1:1 alignment file (default: STDIN)
  --query_FASTA,-q      Path to FASTA of reference used as query by LAST
  --pseudoref,-p        Path to pseudoreference FASTA in query space
  --prefix,-o           Prefix for output files (required)
  --debug,-d            Output debugging information to STDERR
  --version,-v          Output version string

=head1 DESCRIPTION

This script works similarly to groundTruthFromMAF.pl, except it
projects the genotypes encoded in a pseudoref in 
This script generates INSNP files of the ground truth variant calls
in each species' coordinate space based on a MAF file of 1:1 pairwise
alignments produced by LAST. The scaffold IDs in the MAF file should
follow the format [species prefix].[scaffold name], where the species
prefix and scaffold name both are alphanumeric (plus underscores). As
long as the scaffold IDs follow this convention, and you provide the
species prefixes precisely as found in the MAF, the script will work
as intended.
The script will also generate BED files to indicate the extents of each
aligned segment in each appropriate coordinate space.

For example:

=begin text

a score=10 mismap=1e-05
s DyakTai18E2.2L 0 10 + 29953808 AATAACGGCT
s DyakNY73PB.2L  0 10 + 24234981 AACAACGGGT
p                                !!!!!!!!!!
p                                !!!!!!!!!!

=end text

should output in DyakTai18E2_INSNP.tsv:

=begin text

2L	3	T	C
2L	9	C	G

=end text

and in DyakNY73PB_INSNP.tsv:

=begin text

2L	3	C	T
2L	9	G	C

=end text

as well as the following BED files:

DyakTai18E2_aligned_regions.bed:

=begin text

2L	0	10

=end text

DyakNY73PB_aligned_regions.bed:

=begin text

2L	0	10

=end text

=cut

my $help = 0;
my $man = 0;
my $maf_path = "STDIN";
my $query_path = "";
my $pseudoref_path = "";
my $pseudoref_prefix = "";
my $debug = 0;
my $dispversion = 0;
GetOptions('input_maf|i=s' => \$maf_path, 'query_path|q=s' => \$query_path, 'pseudoref|p=s' => \$pseudoref_path, 'prefix|o=s' => \$pseudoref_prefix, 'debug|d+' => \$debug, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

print STDERR "pseudoref prefix is required\n" if $pseudoref_prefix eq "";
pod2usage(-exitval => 1, -verbose => 2, -output => \*STDERR) if $pseudoref_prefix eq "";

#Open the MAF file, or set it up to be read from STDIN:
print STDERR "Opening MAF file\n" if $debug;
my $maf_fh;
if ($maf_path ne "STDIN") {
   unless(open($maf_fh, "<", $maf_path)) {
      print STDERR "Error opening MAF file: ${maf_path}.\n";
      exit 2;
   }
} else {
   open($maf_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $maf_fh so we can seamlessly handle piping
}

#Extract the positional arguments, which are prefixes for the species:
my @prefixes = @ARGV;

#Make a hash so we can quickly check if the species from an "s" line matches:
print STDERR "Constructing prefix hash\n" if $debug;
my %prefix_hash = map { $_ => undef } @prefixes;

if (scalar(@prefixes) != 2) {
   close($maf_fh);
   print STDERR "It appears you have provided more than 2 prefixes. We only support pairwise alignments for the moment.\n" if scalar(@prefixes) > 2;
   print STDERR "It appears you have provided less than 2 prefixes. We only support pairwise alignments for the moment.\n" if scalar(@prefixes) < 2;
   exit 3;
}

#Give them understandable variable names:
my ($query_spp, $target_spp) = @prefixes[0..1];
print STDERR "Using ${query_spp} as query and ${target_spp} as target.\n" if $debug;

#Read in the query FASTA:
print STDERR "Reading in the query FASTA ${query_path}\n" if $debug;
my %query_seqs = ();
my $query_fh;
unless(open($query_fh, "<", $query_path)) {
   print STDERR "Error opening query FASTA ${query_path}\n";
   exit 4;
}
my $seq_holder = "";
my $header = "";
while (my $line = <$query_fh>) {
   chomp $line;
   if ($line =~ /^>/) { #Header line
      if ($header ne "") { #Not the beginning of the file, so store previous
         print STDERR "Adding scaffold ${header} to ${query_spp} sequences.\n" if $debug > 1;
         $query_seqs{$header} = $seq_holder;
         $seq_holder = "";
      }
      #We expect perfect matches of headers between sequences, including metadata
      $header = substr($line, 1);
   } else {
      $seq_holder .= $line;
   }
}
#Load in the final sequence:
if ($header ne "" and $seq_holder ne "") {
   print STDERR "Adding scaffold ${header} to ${query_spp} sequences.\n" if $debug > 1;
   $query_seqs{$header} = $seq_holder;
}
close($query_fh);

#Read in the pseudoref FASTA:
print STDERR "Reading in the pseudoref FASTA ${pseudoref_path}\n" if $debug;
my %pseudoref_seqs = ();
my $pseudoref_fh;
unless(open($pseudoref_fh, "<", $pseudoref_path)) {
   print STDERR "Error opening pseudoref FASTA ${pseudoref_path}\n";
   exit 5;
}
$header = "";
$seq_holder = "";
while (my $line = <$pseudoref_fh>) {
   chomp $line;
   if ($line =~ /^>/) { #Header line
      if ($header ne "") { #Not the beginning of the file, so store previous
         print STDERR "Adding scaffold ${header} to pseudoref sequences.\n" if $debug > 1;
         $pseudoref_seqs{$header} = $seq_holder;
         $seq_holder = "";
      }
      #We expect perfect matches of headers between sequences, including metadata
      $header = substr($line, 1);
   } else {
      $seq_holder .= $line;
   }
}
#Load in the final sequence:
if ($header ne "" and $seq_holder ne "") {
   print STDERR "Adding scaffold ${header} to pseudoref sequences.\n" if $debug > 1;
   $pseudoref_seqs{$header} = $seq_holder;
}

close($pseudoref_fh);

#Establish file handles for the output INSNPs and BEDs for each prefix:
print STDERR "Opening output INSNP and BED handles\n" if $debug;
my $insnp_fh;
my $bed_fh;
unless(open($insnp_fh, ">", "${pseudoref_prefix}_on_${target_spp}_INSNP.tsv")) {
   print STDERR "Could not open INSNP file for ${query_spp} on ${target_spp} for writing.\n";
   close($maf_fh);
   exit 4;
}
unless(open($bed_fh, ">", "${pseudoref_prefix}_on_${target_spp}_aligned_regions.bed")) {
   print STDERR "Could not open BED file for ${query_spp} on ${target_spp} for writing.\n";
   close($maf_fh);
   exit 4;
}

sub revcomp($) {
   my $input_sequence = shift @_;
   my $reverse_sequence = reverse $input_sequence; #Reverse
   $reverse_sequence =~ tr/AaCcGgTtRrYySsWwKkMmBbDdHhVvNn-/TtGgCcAaYyRrSsWwMmKkVvHhDdBbNn-/; #Complement
   return $reverse_sequence;
}

sub stranded_allele($$) {
   my $input_allele = shift @_;
   my $strand = shift @_;
   return $input_allele unless $strand eq "-";
   return revcomp($input_allele) if $strand eq "-";
}

#Added function for extracting a region from the FASTAs:
sub extractRegion {
   my %seqs = %{shift @_};
   my $scaf = shift @_;
   my $start = shift @_;
   my $len = shift @_;

   if (exists($seqs{$scaf})) {
      return substr($seqs{$scaf}, $start, $len);
   } else {
      return undef;
   }
}

#These two functions were taken from addAlignedGaps.pl:
sub gapPositions($) {
   #Run-length encode the gaps so that inserting terminal gaps is possible
   # with splice()
   my $seq = shift @_;
   my @gap_positions = ([0,0]); #This is actually an array of tuples of start and length of each gap
   my @bases = split //, $seq;
   my $seqlen = scalar(@bases);
   for (my $i = 0; $i < $seqlen; $i++) {
      if ($bases[$i] eq "-") {
         if ($gap_positions[$#gap_positions][0]+$gap_positions[$#gap_positions][1] == $i) {
            #Extend the previous gap by one:
            my @revised_gap = @{$gap_positions[$#gap_positions]};
            $revised_gap[1]++;
            splice(@gap_positions, $#gap_positions, 1, \@revised_gap);
         } else {
            #Create a new gap of length 1:
            push @gap_positions, [$i, 1];
         }
      }
   }
   return \@gap_positions;
}

sub addGaps($$) {
   my $seq = shift @_;
   my @gap_sequence = @{shift @_;};
   my @gapped_sequence = split //, $seq; #We're going to splice in gaps *forwards*
   #Why forwards? Because our tuples are *in* aligned coordinate space, not unaligned coordinate space
   for (my $i = 0; $i <= $#gap_sequence; $i++) {
      print STDERR "splice will give warning for ", $gap_sequence[$i][0], " gap of length ", $gap_sequence[$i][1], " because array is of size ", scalar(@gapped_sequence), "\n" if $gap_sequence[$i][0] > scalar(@gapped_sequence);
      splice(@gapped_sequence, $gap_sequence[$i][0], 0, ('-')x$gap_sequence[$i][1]);
   }
   return join("", @gapped_sequence);
}

#Hash of hashes to store the elements of a MAF alignment:
my %aligned_seqs = ();
#Hash of offsets to stay in the appropriate coordinate space:
my %offsets = ();
#Array of species in the alignment:
my @species_aligned = ();
#Number of p lines we've encountered, which should equal the number of s lines:
my $num_p_lines = 0;

#Keep track of the total number of alignments for debugging:
my $num_a_lines = 0;

print STDERR "Parsing MAF ${maf_path}\n" if $debug;
while (my $line = <$maf_fh>) {
   chomp $line;
   next if $line =~ /^#/ or $line eq ""; #Skip header and empty lines
   my @maf_elems = split /\s+/, $line; #MAF uses padded spacing
   @species_aligned = () if $maf_elems[0] eq "a"; #Reset species on new alignment
   $num_p_lines = 0 if $maf_elems[0] eq "a"; #Reset p count on new alignment
   $num_a_lines++ if $maf_elems[0] eq "a";
   print STDERR "Parsed ${num_a_lines} alignments\n" if $debug and ${num_a_lines} % 1000 == 1 and $maf_elems[0] eq "a";
   if ($maf_elems[0] eq "s") { #For aligned Sequence records, store the details
      my ($species, $scaffold) = split /\./, $maf_elems[1], 2;
      unless (exists($prefix_hash{$species})) {
         print STDERR "Alignment has species not found in provided prefixes: ${line}\n";
         close($maf_fh);
         exit 5;
      }
      push @species_aligned, $species;
      my @seq_arr = split //, uc($maf_elems[6]); #Splitting string into array makes finding differences easier
      $aligned_seqs{$species} = {'scaffold' => $scaffold,
         'start' => $maf_elems[4] eq "-" ? $maf_elems[5]-$maf_elems[2] : $maf_elems[2]+1,
         'size' => $maf_elems[3],
         'strand' => $maf_elems[4],
         'seq' => \@seq_arr};
   }
   if ($maf_elems[0] eq "p") { #Identify SNPs once we're in the Probability lines
      $num_p_lines++;
      next unless $num_p_lines == scalar(@species_aligned); #Skip unless we're on the last p line
      #Make a hash of the aligned species so we can take set differences:
      my %aligned_spp = map { $_ => undef } @species_aligned;
      #Skip the alignment and give a warning if Species2 wasn't found:
      unless (exists($aligned_spp{$target_spp})) {
         print STDERR "Species2 prefix ${target_spp} not found in alignment ${num_a_lines}, skipping.\n";
         next;
      }

      #Sequence length for each "s" record should be identical, so just take
      # the length from the first species in the alignment:
      my $aln_length = scalar(@{$aligned_seqs{$species_aligned[0]}{'seq'}});

      #Initialize the offsets for each coordinate space:
      for my $species (@species_aligned) {
         $offsets{$species} = 0;
      }

      #Construct the BED line for this alignment for the target species:
      my $scaffold = $aligned_seqs{$target_spp}{'scaffold'};
      my $BEDstart = $aligned_seqs{$target_spp}{'strand'} eq "-" ? $aligned_seqs{$target_spp}{'start'}-$aligned_seqs{$target_spp}{'size'} : $aligned_seqs{$target_spp}{'start'}-1;
      my $BEDend = $aligned_seqs{$target_spp}{'strand'} eq "-" ? $aligned_seqs{$target_spp}{'start'} : $aligned_seqs{$target_spp}{'start'}+$aligned_seqs{$target_spp}{'size'}-1;
      print $bed_fh join("\t", $scaffold, $BEDstart, $BEDend), "\n";
      #Determine the interval scaffold and start for the query species:
      my $query_scaffold = $aligned_seqs{$query_spp}{'scaffold'};
      my $query_BEDstart = $aligned_seqs{$query_spp}{'strand'} eq "-" ? $aligned_seqs{$query_spp}{'start'}-$aligned_seqs{$query_spp}{'size'} : $aligned_seqs{$query_spp}{'start'}-1;
      #We assume that the Species2 prefix belongs to your query species
      #Extract the aligned Species2 sequence as a string, not an array:
      my $query_aln = join("", @{$aligned_seqs{$query_spp}{'seq'}});
      #Now identify the positions of gaps in the alignment:
      my @query_gaps = @{gapPositions($query_aln)};
      #Determine the ungapped sequence length by subtracting all gaps from aln_length:
      my $ungapped_length = $aln_length;
      for my $gap_ref (@query_gaps) {
         $ungapped_length -= $gap_ref->[1];
      }
      #Extract the query region and the pseudoref region:
      my $query_seq = extractRegion(\%query_seqs, $query_scaffold, $query_BEDstart, $ungapped_length);
      unless (defined($query_seq)) {
         print STDERR "Something went wrong with extraction, scaffold ${query_scaffold} doesn't exist for ${query_spp}.\n";
         exit 6;
      }
      my $pseudoref_seq = extractRegion(\%pseudoref_seqs, $query_scaffold, $query_BEDstart, $ungapped_length);
      unless (defined($pseudoref_seq)) {
         print STDERR "Something went wrong with extraction, scaffold ${query_scaffold} doesn't exist for the pseudoref.\n";
         exit 7;
      }
      #Revcomp if necessary:
      $query_seq = revcomp($query_seq) if $aligned_seqs{$query_spp}{'strand'} eq '-';
      $pseudoref_seq = revcomp($pseudoref_seq) if $aligned_seqs{$query_spp}{'strand'} eq '-';

      #Add gaps to query and pseudoref, and split into arrays:
      my $gapped_query = addGaps($query_seq, \@query_gaps);
      my @aligned_query = split //, $gapped_query;
      if (uc($query_aln) ne uc($gapped_query)) {
         print STDERR "Something went wrong with the extraction or gap addition, these two don't match...\n";
         close($maf_fh);
         exit 8;
      }
      my $gapped_pseudoref = addGaps($pseudoref_seq, \@query_gaps);
      my @aligned_pseudoref = split //, $gapped_pseudoref;
      #As a test, add gaps into the query sequence and compare to the
      # aligned query sequence:
      #Also include the pseudoref in the comparison
      if ($debug > 1) {
         print STDERR "${scaffold}:${BEDstart} vs. ${query_scaffold}:${query_BEDstart}\n" if $debug > 1;
         print STDERR join("", @{$aligned_seqs{$target_spp}{'seq'}}), "\n" if $debug > 1;
         for (my $i = 0; $i < $aln_length; $i++) {
            print STDERR "|" if uc($aligned_seqs{$query_spp}{'seq'}[$i]) eq uc($aligned_seqs{$target_spp}{'seq'}[$i]) and $debug > 1;
            print STDERR "X" if uc($aligned_seqs{$query_spp}{'seq'}[$i]) ne uc($aligned_seqs{$target_spp}{'seq'}[$i]) and $debug > 1;
         }
         print STDERR "\n" if $debug > 1;
         print STDERR "${query_aln}\n" if $debug > 1;
         for (my $i = 0; $i < $aln_length; $i++) {
            print STDERR "|" if uc($aligned_seqs{$query_spp}{'seq'}[$i]) eq uc($aligned_query[$i]) and $debug > 1;
            print STDERR "X" if uc($aligned_seqs{$query_spp}{'seq'}[$i]) ne uc($aligned_query[$i]) and $debug > 1;
         }
         print STDERR "\n" if $debug > 1;
         print STDERR "${gapped_query}\n" if $debug > 1;
         for (my $i = 0; $i < $aln_length; $i++) {
            print STDERR "|" if uc($aligned_query[$i]) eq uc($aligned_pseudoref[$i]) and $debug > 1;
            print STDERR "X" if uc($aligned_query[$i]) ne uc($aligned_pseudoref[$i]) and $debug > 1;
         }
         print STDERR "\n" if $debug > 1;
         print STDERR "${gapped_pseudoref}\n" if $debug > 1;
      }
      #Now proceed with the normal groundTruthFromMAF procedure in the target_spp coordinate space:
      for (my $i = 0; $i < $aln_length; $i++) {
         #Find SNP, assuming there are 2 spp in the alignment:
         print STDERR join(",", uc($aligned_seqs{$target_spp}{'seq'}[$i]), uc($aligned_pseudoref[$i]), "all"), "\n" if $debug > 2;
         if (uc($aligned_pseudoref[$i]) ne uc($aligned_seqs{$target_spp}{'seq'}[$i]) and $aligned_pseudoref[$i] ne '-' and $aligned_seqs{$target_spp}{'seq'}[$i] ne '-') {
            #Print the SNP to the appropriate INSNP:
            my $insnp_position = $aligned_seqs{$target_spp}{'strand'} eq "-" ? $aligned_seqs{$target_spp}{'start'}-$offsets{$target_spp} : $aligned_seqs{$target_spp}{'start'}+$offsets{$target_spp};
            my $query_allele = stranded_allele($aligned_pseudoref[$i], $aligned_seqs{$target_spp}{'strand'});
            my $target_allele = stranded_allele($aligned_seqs{$target_spp}{'seq'}[$i], $aligned_seqs{$target_spp}{'strand'});
            print STDERR join(",", uc($target_allele), uc($query_allele), "snp"), "\n" if $debug > 2;
            print $insnp_fh join("\t", $aligned_seqs{$target_spp}{'scaffold'}, $insnp_position, $target_allele, $query_allele), "\n";
         }
         for my $species (@species_aligned) {
            $offsets{$species}++ unless $aligned_seqs{$species}{'seq'}[$i] eq "-";
         }
      }
   }
}
close($maf_fh);
print STDERR "Done parsing all ${num_a_lines} alignments from MAF file\n" if $debug;

print STDERR "Closing INSNP and BED files\n" if $debug;

close($insnp_fh);
close($bed_fh);

exit 0;
