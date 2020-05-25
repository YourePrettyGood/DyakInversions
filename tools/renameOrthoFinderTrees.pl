#!/usr/bin/perl
use POSIX;
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $SCRIPTNAME = "renameOrthoFinderTrees.pl";
my $VERSION = "1.0";
#Version 1.0 derived from remapOrthogroups.pl
#Applies to OrthoFinder Newick trees, where the transcript names are generic
# and need to be mapped back to their original names

=pod

=head1 NAME

renameOrthoFinderTrees.pl - Converts IDs in OrthoFinder Newick trees to original annotation IDs

=head1 SYNOPSIS

renameOrthoFinderTrees.pl [options]

 Options:
  --help,-h,-?         Display this help documentation
  --prot_id_map,-i     Input concatenated map of OrthoFinder protein IDs to
                       original annotation/transcript IDs
                       (column 1 is OrthoFinder ID, column 2 is transcript ID)
  --tree,-t            Path to Newick tree produced by OrthoFinder
  --debug,-d           Output extra information to STDERR

=head1 DESCRIPTION

renameOrthoFinderTrees.pl converts the protein IDs in any supplied
Newick tree from OrthoFinder back into transcript IDs from the
original annotations using a concatenated map produced by combining
the .map files made by prep_proteomes.mk's proteome target.

=cut


my $prot_id_map = "";
my $treefile = "";
my $help = 0;
my $man = 0;
my $debug = 0;
my $dispversion = 0;

GetOptions('prot_id_map|i=s' => \$prot_id_map, 'tree|t=s' => \$treefile, 'debug|d+' => \$debug, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

my $mapfh;
if (! -e ${prot_id_map}) {
   print STDERR "Unable to open protein ID map file ${prot_id_map} -- does not exist.\n";
   exit 1;
}
open($mapfh, "<", $prot_id_map) or die "Unable to open protein ID map file ${prot_id_map}";
my %prot_map = ();
while (my $line = <$mapfh>) {
   chomp $line;
   my @mapelems = split /\t/, $line;
   #Map is from OrthoFinder input ID (key) to original annotation/transcriptome ID (value):
   $prot_map{${mapelems[1]}} = ${mapelems[0]};
   print STDERR "Added ", join(" => ", ${mapelems[1]}, ${mapelems[0]}), " to OrthoFinderID=>protID map\n" if $debug > 1;
}
close($mapfh);

sub renameTree {
   my $tree = shift @_;
   my %prot_map = %{shift @_};
   my $debug = shift @_;
   
   if ($tree =~ /,/) { #Recurse, and merge
      my @branches = split /,/, $tree;
      
   } else { #Base case
      
   }
   return(join(",", @branches));
}

my $treefh;
if (! -e ${treefile}) {
   print STDERR "Unable to open OrthoFinder Newick tree ${treefile} -- does not exist.\n";
   exit 2;
}
open($treefh, "<", $treefile) or die "Unable to open OrthoFinder Newick tree ${treefile}";

#Parse the tree, performing substitutions as needed:
my $tree = <$treefh>;
chomp $tree;
my $renamed_tree = renameTree($tree, \%prot_map, $debug);
print $renamed_tree, ";\n";
close(treefh);
