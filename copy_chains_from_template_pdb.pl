#! /usr/bin/perl
# Some modifications to a PDB will strip out the chains. This will add them
# in and keep the same residue numbering to boot!
#
# Also possible that there are some mutations. The final (optional) arg is a
# comma-separate list of mutations that matches the following regex:
#   \w\d+\w{3}-\w{3}(:?,\w\d\w{3]-\w{3})
# where each pattern is CHAIN RESNUM RESID - RESID, e.g. for the following two
# changes:
#  - MET at position 147 of chain C changed to THR
#  - ASP at position 148 of chain C changed to ILE
# use the following:
#  C147MET-THR,C148ASP-ILE
#
# For best performance (i.e. it won't break), the residues in both files muste
# be in the same order.

use strict;
use Data::Dumper qw(Dumper);

my $tmpl = shift;
my $inp = shift;
my $outp = shift; # Can be same as `inp` but must be different from `tmpl`.
my $muts = shift; 

die "usage <orig> <input> <output> [<mutation_str>]"
if !$tmpl || !$inp || !$outp;

# Lookup table of acceptable translations.
my %lookup_table;
$lookup_table{"HIS"}{"HIE"} = 1;
$lookup_table{"HIS"}{"HSD"} = 1;
$lookup_table{"ASP"}{"ASN"} = 1;

# Lookup table of accetped mutations.
my %mutations;
my @mut_strs = split ',', $muts;
for (@mut_strs) {
  if (!/(\w)(\d+)(\w{3})-(\w{3})/) {
    print "ERROR: mutation doesn't match the regex: $_\n";
  }
  $mutations{"$1$2$3"} = $4;
}
print "Possible mutations:\n";
for my $k (keys %mutations) {
  print "  $k -> $mutations{$k}\n";
}


sub sameRes {
  # <TEMPL_CHAIN> <TEMPL_RESNUM> <TEMPL_RES> <TEST_RES>
  my ($c, $x, $t, $m) = @_;
  #print "Testing [$t] vs [$m]\n";
  return 1 if ($t eq $m);
  return 1 if (defined($lookup_table{$t}{$m}));

  # Check and see if it's an allowed mutation.
  # Spaces in TEMPL_RESNUM make this break.
  (my $xx = $x) =~ s/^\s*(\d+)\s*$/$1/;
  #print "Checking mutations{$c|$xx|$t} = $m\n";
  return 1 if (defined($mutations{"$c$xx$t"}) and $mutations{"$c$xx$t"} eq $m);

  return 0;
}

print "Reading input from [$inp]\n";
open (my $INP, '<', $inp) or die "Could not open input filename [$inp]: $!\n";
my @input = <$INP>;
close($INP);

sub get_next_res {
  my ($INF, $curres, $curnum) = @_;
  my $newres = 'NA';

  while(<$INF>) {
    next if $_ !~ /^(ATOM)/;
    my $thisres = substr($_, 17, 3);
    my $thisnum = substr($_, 22, 4);
    my $thischain = substr($_, 21, 1);
    #print "Looking at [$thisres] $_\n";
    if ($thisres ne $curres or $thisnum ne $curnum) {
      return ($thisres, $thisnum, $thischain);
    }
  }

  return ($newres, '-1', '0');
}

my $midx = 0; # Current input index.

print "Reading template from [$tmpl]\n";
my $curres = 'NA';
my $curnum = '-1';
my $curchain = '0';

open (my $TMPL, "<", $tmpl) or die "Could not open orig filename [$tmpl]: $!\n";
open (my $OUTF, ">", $outp) or die "Could not open output file [$outp]: $!\n";
while(!eof($TMPL) and $midx <= $#input) {
  while ($midx <= $#input and $input[$midx] !~ /^(ATOM)/) {
    print $OUTF $input[$midx];
    $midx++;
  }

  # Get the current residue name and number.
  ($curres, $curnum, $curchain) = &get_next_res($TMPL, $curres, $curnum);
  #print "Found curres $curres:$curchain:$curnum at midx=$midx\n";


  last if ($curres eq 'NA');
  my $inpres = substr $input[$midx], 17, 3;
  my $inpnum = substr $input[$midx], 22, 4;
  if (!&sameRes($curchain, $curnum, $curres, $inpres)) {
    # Failed!
    print "ERROR: Could not identify correct residue! ",
    "Found inp[X:$inpres:$inpnum], expected cur[$curchain:$curres:$curnum]\n";
    exit -1;
  }
  my $nextinpnum = $inpnum;
  while (&sameRes($curchain, $curnum, $curres, $inpres) && $inpnum eq $nextinpnum) {
    # Also modify the residue number and chain to match.
    substr($input[$midx], 22, 4) = $curnum;
    substr($input[$midx], 21, 1) = $curchain;
    # Then print it out.
    #print "Printing [$inpres:$inpnum]=[$curres:$curnum] -> [$input[$midx]]\n";
    print $OUTF $input[$midx++];

    last if $midx == $#input; # Quit if we've hit the end of the input file.
    $inpres     = substr $input[$midx], 17, 3;
    $nextinpnum = substr $input[$midx], 22, 4;
  }
}
print "Finished with file\n";
close($OUTF);
close($TMPL);
