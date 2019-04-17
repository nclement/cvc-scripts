#! /usr/bin/perl
# Take two pdbs and put them in a file, then minimize them. Amber strips out
# the chains and messes with the numbers, so this un-splits the minimized PDB
# to match the originals.
#
# For best performance (i.e. it won't break), specify the input proteins in the
# same fashion as you did when you combined them for the minimized bit.

use strict;
use Data::Dumper qw(Dumper);

my $inp_l = shift;
my $inp_r = shift;
my $minned = shift;
my $outp_l = shift;
my $outp_r = shift;

die "usage <orig L> <orig R> <minimized tog> [<output L> <output R>]"
    if !$inp_l || !$inp_r || !$minned;

# Lookup table of acceptable translations.
my %lookup_table;
$lookup_table{"HIS"}{"HIE"} = 1;
$lookup_table{"HIS"}{"HSD"} = 1;

sub sameRes {
  my ($t, $m) = @_;
  #print "Testing [$t] vs [$m]\n";
  return 1 if ($t eq $m);
  return 1 if (defined($lookup_table{$t}{$m}));
  return 0;
}

if (!$outp_l && !$outp_r) {
  ($outp_l = $inp_l) =~ s/(...._l_.\.pdb)/lc($1)/e;
  ($outp_r = $inp_r) =~ s/(...._r_.\.pdb)/lc($1)/e;
}
print "Writing to $outp_l and $outp_r\n";

open (my $MINNED, '<', $minned) or die "Could not open minned filename [$minned]: $!\n";
my @minned = <$MINNED>;
close($MINNED);

sub get_next_res {
  my ($INF, $curres, $curnum) = @_;
  my $newres = 'NA';

  while(<$INF>) {
    next if $_ !~ /^(ATOM|TER)/;
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

my $midx = 0; # Current minned index.

sub process_file {
  my $inf = shift;
  my $outf = shift;

  print "Reading from [$inf]\n";
  my $curres = 'NA';
  my $curnum = '-1';
  my $curchain = '0';

  open (my $INF, "<", $inf) or die "Could not open orig filename [$inf]: $!\n";
  open (my $OUTF, ">", $outf) or die "Could not open output file [$outf]: $!\n";
  while(!eof($INF) and $midx <= $#minned) {
    while ($midx <= $#minned and $minned[$midx] !~ /^(ATOM|TER)/) {
      print $OUTF $minned[$midx];
      $midx++;
    }

    # Get the current residue name and number.
    ($curres, $curnum, $curchain) = &get_next_res($INF, $curres, $curnum);
    print "Found curres $curres:$curchain:$curnum at midx=$midx\n";


    last if ($curres eq 'NA');
    my $minres = substr $minned[$midx], 17, 3;
    my $minnum = substr $minned[$midx], 22, 4;
    if (!&sameRes($curres, $minres)) {
      # Failed!
      print "ERROR: Could not identify correct residue! ",
            "Found min[$minres:$minnum], expected cur[$curres:$curnum]\n";
      exit;
    }
    my $nextminnum = $minnum;
    while (&sameRes($curres, $minres) && $minnum eq $nextminnum) {
      # Also modify the residue number and chain to match.
      substr($minned[$midx], 22, 4) = $curnum;
      substr($minned[$midx], 21, 1) = $curchain;
      # Then print it out.
      #print "Printing [$minres:$minnum]=[$curres:$curnum] -> [$minned[$midx]]\n";
      print $OUTF $minned[$midx++];

      last if $midx == $#minned; # Quit if we've hit the end of the minimized file.
      $minres     = substr $minned[$midx], 17, 3;
      $nextminnum = substr $minned[$midx], 22, 4;
    }
  }
  print "Finished with file\n";
  close($OUTF);
  close($INF);
}

&process_file($inp_l, $outp_l);
&process_file($inp_r, $outp_r);
# Now, process R
