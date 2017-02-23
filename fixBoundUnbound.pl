#! /usr/bin/perl
####################
# Sometimes, the benchmarks don't have the same residues in both files.
# For this reason, this script will grab both input files and produce two
# output files that have the same number and type of residues. You should
# probably run fixAtomNum.pl after you've finished this script.

my $inpdb1 = shift;
my $inpdb2 = shift;
my $outpdb1 = shift;
my $outpdb2 = shift;

open (my $IN1, "<", $inpdb1) or die "Could not open input 1 [$inpdb1]: $!\n";
open (my $IN2, "<", $inpdb2) or die "Could not open input 2 [$inpdb2]: $!\n";
open (my $OUT1, ">", $outpdb1) or die "Could not open output 1 [$outpdb1]: $!\n";
open (my $OUT2, ">", $outpdb2) or die "Could not open output 2 [$outpdb2]: $!\n";

my %res1 = ();
for (<$IN1>) {
  next if !/^ATOM /;
  my $resn = substr($_, 17, 3);
  my $resi = substr($_, 22, 4);

  # include the residue number so we can sort later.
  $res1{"${resi}_$resn"} .= $_;
}
my %res2 = ();
for (<$IN2>) {
  next if !/^ATOM /;
  my $resn = substr($_, 17, 3);
  my $resi = substr($_, 22, 4);

  # include the residue number so we can sort later.
  $res2{"${resi}_$resn"} .= $_;
}

# Now, create the outputs
for (sort keys %res1) {
  if (defined $res2{$_}) {
    print $OUT1 $res1{$_};
    print $OUT2 $res2{$_};
  }
}
