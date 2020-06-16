#! /usr/bin/perl

use strict;

my %lookup_table = (
  "ALA"=>"A",
  "ARG"=>"R",
  "ASN"=>"N",
  "ASP"=>"D",
  "ASX"=>"B",
  "CYS"=>"C",
  "GLU"=>"E",
  "GLN"=>"Q",
  "GLX"=>"Z",
  "GLY"=>"G",
  "HIS"=>"H",
  "ILE"=>"I",
  "LEU"=>"L",
  "LYS"=>"K",
  "MET"=>"M",
  "PHE"=>"F",
  "PRO"=>"P",
  "SER"=>"S",
  "THR"=>"T",
  "TRP"=>"W",
  "TYR"=>"Y",
  "VAL"=>"V",
);
sub get_residue {
  my $r = shift;
  return $lookup_table{$r} if defined $lookup_table{$r};
  return "_XXX $r XXX_";
}

my $resi_ = -1;
my $icode_ = "";
my $chain_ = "";

while (<>) {
  next if (!/^ATOM/);
  my $resi = substr $_, 22, 4;
  my $resn = substr $_, 17, 3;
  my $icode = substr $_, 26, 1;
  my $chain = substr $_, 21, 1;

  print "/" if ($chain_ ne "" and $chain ne $chain_);
  $chain_ = $chain;

  # From the standard:
  #   "The combination of residue numbering and insertion code defines the unique
  #   residue"
  if ("$resi$icode" ne "$resi_$icode_") {
    print "$chain $resi ", get_residue($resn), "\n";
  }
  $resi_ = $resi;
  $icode_ = $icode;
}
print "\n";
