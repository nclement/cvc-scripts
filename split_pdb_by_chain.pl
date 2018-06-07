#! /usr/bin/perl

###################################################################
## Will split a PDB by chains. Kind of crude, so see if it works? #
###################################################################

use strict;

my $pdb = $ARGV[0];
die "Usage: split_pdb_by_chain.pl <pdb>\n" if (not defined $pdb);
$pdb =~ /^(.*)\.pdb/;
my $pdb_base = $1;

my $chain = "";
my $pdb_str = "";
my @valid_lines = qw(ATOM HETATM);
LINE: while(<>) {
  if (not (/^ATOM/ or /^HETATM/)) {
    if (/^TER/ or /^END/) {
      $pdb_str .= $_;
    }
    next;
  }

  my $tch = substr $_, 21, 1;
  $chain = $tch if $chain eq "";
  if ($tch ne $chain) {
    my $fn = $pdb_base . "_$chain.pdb";
    print "Writing to $fn\n";
    open(my $OUTF, ">", $fn);
    print $OUTF $pdb_str;
    close $OUTF;
    $pdb_str = "";
    $chain = $tch;
  }
  $pdb_str .= $_;
}

my $fn = $pdb_base . "_$chain.pdb";
print "Writing to $fn\n";
open(my $OUTF, ">", $fn);
print $OUTF $pdb_str;
close $OUTF;
$pdb_str = "";
