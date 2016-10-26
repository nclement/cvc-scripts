#! /usr/bin/perl

use strict;

my $infile = shift;
open INF, $infile or die "ERROR: Could not open input file [$infile]: $!\n";

$infile =~ /^(.*)\.pdb/;
my $base = $1;

my $prev_domain = "";

while(<INF>) {
  if (/^TER\s+/) {
    print OUTF $_;
    #print $_;
  }
  next if !/^ATOM/;

  /^.{21}(.)/;
  if ($prev_domain ne $1) {
    $prev_domain = $1;
    close(OUTF);
    print "$base-$1.pdb\n";
    open OUTF, ">$base-$1.pdb";
  }
  print OUTF $_;
  
}
