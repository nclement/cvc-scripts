#!/usr/bin/perl
# Script for converting from HingeProt fcc.txt file to the kind of file we'd
# like for the dumb sampleHinges script.

my $fcc_txt = shift;
my $levels = shift; # How many levels should we do?
my $num_samples = shift; # How many samples should this file have?

open(my $FCC, "<", $fcc_txt) or die "Could not open fcc.txt file: $!";

my %hinges = ();

while(<$FCC>) {
  /(\d+) [nc] (\d+) (\S+)/;
  next if $1 > $levels;

  # Add up all the different movements.
  $hinges{$2} += $3;
}

my $num_hinges = keys %hinges;
print <<EOM;
1
$num_samples
-5 5
$num_hinges
EOM

for my $k (keys %hinges) {
  print "$k $k\n$k $hinges{$k}\n";
}
