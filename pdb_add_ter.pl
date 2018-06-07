#! /usr/bin/perl

use strict;

my @lines = <>;

if (grep /^TER/, @lines) {
  print @lines;
  exit;
}
# Otherwise, add the TER at the end.
# ATOM   2164  OE2 GLU R1628     -87.350  37.795   9.083  1.00 40.00           O  
# becomes
# TER    2165      GLU R1628
my $last_str = $lines[$#lines];
my $atomi = int(substr $last_str, 6, 11) + 1;
my $res = substr $last_str, 17, 3;
my $chain = substr $last_str, 21, 1;
my $resi = substr $last_str, 22, 4;
print @lines;
printf "TER   %5d      %s %s%s\n", $atomi, $res, $chain, $resi;
