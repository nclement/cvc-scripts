#! /usr/bin/perl
#########################
# Will create a suitable R script from the hbond output file.
#
# Input .hbond file from hbondInterface.sh output
# Writes to stdout
########################

use strict;

my $hbond = shift;
# Get the name of the input
my $fname = "hbond_res.pdf";
if ($hbond =~ /(.*).hbond/) {
  $fname = "$1_hbond.pdf";
}

open(my $HBIN, "<", $hbond) or die "Error: Could not open hbond input file [$hbond]: $!\n";

# Read the header.
my $hline = <$HBIN>;
$hline =~ /(\d+) (\d+)/;
my $atoms_a = $1;
my $atoms_b = $2;

print("require(ggplot2)\n");
print("df <- data.frame(donor=numeric(), acceptor=numeric(), energy=numeric(), from=character(), to=character(), interface=logical())\n");
while(<$HBIN>) {
  /(\d+) (\d+) (\S+)/;
  my $donor = $1;
  my $acctr = $2; # acceptor
  my $engy = $3;   # energy

  my $from = $donor <= $atoms_a ? "A" : "B";
  my $to = $acctr <= $atoms_a ? "A" : "B";
  # Need to handle the first thing
  print("df <- rbind(df, data.frame(donor=$donor, acceptor=$acctr, energy=$engy, from='$from', to='$to', ");
  if ($from ne $to) {
    print("interface=T");
  } else {
    print("interface=F");
  }
  print("))\n");
}

print("q <- ggplot(df, aes(donor, acceptor, color=energy)) + \n");
print("     geom_point()\n");
print("q <- q + theme(axis.text.x=element_blank(), axis.text.y=element_blank())\n");
print("q <- q + geom_hline(yintercept=$atoms_a, col='black', lty=2) + \n");
print("         geom_vline(xintercept=$atoms_a, col='black', lty=2)\n");
print("ggsave('$fname', q, height=7, width=8)\n");
