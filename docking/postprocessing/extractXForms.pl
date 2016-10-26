#! /usr/bin/perl
# Will extract a number of xforms from an F2Dock output file.
# Can then be used in the conformationGenerator program.

my $numConf = shift;

die "usage extractXForms.pl <# XForms> <F2Dock.txt>\n" if (@ARGV < 1);
my @lines = <>;

print "$numConf\n";
my $count = 0;
for (reverse @lines) {
  next if /^#/;
  last if ++$count == $numConf;

  my @elems = split;
  print join(" ",@elems[15 .. 26]), "\n";
}
