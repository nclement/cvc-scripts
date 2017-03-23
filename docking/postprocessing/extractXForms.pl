#! /usr/bin/perl
# Will extract a number of xforms from an F2Dock output file.
# Can then be used in the conformationGenerator program.

my $numConf = shift;

die "usage extractXForms.pl <# XForms> <F2Dock.txt> (redirects to stdout)\n" if (@ARGV < 1);
my @lines = <>;

my $format = 

print "$numConf\n";
my $count = 0;
my $length = 0;
for (reverse @lines) {
  if (/OUTPUT FORMAT: (\d+)/) {
    $length
  }
  next if /^#/;
  last if ++$count == $numConf;

  my @elems = split;
  # Old version of F2Dock
  if (@elems == 54) {
    print join(" ",@elems[40 .. 51]), "\n";
  } else {
    print join(" ",@elems[15 .. 26]), "\n";
  }
}
