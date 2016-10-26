#! /usr/bin/perl

use strict;

# Gets summary stats for a given column in a file
my $col = shift;
my $fn = @ARGV[0];

my @data = ();
while(<>) {
  my @elems = split;
  push @data, $elems[$col];
}

my $count = @data;
my $avg = &avg(\@data);
my $sd = &sd($avg, \@data);
print "$fn $count $avg $sd\n";

sub avg {
  my $dat = shift;
  my $count = 0;
  my $total = 0;
  foreach (@$dat) {
    $count++;
    $total += $_;
  }

  return 'nan' if $count == 0;
  return $total/$count;
}

sub sd {
  my $avg = shift;
  my $dat = shift;

  return 'nan' if $avg eq 'nan';
  # SD of 0 if only one point
  return 0 if @$dat == 1;

  my $sqtotal = 0;
  foreach(@$dat) {
    $sqtotal += ($avg - $_) ** 2;
  }
  my $std = ($sqtotal / (@$dat - 1)) ** 0.5;
  return $std;
}
