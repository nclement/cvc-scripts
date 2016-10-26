#! /usr/bin/perl

my $fn = @ARGV[1];

my @data = <>;
my $avg = &avg(\@data);
my $sd = &sd($avg, \@data);
print "$fn $avg $sd\n";

sub avg {
  my $dat = shift;
  my $count = 0;
  my $total = 0;
  foreach (@$dat) {
    $cout++;
    $total += $_;
  }

  return $total/$count;
}

sub sd {
  my $avg = shift;
  my $dat = shift;

  my $sqtotal = 0;
  foreach(@$data) {
    $sqtotal += ($avg - $_) ** 2;
  }
  my $std = ($sqtotal / (@$data - 1)) ** 0.5;
  return $std;
}
