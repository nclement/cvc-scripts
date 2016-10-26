#! /usr/bin/perl

my $count = 0;
my $total = 0;
while(<>) {
	$total += $_;
	$count++;
}
print $total/$count, "\n";
