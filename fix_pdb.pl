#! /usr/bin/perl

use strict;

# Must provide the chain
my $fakechain = @ARGV > 1 ? shift : "A";

my $float = qr/-?\d+\.\d+/;
my $regex = qr/^(ATOM|HETATM) #ATOM
		\s*(\d+) #POSITION
		\s+(\S+) #ATOM NAME
		\s+([A-Z]+) #RES NAME
		\s+([A-Z])?	#CHAN ID (optional)
		\s+(\d+)	#RES NUMBER
		\s*($float) # X CORD
		\s*($float) # Y CORD
		\s*($float) # Z CORD
		\s*($float) # OCCUPANCY
		\s*($float) # TEMP FACTOR
		(.*)$/x;
while(<>) {
	if (!/$regex/) {
		print "$_";
		next;
	}
	my @matches = ($_ =~ /$regex/);
	$matches[4] |= $fakechain;

	printf "%-6s%5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n", @matches;
}
