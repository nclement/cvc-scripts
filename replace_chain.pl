#! /usr/bin/perl
# Replaces the chain of the PDB with the one supplied.

use strict;

my $chain = shift;
# Second arg is PDB

while(<>) {
	if (!/^(ATOM|HETATM)/) {
		print $_;
		next;
	}
	my $line = $_;
  substr($line, 21, 1) = $chain;
  print $line;
}
