#! /usr/bin/perl
##
# Some of the PDBs that are produced don't have residue identifyers, and that
# causes problems in Amber. This takes an original PDB and a simulated PDB and
# makes sure the first 23 characters (which includes the residue numbers) are
# identical.
##

use strict;

my $origPDB = shift;
my $simPDB = shift;
my $samechars = 23;

open my $OP, "<", $origPDB or die "Could not open original PDB: $!\n";
open my $SP, "<", $simPDB or die "Could not open simulated PDB: $!\n";

my @OLines = <$OP>;
my @SLines = <$SP>;
close ($OP);
close ($SP);


my $l_lines = 0;
my $o_lines;
for (0..$#OLines) {
	$o_lines = 0.0+$_;
	next if $OLines[$o_lines] !~ /^ATOM/;
	my @orig = split //, $OLines[$o_lines];
	die "Error: Lines are not equal" if $l_lines > $#SLines;
	my @sim = split //, $SLines[$l_lines];

	splice @sim, 21, 1, $orig[21];
	@SLines[$l_lines++] = join '', @sim;
}

open $SP, ">", $simPDB;
for (@SLines) {
	print $SP $_;
}
close $SP;
