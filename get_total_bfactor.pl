#! /usr/bin/perl
# ###############################
# Calculates the total B-factor for the given PDB file
# #############################

my $pdb = shift;

open(my $PDB, "<", $pdb) or die "Could not open PDB file [$pdb]: $!\n";

my $count = 0;
my $total = 0;
while(<$PDB>) {
	next if !/^ATOM/;

	my $b = substr($_, 60, 6);
	#print "$b\n";
	$total += $b;
	$count++;
}

printf "%.2f\t%d", $total, $count;
