#! /usr/bin/perl

my $pdb = shift;

open(my $PDBIN, "<", $pdb) or die "Could not open $pdb: $!\n";

my @PDB = <$PDBIN>;
close($PDBIN);

# Open it up again as an output stream
open(my $PQROUT, ">", $pdb);
my $atom
for (@PDB) {

}
