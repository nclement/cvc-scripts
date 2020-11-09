#! /usr/bin/perl
# Fixes the atom count so it starts at 1 for the given file.
# I don't know why, but there are some serious bugs when running the F3Dock
# sampleHinges method without doing this first. Could be a version problem?
#
# Also fixes the residue count, unless a second argument is given.
#
# usage: fix_pdb_atomcount.pl <PDB_INPUT>
#
# Prints to STDOUT

use strict;

my $atom_count = 1;

my $resi_count = 0;
my $count = 0;
my $prev_res = -1;
my $prev_insertion = " ";
my $residue_output = "";

my $INF;
# It's possible that we're just reading from STDIN, so we need to handle this.
if ( -t STDIN ) {
  my $inf = shift;
  open($INF, "<", $inf) or die "Error: Could not open input file [$inf]: $!";
} else {
  $INF = \*STDIN;
}
my $fix_residue = shift;
$fix_residue = $fix_residue ? 0 : 1;


while(<$INF>) {
	if (!/^(ATOM|HETATM)/) {
    #print $_;
    $residue_output .= $_;
    $resi_count++;
		next;
	}
	my $line = $_;
	my $resi = substr($line, 22, 4);
  # If there's something in the insertion column, it means it's been inserted
  # to preserve protein residue numbering.
	my $insertion = substr($line, 26, 1);

	if ($resi != $prev_res || $insertion != $prev_insertion) {
		# Don't print out residues that only have one atom
		if ($resi_count > 1) {
			print $residue_output;
		}
		$residue_output = "";
		$resi_count = 0;
		$count++;
		$prev_res = $resi;
    $prev_insertion = $insertion;
	}
  substr( $line, 5, 6 ) = sprintf("% 6d", $atom_count++);
  if ($fix_residue) {
    substr( $line, 22, 4 ) = sprintf("%4d", $count);
  }
	$residue_output .= $line;
	$resi_count++;
}
if ($resi_count > 1) {
	print $residue_output;
}
