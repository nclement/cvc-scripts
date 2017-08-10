#! /usr/bin/perl

# This script will print a table of the pairwise RMSD of all proteins in the
# given file. Should be easy to slurp into R and do 'image' on the resulting
# matrix.
#
# Can specify the interface atoms if it is desired to obtain the RMSD only
# along these atoms.
#
# If PDB_fn has more than one column, only the first column will be used.

use strict;
use Parallel::ForkManager;
use Data::Dumper;

my $PDB_fn = shift;
my $interface_atoms = shift; # Optional

my $max_proc = 32;
my $rmsd_prog = "/work/01872/nclement/scripts/getcRMSD_pymol_atomlist.sh";

die "usage: pairwise_rmsd <PDB_FN>" if !$PDB_fn;

# Default is all atoms.
$interface_atoms = '*' if !$interface_atoms;

open (my $PDBS, "<", $PDB_fn) or die "Could not open pdb file: $!\n";

my @pdbs = ();
my @rmsds = ();
my @crmsds = ();

while(<$PDBS>) {
  chomp;
  /^(\S+)/;
  push @pdbs, $1;
}

my $max_so_far = 0;
# Start the parallel manager
my $pm = Parallel::ForkManager->new($max_proc);
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result_ref) = @_;
    #print STDERR Dumper($result_ref), "\n";
    $rmsds[$result_ref->{i}][$result_ref->{j}] = $result_ref->{rmsd};
    $crmsds[$result_ref->{i}][$result_ref->{j}] = $result_ref->{crmsd};
    if ($result_ref->{j} > $max_so_far) {
      $max_so_far = $result_ref->{j};
      &printProgress($result_ref->{i}, $max_so_far, scalar(@pdbs));
    }
  });

# Print blank lines for starting
print STDERR "Progress:\n\n";
printProgress(0, 0, scalar(@pdbs));
for (my $i = 0; $i < @pdbs; $i++) {
  $max_so_far = 0;
  # Fancy stuff for pretty printing
#   for (my $j = 0; $j < @pdbs + 10; $j++) {
#     print STDERR " ";
#   }
#  print STDERR "\r$i/", scalar(@pdbs);

  push @rmsds, [];
  push @crmsds, [];

  # Add the RMSD from each of the previous proteins.
  for (my $j = 0; $j < $i; $j++) {
    push @{$rmsds[$i]}, $rmsds[$j][$i];
    push @{$crmsds[$i]}, $crmsds[$j][$i];
  }
  # Add zero for yourself.
  push @{$rmsds[$i]}, 0;
  push @{$crmsds[$i]}, 0;
  INNER_J:
  for (my $j = $i + 1; $j < @pdbs; $j++) {
    #print "$i $j\n";
    # Add a hold value.
    push @{$rmsds[$i]}, 0;
    push @{$crmsds[$i]}, 0;
    
    # Start the thread.
    $pm->start($i*scalar(@pdbs) + $j) and next INNER_J;
    #$pm->start($j) and next INNER_J;
    #$pm->start and next INNER_J;
    my $rmsd = `$rmsd_prog $pdbs[$i] $pdbs[$j] '*'`;
    my $crmsd = `$rmsd_prog $pdbs[$i] $pdbs[$j] $interface_atoms`;
    chomp $rmsd;
    chomp $crmsd;
    #print STDERR "$rmsd_prog $pdbs[$i] $pdbs[$j] $interface_atoms\n";
    #my %this_rmsd = (i=>$i, j=>$j, rmsd=>$rmsd);
    #print "Here, thing is:", Dumper(%this_rmsd), "\n";
    #$pm->finish(0, \%this_rmsd);
    $pm->finish(0, {i=>$i, j=>$j, rmsd=>$rmsd, crmsd=>$crmsd});
  }
  $pm->wait_all_children;
}
print STDERR "\n";

for (my $i = 0; $i < @pdbs; $i++) {
  print "RMSD $pdbs[$i] ";
  for (my $j = 0; $j < @pdbs; $j++) {
    print "$rmsds[$i][$j] ";
  }
  print "\n";
}
for (my $i = 0; $i < @pdbs; $i++) {
  print "cRMSD $pdbs[$i] ";
  for (my $j = 0; $j < @pdbs; $j++) {
    print "$crmsds[$i][$j] ";
  }
  print "\n";
}

sub printProgress {
  my $i = shift;
  my $j = shift;
  my $total = shift;

  # Total window is 80, and we have:
  #  - 10 for beginning
  #  - 4 for end
  # Should give us 66-2 for middle
  my $width = 79 - 16;

  my $progress_i = 1.0 * $i / $total;
  my $progress_j = 1.0 * $j / $total;
  printf STDERR "\e[1A\r%4d/%4d [%*s]%3d%%\n", $i, $total, -$width, "#" x int($progress_i*$width+0.5), $progress_i*100;
  printf STDERR "\r%4d/%4d [%*s]%3d%%", $j, $total, -$width, "#" x int($progress_j*$width+0.5), $progress_j*100;
}
