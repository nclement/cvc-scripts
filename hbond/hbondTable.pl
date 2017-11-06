#! /usr/bin/perl
#########################
# Will create a suitable R script from the hbond output file.
#
# Input files from hbondInterface.sh output
# Writes to stdout
########################

use strict;
use Math::Trig;

my $name = shift;  # Give it some name
my $hbond = shift; # .hbond
my $pdb = shift;   # .pdb
my $mol2 = shift;  # .mol2
my $no_header = shift; # Some non-zero value here if the header should not be printed out.

die "usage: hbondTable.pl <name> <*.hbond> <*.pdb> <*.mol2>\n" if (!$hbond || !$pdb || !$mol2) ;


# Determines if this is a backbone or not.
my %backboneAtoms = (
  " C  " => 1,
  " O  " => 1,
  " N  " => 1,
  " CA " => 1 );


my %bonds = &readMol2($mol2);

open(my $PDBIN, "<", $pdb) or die "Error: Could not open pdb input file [$pdb]: $!\n";
my %aStat = ();
while(my $line = <$PDBIN>) {
  next if (!($line =~ /^ATOM/));
  my $atomi = int(substr( $line, 6, 5 ));
  my $aname = substr( $line, 12, 4);
  my $aa = substr( $line, 17, 3 );
  my $chain = substr( $line, 21, 1);
  my $resi = substr( $line, 22, 4 );
  my $x = 1.0 * substr($line, 30, 8);
  my $y = 1.0 * substr($line, 38, 8);
  my $z = 1.0 * substr($line, 46, 8);

  $aStat{$atomi} = {
    aname=>$aname,
    aa=>$aa,
    chain=>$chain,
    resi=>$resi,
    x=>$x,
    y=>$y,
    z=>$z
  };
}

# Print out the header
if (!$no_header) {
  print "name ";
  print "Didx Hidx Aidx ABidx AB2idx ";
  print "DaType HaType AaType ABaType AB2aType ";
  print "Dresi Hresi Aresi ABresi AB2resi ";
  print "Dresn Hresn Aresn ABresn AB2resn ";
  print "Dchain Hchain Achain ABchain AB2chain ";
  print "Dbackbone Abackbone ";
  print "AHbondCount AneighborCount ";
  print "bondType bondType2 ";
  print "hbe_type hbw_type ";
  print "energy distance theta psi psi2 interface\n";
}

open(my $HBIN, "<", $hbond) or die "Error: Could not open hbond input file [$hbond]: $!\n";

# Read the header.
my $hline = <$HBIN>;
$hline = <$HBIN> while($hline =~ /^#/);
$hline =~ /(\d+) (\d+)/;
my $atoms_a = $1;
my $atoms_b = $2;

# Array of hydrogen bonds.
my @HBonds = ();
# Counts of numbers of times each atom is involved.
my %HbDonor = ();
my %HbAcctr = ();

DONOR_LIST:
while(<$HBIN>) {
  # Ignore comments.
  next if /^#/;
  # D_ATOMI D_RESI D_TYPE A_ATOMI A_RESI A_TYPE HBE_TYPE HB_WEIGHT_TYPE ENERGY
  /(\d+) +\d+ +\S+ +(\d+) +\d+ +\S+ +(\d+) +(\d+) +(\S+)/;
  my $donor = $1;
  my $acctr = $2;  # acceptor
  my $hbe_type = $3; # hbond energy type
  my $hbw_type = $4; # hbond weight type
  my $engy = $5;   # energy

  my $from = $donor <= $atoms_a ? "A" : "B";
  my $to = $acctr <= $atoms_a ? "A" : "B";
  my $interface = 0;
  if ($from ne $to) {
    $interface = 1;
  }
  if (!defined($aStat{$donor})) {
    print STDERR "Error: donor $donor not defined!!\n";
  } elsif (!defined($aStat{$acctr})) {
    print STDERR "Error: acceptor $acctr not defined!!\n";
  } else {
    # Count how many times this atom is involved in bonding.
    $HbDonor{$donor}++;
    $HbAcctr{$acctr}++;
    # Donor will always be hydrogen.
    if (!defined($bonds{$donor})) {
      print STDERR "Error: hydrogen $donor($aStat{$donor}{aname}) not bonded to anything!!\n";
      next DONOR_LIST;
    }
    if (!defined($bonds{$acctr})) {
      print STDERR "Error: acceptor $acctr($aStat{$acctr}{aname}) not bonded to anything!!\n";
      next DONOR_LIST;
    }

    if (keys $bonds{$acctr} > 2) {
      print STDERR "WARN: bond $acctr ($aStat{$acctr}{aa}:$aStat{$acctr}{aname}) has ", scalar keys $bonds{$acctr}, " bonds!";
      print STDERR " just using the first two.\n";
    }
    if (keys $bonds{$acctr} < 1) {
      print STDERR "WARN: bond $acctr ($aStat{$acctr}{aa}:$aStat{$acctr}{aname}) has ", scalar keys $bonds{$acctr}, " bonds!";
      print STDERR " must have at least one.\n";
    }
    # Get our atom numbers.
    my $D = (keys $bonds{$donor})[0];
    my $AB = (keys $bonds{$acctr})[0];
    my $AB2 = -1;
    if (keys $bonds{$acctr} > 1) {
      $AB2 = (keys $bonds{$acctr})[1];
    }
    # Ignore all hydrogen friends.
    if ($AB2 != -1 and $aStat{$AB2}{aname} =~ /.*H.*/) {
      $AB2 = -1;
    }
    if ($aStat{$AB}{aname} =~ /.*H.*/) {
      $AB = $AB2; # This had better exist!!
      $AB2 = -1;
      # Get the next available for AB2
      my $i = 2;
      while ($i < keys $bonds{$acctr} && $aStat{(keys $bonds{$acctr})[$i]}{anam} =~ /.*H.*/) {
        $AB2 = (keys $bonds{$acctr})[$i];
      }
    }
    if ($AB == -1) {
      print STDERR "Error: donor ($aStat{$acctr}{aname}) only bound to hydrogens ", keys $bonds{$acctr}, "\n";
      next DONOR_LIST;
    }
    my $bond = $bonds{$acctr}{$AB};
    my $bond2 = $AB2 == -1 ? "NA" : $bonds{$acctr}{$AB2};

    push @HBonds, [ $D, $donor, $acctr, $AB, $AB2, $engy, $interface, $bond, $bond2, $hbe_type, $hbw_type ];

  }
}

for (@HBonds) {
  (my $D, my $H, my $acctr, my $AB, my $AB2, my $energy, my $interface, my $bond, my $bond2, my $hbe_type, my $hbw_type) = @{$_};
    # Now, print the statistics.
    print "$name $D $H $acctr $AB $AB2 ";
    # Atom name
    print "$aStat{$D}{aname} $aStat{$H}{aname} $aStat{$acctr}{aname} $aStat{$AB}{aname} ";
    $AB2 == -1 ? print "NA " : print "$aStat{$AB2}{aname} ";
    # Residue number
    print "$aStat{$D}{resi} $aStat{$H}{resi} $aStat{$acctr}{resi} $aStat{$AB}{resi} ";
    $AB2 == -1 ? print "NA " : print "$aStat{$AB2}{resi} ";
    # Residue name
    print "$aStat{$D}{aa} $aStat{$H}{aa} $aStat{$acctr}{aa} $aStat{$AB}{aa} ";
    $AB2 == -1 ? print "NA " : print "$aStat{$AB2}{aa} ";
    # Chain
    print "$aStat{$D}{chain} $aStat{$H}{chain} $aStat{$acctr}{chain} $aStat{$AB}{chain} ";
    $AB2 == -1 ? print "NA " : print "$aStat{$AB2}{chain} ";
    # Is Backbone
    print &isBackbone($aStat{$D}{aname}), " ", &isBackbone($aStat{$acctr}{aname}), " ";
    # Number of bonds for acceptor.
    print "$HbAcctr{$acctr} ", scalar keys $bonds{$acctr}, " ";
    print "$bond $bond2 ";
    print "$hbe_type $hbw_type ";
    #      energy         distance                                        theta                  psi
    printf "%.5f %.5f %.5f %.5f ",  $energy , &dist($H, $acctr), &angle($D, $H, $acctr), &angle($H, $acctr, $AB);
    #printf "%.5f %.5f %.5f %.5f ",  $energy , &dist($H, $acctr), &angle($acctr, $H, $D), &angle($AB, $acctr, $H);
    if ($AB2 != -1) {
      printf "%.5f ", &angle($H, $acctr, $AB2);
      #printf "%.5f ", &angle($AB2, $acctr, $H);
    } else {
      print "NA ";
    }
    print "$interface\n";
}

sub isBackbone {
  my $atom = shift;
  return defined $backboneAtoms{$atom} ? 1 : 0;
}

# Computes the distance between two atoms specified by their index.
sub dist {
  my $donor = shift;
  my $acceptor = shift;
  return sqrt(($aStat{$donor}{x} - $aStat{$acceptor}{x})**2 + 
              ($aStat{$donor}{y} - $aStat{$acceptor}{y})**2 + 
              ($aStat{$donor}{z} - $aStat{$acceptor}{z})**2);
}

# Computes the angle a-b-c
sub angle {
  my $a = shift;
  my $b = shift;
  my $c = shift;

  my @AB = &subAtm($a, $b);
  my @BC = &subAtm($c, $b);

  return &angle3d(\@AB, \@BC);
}

sub subAtm {
  my $a = shift;
  my $b = shift;

  return ($aStat{$a}{x} - $aStat{$b}{x},
          $aStat{$a}{y} - $aStat{$b}{y},
          $aStat{$a}{z} - $aStat{$b}{z});
}

# Returns the angle ABC, defined by vectors AB and BC
sub angle3d {
  my $AB = shift;
  my $BC = shift;

  # Normalize
  my $lAB = sqrt(@{$AB}[0]**2 + @{$AB}[1]**2 + @{$AB}[2]**2);
  my $lBC = sqrt(@{$BC}[0]**2 + @{$BC}[1]**2 + @{$BC}[2]**2);

  # Get dot product.
  my $dotp = @{$AB}[0]*@{$BC}[0] + @{$AB}[1]*@{$BC}[1] + @{$AB}[2]*@{$BC}[2];
  my $deg = rad2deg(&acos($dotp / ($lAB * $lBC)));
#   print STDERR "(a=$d1, b=$a, c=$d2)\n";
#   print STDERR "a <- c($aStat{$d1}{x}, $aStat{$d1}{y}, $aStat{$d1}{z})\n";
#   print STDERR "b <- c($aStat{$a}{x}, $aStat{$a}{y}, $aStat{$a}{z})\n";
#   print STDERR "c <- c($aStat{$d2}{x}, $aStat{$d2}{y}, $aStat{$d2}{z})\n";
#   print STDERR "acos((a-b) %*% (b-c) / (sum((a-b)^2) * sum((b-c)^2))) * 180 / pi\n";
#   print STDERR "  I get $deg\n";
  return $deg;
}

sub acos { atan2( sqrt(1 - $_[0] * $_[0]), $_[0] ) }

sub dotProduct {
  my $AB = shift;
  my $BC = shift;
}
# Cross product between two vectors of size 3
sub crossProd3d {
  my $a = shift;
  my $b = shift;
  return (@{$a}[1] * @{$b}[2] - @{$a}[2] * @{$b}[1],
          @{$a}[2] * @{$b}[0] - @{$a}[0] * @{$b}[2],
          @{$a}[0] * @{$b}[1] - @{$a}[1] * @{$b}[0]);
}

# Computes the torsion angle between x-A-B-y, along the bond A-B
sub torsion {
  my $x = shift;
  my $A = shift;
  my $B = shift;
  my $y = shift;

  # Get the vectors from atoms.
  my @vprev = &subAtm($x, $A);
  my @vbond = &subAtm($A, $B);
  my @vnext = &subAtm($B, $y);

  # Get the normalized planes between them.
  my @plane1n = &crossProd3d(\@vbond, \@vprev);
  my @plane2n = &crossProd3d(\@vnext, \@vbond);
  &normalize(@plane1n);
  &normalize(@plane2n);

  # Get the angle between the two planes.
  my $torsion = &acos(&dotProd(\@plane1n, \@plane2n));
  
  # See if we should flip the angle.
  my @temp_v1 = &crossProd3d(\@plane1n, \@vbond);
  if (&dotProd(\@temp_v1, \@plane2n) < 0) {
    $torsion *= -1;
  }

  return $torsion;
}

sub normalize {
  my $vec = shift;

  print STDERR "@{$vec}\n";
  my $n = sqrt(@{$vec}[0]**2 + @{$vec}[1]**2 + @{$vec}[2]**2);
  @{$vec}[0] /= $n;
  @{$vec}[1] /= $n;
  @{$vec}[2] /= $n;
}

sub readMol2 {
  my $mol2file = shift;

  open(my $M2, "<", $mol2file) or die "Error: Could not open mol2 file [$mol2file]: $!\n";

  # Read until you find @<TRIPOS>BOND
  while(<$M2>) { last if /@<TRIPOS>BOND/; }

  my %bonds = ();
  # Does this end??
  while(<$M2>) {
    # Break if we find something empty or another TRIPOS
    last if /@<TRIPOS>/;
    last if /^$/;

    # The format of each line is:
    #   BOND_ID START_ATOM END_ATOM BOND_TYPE
    # where BOND_TYPE is one of:
    #  - 1 = single 
    #  - 2 = double 
    #  - 3 = triple 
    #  - am = amide 
    #  - ar = aromatic 
    #  - du = dummy 
    #  - un = unknown (cannot be determined from the parameter tables) 
    #  - nc = not connected
    /(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/;
    my $from = $2;
    my $to = $3;
    my $type = $4;

    # Set a bond both ways.
    $bonds{$from}{$to} = $type;
    $bonds{$to}{$from} = $type;
  }

  return %bonds;
}
