#! /usr/bin/perl -w
##########################################################
# pdbchm.pl
# structure charmm minimization on input pdb.
##########################################################
use strict;
use File::Basename;
use Getopt::Long;

my $home = $ENV{'HOME'};

my $wdir="";
my $charmm="/org/centers/cvc/software/share/usr.linux.x86_64/c36b2/exec/gnu/charmm";
my $ertf="pdbamino.rtf";
my $eprm="parm.prm";
my $nsteps=200;
my $minimize=1;
my $eloops="";
my $edisu="";
my $smod="";
my $clean=1;
my $osuffix="_cmin";

GetOptions
(
  "loops=s" => \$eloops,
  "disu=s" => \$edisu,
  "smod=s" => \$smod,
  "rtf=s" => \$ertf,
  "prm=s" => \$eprm,
  "wdir=s" => \$wdir,
  "charmm=s" => \$charmm,
  "nsteps=i" => \$nsteps,
  "dont-clean" => sub { $clean = 0 },
  "dont-minimize" => sub { $minimize = 0 },
  "osuffix:s" => \$osuffix,
);

my $debug=0;

if (scalar @ARGV < 2)
{
  print_short_args ();
  exit;
}
my $pdb=shift @ARGV;
$pdb=lc($pdb);
my @chains = @ARGV;

die "RTF not found\n" unless (-e $ertf);
die "PRM not found\n" unless (-e $eprm);

# check name
my ($pdbname, $pdbdir, $pdbsuffix) = fileparse ($pdb, qr/\.[^.]*/);
if ($pdbsuffix eq ".pdb")
{
  $pdb=substr ($pdb, 0, length($pdb)-4);
}

# deal vith ? chain id
my ($cho, $lin);
my $ncho=-1;
my @nk;
my $mi;
foreach $cho (@chains)
{
  $ncho=$ncho+1;
  if($cho eq "?")
  {
    @nk=(@nk,$ncho);
  }
}

if(@nk)
{
  for($mi=@nk-1; $mi>=0; $mi--)
  {
    splice(@chains, $nk[$mi],1);
  }

  my @filss = `ls $wdir$pdb-?.*.pdb`;
  chomp @filss;
  @filss = sort { substr($a, -9) cmp substr($b, -9) } @filss;
  foreach $lin(@filss)
  {
    my $nu=length($lin)-10;
    my $no=length($wdir.$pdb."-");
    $lin=substr ($lin, $no, $nu-$no+1);
  }
  my %chainh = map {$_, 1} @chains; # get unique elements
  for (@filss)
  {
    next if $chainh{$_}++;
    push @chains, $_;
  }
}

my $nloo=0;
my @loo;

if($eloops)
{
  @loo=split(/,/, $eloops);
  $nloo=@loo;
  if($nloo-int($nloo/3)*3!=0)
  {
    print "inconsistent definition of loops\n";
    print "$eloops\n";
    exit
  }
  $nloo=$nloo/3;
}

my $i;
my @se;
my @fr;
my @sr;
for($i=0; $i<$nloo; $i++)
{
  @se=(@se,($smod.(shift @loo)));
  @fr=(@fr,shift @loo);
  @sr=(@sr,shift @loo);
}

my $ndisu=0;
my @disul;

if($edisu)
{
  @disul=split(/,/, $edisu);
  $ndisu=@disul;
  if($ndisu-int($ndisu/4)*4!=0)
  {
    print "inconsistent definition of disulfide pairs\n";
    print "$eloops\n";
    exit
  }
  $ndisu=$ndisu/4;
}

my @sds1;
my @sds2;
my @rds1;
my @rds2;
for($i=0; $i<$ndisu; $i++)
{
  @sds1=(@sds1,shift @disul);
  @rds1=(@rds1,shift @disul);
  @sds2=(@sds2,shift @disul);
  @rds2=(@rds2,shift @disul);
}

my $rtf=$ertf;
my $par=$eprm;
my $oscript="$wdir$$.inp";

my $fil;
my $cha;

open OSCRIPT, "> $oscript" or die "couldn't open file: $oscript\n";
print OSCRIPT "* read protein\n";
print OSCRIPT "* multiple chains with breaks\n";
print OSCRIPT "*\n";
print OSCRIPT "\n";

print OSCRIPT "bomlev -2\n";
print OSCRIPT "! RTF AND PARAM FILEiS\n";
print OSCRIPT "open read card unit 2 name $rtf\n";
print OSCRIPT "read rtf card unit 2\n";
print OSCRIPT "close unit 2\n";
print OSCRIPT "open read card unit 2 name $par\n";
print OSCRIPT "read param card unit 2\n";
print OSCRIPT "close unit 2\n";
print OSCRIPT "\n";

print OSCRIPT "! SETUP UNIV COORD READER TO READ PDB FORMAT, PUT B FACTOR IN WMAIN\n";
print OSCRIPT "read univ\n";
print OSCRIPT "* modified pdb setup\n";
print OSCRIPT "*\n";
print OSCRIPT "UNKN\n";
print OSCRIPT "SEGID 1 3\n";
print OSCRIPT "RESID 24 4\n";
print OSCRIPT "TYPE 13 4\n";
print OSCRIPT "RESN 18 3\n";
print OSCRIPT "x 31 8\n";
print OSCRIPT "y 39 8\n";
print OSCRIPT "z 47 8\n";
print OSCRIPT "w 61 6\n";
print OSCRIPT "end\n";
print OSCRIPT "\n";

print OSCRIPT "! GENERATE CHAINS\n";
foreach $cho (@chains){
  $cho=lc($cho);
  my @files = `ls $wdir$pdb-$cho.*.pdb`;
  chomp @files;
  my $nfil=0;
  $cha=$smod.$cho;
  foreach $fil (@files){
    $nfil=$nfil+1;
    print OSCRIPT "open unit 10 read card name $fil\n";
    print OSCRIPT "read sequ pdb unit 10\n";
    my $seg;
    if($nfil>1){$seg=$cha."2";}
    else{$seg=$cha;}
    print OSCRIPT "gener $seg setup first none  last none warn\n";
    if($nfil>1){
      print OSCRIPT "join $cha $seg\n";
    }
    print OSCRIPT "close unit 10\n";
    print OSCRIPT "\n";
  }
}

print OSCRIPT "! build disulfide bonds\n";
for($i=0; $i<$ndisu; $i++)
{
  print OSCRIPT "patch disu $sds1[$i] $rds1[$i] $sds2[$i] $rds2[$i] setup sort warn\n";
}
print OSCRIPT "\n";

my $count = 0;

my $ngap =0;
my @sega;
my @fga;
my @sga;
my $fter;

foreach $cho (@chains){
  $cho=lc($cho);
  my @files = `ls $wdir$pdb-$cho.*.pdb`;
  chomp @files;
  my $lg="     ";
  $cha=$smod.$cho;
  foreach $fil (@files){
    open PDB, "< $fil" or die "couldn't open file $fil: $!\n";
    my @pdblines = <PDB>;
    chomp @pdblines;
    close PDB;
    my $line=shift @pdblines;
    my $sn= substr ($line, 0, 3);
    my $fres = substr ($line, 22, 4);
    my $fresi= substr ($line, 22, 5);
    $line=pop @pdblines;
    my $lres = $fres;
    my $lresi= $fresi;
    if(@pdblines > 0)
    {
      $line=pop @pdblines;
      $lres = substr ($line, 22, 4);
      $lresi= substr ($line, 22, 5);
    }
    if($lg ne "     ")
    {
      $fga[$ngap]=$lg;
      $sga[$ngap]=$fresi;
      $sega[$ngap]=$cha;
      $ngap=$ngap+1;
    }
    else
    {
      $fter=$fresi;
    }
    $lg=$lresi;

    $count = $count+1;
    my $offs = $count - $fres;
    my $tmp = $fres + $offs;
    $count = $lres + $offs;
    print OSCRIPT "rename segid $sn select segid $cha end\n";
    print OSCRIPT "open unit 10 read card name $fil\n";
    print OSCRIPT "read coor univ resi unit 10 select ires $tmp:$count end\n";
    print OSCRIPT "close unit 10\n";
    print OSCRIPT "\n";
    print OSCRIPT "rename segid $cha select segid $sn end\n";
  }
  print OSCRIPT " define some select segid $cha .and. resid $fter .and. Type P end\n";
  print OSCRIPT "if ?NSEL .GT. 0 patch 5TER $cha $fter setup warn\n";

}


print OSCRIPT "AUTO ANGL DIHE !update angles, dihedrals\n";
print OSCRIPT "ic purge   ! CLEANUP IC TABLE\n";
print OSCRIPT "ic param   ! GET MISSING BONDS AND ANGLES FROM PARAMETER FILE\n";
print OSCRIPT "ic build   ! PLACE ANY MISSING COORDS, E.G. TERMINAL O ON CO2-\n";
print OSCRIPT "\n";

print OSCRIPT "! CHECK FOR MISSING HEAVY ATOM COORDS\n";
print OSCRIPT "define testa sele ( .not. type H* ) .and. ( .not. init ) show end\n";
print OSCRIPT "\n";

print OSCRIPT "! PATCHES\n";
print OSCRIPT "\n";

for($i=0; $i<$ngap; $i++)
{
  print OSCRIPT "define some1 sele segid $sega[$i] .and. resi $fga[$i] .and. type C show end\n";
  print OSCRIPT "set reon ?SELRESN\n";
  print OSCRIPT "define some2 sele segid $sega[$i] .and. resi $sga[$i] .and. type N show end\n";
  print OSCRIPT "set retw ?SELRESN\n";
  print OSCRIPT "\n";

  print OSCRIPT "coor dist cut 1.4 sele some1 end sele some2 end\n";
  print OSCRIPT "IF ?NPAIR EQ 0 GOTO OUT$i\n";
  print OSCRIPT "\n";

  print OSCRIPT "if reon EQ GLY if retw eq GLY GOTO JOGG$i\n";
  print OSCRIPT "if reon EQ GLY if retw eq PRO GOTO JOGP$i\n";
  print OSCRIPT "if reon EQ GLY GOTO JOGA$i\n";
  print OSCRIPT "if retw EQ PRO GOTO JOAP$i\n";
  print OSCRIPT "if retw EQ GLY GOTO JOAG$i\n";
  print OSCRIPT "\n";

  print OSCRIPT "patch JOAA  $sega[$i]  $fga[$i] $sega[$i] $sga[$i] setup warn\n";
  print OSCRIPT "GOTO OUT$i\n";
  print OSCRIPT "\n";

  print OSCRIPT "label JOAG$i\n";
  print OSCRIPT "patch JOAG  $sega[$i]  $fga[$i] $sega[$i] $sga[$i] setup warn\n";
  print OSCRIPT "GOTO OUT$i\n";
  print OSCRIPT "\n";

  print OSCRIPT "label JOAP$i\n";
  print OSCRIPT "patch JOAP  $sega[$i]  $fga[$i] $sega[$i] $sga[$i] setup warn\n";
  print OSCRIPT "GOTO OUT$i\n";
  print OSCRIPT "\n";

  print OSCRIPT "label JOGA$i\n";
  print OSCRIPT "patch JOGA  $sega[$i]  $fga[$i] $sega[$i] $sga[$i] setup warn\n";
  print OSCRIPT "GOTO OUT$i\n";
  print OSCRIPT "\n";

  print OSCRIPT "label JOGP$i\n";
  print OSCRIPT "patch JOGP  $sega[$i]  $fga[$i] $sega[$i] $sga[$i] setup warn\n";
  print OSCRIPT "GOTO OUT$i\n";
  print OSCRIPT "\n";

  print OSCRIPT "label JOGG$i\n";
  print OSCRIPT "patch JOGG  $sega[$i]  $fga[$i] $sega[$i] $sga[$i] setup warn\n";
  print OSCRIPT "\n";

  print OSCRIPT "label OUT$i\n";
}

print OSCRIPT "! SAVE THE COMPLETED PSF\n";
print OSCRIPT "open unit 1 write card name $wdir$pdb$osuffix.psf\n";
print OSCRIPT "write psf card oldpsf unit 1\n";
print OSCRIPT "* the subunits of $pdb; initial CHARMM psf\n";
print OSCRIPT "*\n";
print OSCRIPT "\n";

print OSCRIPT "! INDICATE LOOPS TO MINIMIZE\n";
for($i=0; $i<$nloo; $i++)
{
  print OSCRIPT "define testd sele (segid $se[$i] .and. resi $fr[$i]:$sr[$i]) show end\n";
  print OSCRIPT "define teste sele testa .or. testd end\n";
  print OSCRIPT "define testa sele teste end\n";
}
print OSCRIPT "\n";

print OSCRIPT "! SAVE THE IC TABLE FILLED WITH XTAL DATA\n";
print OSCRIPT "ic fill\n";
print OSCRIPT "open unit 1 write card name $wdir$pdb\_$$.ic\n";
print OSCRIPT "write ic card unit 1\n";
print OSCRIPT "* the subunits of $pdb\n";
print OSCRIPT "*\n";
print OSCRIPT "\n";

if ($minimize) {
  print OSCRIPT "! USE HBUILD TO REBUILD H ATOMS; SPINS METHYLS, ETC. TO LOCAL MINIMUM\n";
  print OSCRIPT "coor init sele type H* .and. .not. resn HG end\n";
  print OSCRIPT "hbuild sele type H* .and. .not. resn HG end\n";
  print OSCRIPT "\n";
}

print OSCRIPT "! CHECK FOR ANY MISSING COORDS\n";
print OSCRIPT "define testb sele .not. init show end\n";
print OSCRIPT "\n";

print OSCRIPT "! CHECK FOR EMPTINESS\n";
print OSCRIPT "define tot select all end\n";
print OSCRIPT "if ?nsel .eq. 0 stop\n";

print OSCRIPT "define testc sele testa .or. (TYPE H*) show end\n";
print OSCRIPT "CONS HARM SELE resn HG .or. .not. testc  END - \n";
print OSCRIPT "     FORCE 20 MASS \n";
print OSCRIPT "\n";

if ($minimize) {
  #print OSCRIPT "ENER ATOM ACE IEPS 4.0 SEPS 78.0 ALPHa 1.2 SIGMa 3 SWITch - \n";
  print OSCRIPT "ENER ATOM  IEPS 4.0 SEPS 78.0 ALPHa 1.2 SIGMa 3 SWITch - \n";
  print OSCRIPT "       VDIS VSWI CUTNB 15.0 CTONNB 8.0 CTOFNB 13.0\n";
  print OSCRIPT "\n";

  print OSCRIPT "MINIMIZE ABNR NSTEP $nsteps NPRINT 10 \n";
  print OSCRIPT "\n";
}

print OSCRIPT "CONS FIX SELE all end\n";
print OSCRIPT "\n";

#print OSCRIPT "ENER ATOM ACE IEPS 4.0 SEPS 78.0 ALPHa 1.2 SIGMa 3 SWITch - \n";
print OSCRIPT "ENER ATOM IEPS 4.0 SEPS 78.0 ALPHa 1.2 SIGMa 3 SWITch - \n";
print OSCRIPT "      VDIS VSWI CUTNB 15.0 CTONNB 8.0 CTOFNB 13.0\n";
print OSCRIPT "\n";

print OSCRIPT "OPEN UNIT 52 NAME  $wdir$pdb\_$$.ener WRITE FORM\n";
print OSCRIPT "WRITE ENER UNIT 52\n";
print OSCRIPT "CLOSE UNIT 52\n";
print OSCRIPT "\n";


print OSCRIPT "! SAVE THE COMPLETED CHARMM VERSION OF THE PDB XTAL COORDS\n";
print OSCRIPT "open unit 1 write card name $wdir$pdb\_$$.crd\n";
print OSCRIPT "write coor card unit 1\n";
print OSCRIPT "* subunits of $pdb; initial CHARMM coords\n";
print OSCRIPT "*\n";
print OSCRIPT "\n";

print OSCRIPT "open unit 1 write card name $wdir$pdb$osuffix.pdb\n";
print OSCRIPT "title unit 1\n";
print OSCRIPT "* this file generated using CHARMM\n";
print OSCRIPT "write coor pdb offi unit 1\n"; # Write official PDB file
print OSCRIPT "\n";

print OSCRIPT "stop\n";
close OSCRIPT;
my $out = $oscript."\.out";
system "$charmm < $oscript>$out";
if ($clean)
{
  system "rm $oscript";
  system "rm $out";
  system "rm -f $wdir$pdb\_$$.crd";
  system "rm -f $wdir$pdb\_$$.ener";
  system "rm -f $wdir$pdb\_$$.ic";
}



sub print_short_args
{
  printf ("usage: %s options PDB_code ch1 ch2 ...\n", $0);
  printf ("options: -loops seg,first,last,...\n");
  printf ("         -disu seg1,res1,seg2,res2,...\n");
  printf ("         -smod segment_prefix\n");
  printf ("         --wdir wdir\n");
  printf ("         --charmm charmm\n");
  printf ("         --prm prm\n");
  printf ("         --rtf rtf\n");
  printf ("         --nsteps nsteps\n");
  printf ("         --dont-clean\n");
  printf ("         --dont-minimize\n");
  printf ("         --osuffix=s (default: _cmin)\n");
}
