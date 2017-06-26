#! /usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Fatal qw(open close);

my $home   = $ENV{'HOME'};
my $debug  = 0;
my $helpf  = 0;              # help flag
my $renumf = 0;              # renumber flag
my $segidf = 0;              # add segid flag
my $hem    = 0;

GetOptions(
    "h"     => \$helpf,
    "help"  => \$helpf,
    "renum" => \$renumf,
    "segid" => \$segidf,
    "hem"   => \$hem,
);

if ($helpf) {
    print_args();
    exit 1;
}

if ( scalar @ARGV < 1 ) {
    print_short_args();
    print_help_args();
    exit 1;
}

my $pdb = shift @ARGV;

# create names
my ( $pdbname, $pdbdir, $pdbsuffix ) = fileparse( $pdb, qr/\.[^.]*/ );
my $pdbpre = $pdbdir . $pdbname;
$pdbpre=lc($pdbpre);
$pdbsuffix=lc($pdbsuffix);

# read pdb
open my $PDB, '<', $pdb;
my @pdblines = <$PDB>;
chomp @pdblines;
close $PDB;

my $segid       = 0;
my $ofile       = "NULL";
my $chain       = "x";
my $myresi      = 1;
my $prevaltloc  = " ";       # alternate location indicator
my $prevresname = "NULL";
my $prevresi    = -inf;
my $previcode   = "NULL";    # insertion code

for my $line (@pdblines) {
    my $record = substr( $line, 0, 6 );
    if ($debug) {
        print "*** $ofile ***\n";
        print "\n";
        print "$line\n";
    }

    if ( ( $record eq "TER   " || $record eq "END   " ) && $myresi != 1 ) {

        #print "newseg from ter,end: myresi: $myresi\n";
        newseg();
        next;
    }
    if ( $record eq "ATOM  " || $record eq "HETATM" ) {
        my $altloc  = substr( $line, 16, 1 );
        my $resname = substr( $line, 17, 3 );
        $chain = lc( substr( $line, 21, 1 ) );
        my $resi  = substr( $line, 22, 4 );
        my $icode = substr( $line, 26, 1 );

        if ( $resname eq "HOH" || $resname eq "WAT" ) {
            next;
        }
        if ( $resname eq "MSE" )    # modified MET
        {
            $record = "ATOM  ";
        }
        if ( $resname eq "ASH" ) {
            $resname = "ASPH";
            substr( $line, 17, 4 ) = "ASPH";
        }
        if ( $resname eq "GLH" ) {
            $resname = "GLUH";
            substr( $line, 17, 4 ) = "GLUH";
        }
        if ( $resname eq "CYX" ) {
            $resname = "CYS";
            substr( $line, 17, 3 ) = "CYS";
        }
        if ( $resname eq "HIP" ) {
            $resname = "HSC";
            substr( $line, 17, 3 ) = "HSC";
        }
        if ( $resname eq "HID" ) {
            $resname = "HIS";
            substr( $line, 17, 3 ) = "HIS";
        }
        if ( $resname eq "HIE" ) {
            $resname = "HSD";
            substr( $line, 17, 3 ) = "HSD";
        }
        if ( $resname eq "NGL" )    #really NGLY
        {
            $resname = "GLYP";
            substr( $line, 17, 4 ) = "GLYP";
        }
        if ( $resname eq "NPR" )    #really NPRO
        {
            $resname = "PROP";
            substr( $line, 17, 4 ) = "PROP";
        }

        $chain     = "x"    if ( $chain eq " " );
        $prevresi  = $resi  if ( $prevresi == -inf );
        $previcode = $icode if ( $previcode eq "NULL" );

        my $residelta = $resi - $prevresi;
        if ( $residelta > 0 || $previcode ne $icode ) {
            $myresi++;
        }

        #if (abs ($residelta) > 1)
        if ( $residelta < 0 || $residelta > 1 ) {

            #print "newseg from resdelta: $residelta\n";
            newseg();
        }
        if ( $residelta == 0 && $previcode ne $icode ) {

            #print "newseg from icode: $resi$icode\n";
            newseg();
        }
        if ( $record eq "HETATM" && $prevresname ne $resname ) {

            #print "newseg from HETATM new resname: $resname\n";
            newseg();
        }

        # remember values
        $prevresname = $resname;
        $prevresi    = $resi;
        $previcode   = $icode;

        my $hetresname = $resname;
        $hetresname =~ s/\s+//;
        $hetresname = 'x' if ( $hetresname eq "" );
        $hetresname = lc $hetresname;
        $chain      = "h$hetresname" if ( $record eq "HETATM" and not ($hetresname eq "hem" and $hem));
        if ( $ofile eq "NULL" ) {
            ofile();
            system "rm -f $ofile";    # rm ofile for appending
        }
        ofile();

        # modify line
        my $lline = $line;            # local line
        substr( $lline, 0, 6 ) = "ATOM  ";
        if ($renumf) {
            substr( $lline, 22, 4 ) = sprintf( "%4d", $myresi );
        }
        if ($segidf) {
            my $segidi = 72;          # segid index
            if ( length $lline <= $segidi ) {
                $lline = sprintf( "%-80s", $lline );
            }
            substr( $lline, $segidi, 4 ) = sprintf( "%04d", $segid );
        }
        if ( substr( $lline, 12, 4 ) eq "OCT1" ) {
            substr( $lline, 12, 4 ) = " O  ";
        }

        open my $OFILE, '>>', $ofile;
        print $OFILE "$lline\n";
        close $OFILE;
    }
}

# add END to the last file
newseg();

sub newseg {
    if ( $ofile eq "NULL" ) {
        return;
    }

    ### print final line
    open my $OFILE, '>>', $ofile;
    print $OFILE "END\n";
    close $OFILE;

    $segid++;

    # reset vals
    $ofile       = "NULL";
    $myresi      = 1;
    $prevaltloc  = " ";
    $prevresname = "NULL";
    $prevresi    = -inf;
    $previcode   = "NULL";
}

sub ofile {
    $ofile = sprintf "$pdbpre-$chain.%04d$pdbsuffix", $segid;
}

############################
sub print_help_args {
    print STDERR "try '$0 -h' for a list of options\n";
}

sub print_short_args {
    print STDERR "usage: $0 [options] PDB\n";
    print STDERR "process the pdb file\n";
}

sub print_args {
    print_short_args();

    print STDERR "\n";
    print STDERR "Options:\n";
    print STDERR "\n";

    print STDERR "-h, --help    Display this help and exit\n";
    print STDERR "--renum       renumber the residues in the output files\n";
    print STDERR "--segid       add segment ids to the output files\n";
    print STDERR "--hem         keep file name based on chain for the heteroatom residue hem\n";
}
