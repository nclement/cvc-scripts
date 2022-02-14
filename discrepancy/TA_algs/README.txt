This package contains files to compile the Threshold Accepting-based
star discrepancy lower bounds algorithms described in (M. Gnewuch,
M. Wahlström, C. Winzen: A New Randomized Algorithm to Approximate the
Star Discrepancy Based on Threshold Accepting).  Things here are
pretty basic, interface-wise, but at least on a unix-like system the
following should work. (I have not tried, but I expect that either a
modern Mac (Mac OS X) or a Windows machine running Cygwin would work
the same way.)

NB: This package has been modified slightly by Nathan Clement
(nclement@utexas.edu) to work with modern C++ compilers. You may need to edit
the two variables `GXX` and `GPP` to work with your compiler, e.g. consider
setting them to

    GXX=gcc
    GPP=g++


COMPILING

For the exact versions used in the paper, try

    make TA_basic

which creates an executable TA_basic, and

    make TA_improved

which creates two executables TA_improved_delta and
TA_improved_bardelta (because delta and bardelta are optimized
separately in the TA_improved algorithm).  If you want to experiment
with other combinations of settings, all the code is available in the
package. 


RUNNING

All programs read a point set either from a given filename or via
stdin.  The normal usage would be something like the following.

    ./TA_improved_bardelta -iter 100000 -trials 10 (dim) (points) (filename)
    ./TA_improved_delta -iter 100000 -trials 10 (dim) (points) (filename)

or

    generate_pointset | ./TA_improved_bardelta -iter 100000 -trials 10 (dim) (points)
    generate_pointset | ./TA_improved_delta -iter 100000 -trials 10 (dim) (points)

if you have a program that generates a point set on stdout. (Again,
delta and bardelta are computed by different programs for TA_improved,
and will have to be run separately like this.)


The programs currently have a kind of "debug flag" active, meaning
they generate a lot of output (every improving set found).  If you
want to cut this down, look for lines

    #define PRINT_ALL_UPDATES

in the source files and comment them out, i.e., replace by

    // #define PRINT_ALL_UPDATES

That should be all. Hope you finds this useful.

Greetings,
Magnus Wahlström
