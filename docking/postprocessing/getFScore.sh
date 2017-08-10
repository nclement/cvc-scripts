#! /bin/bash

##########################################################################
# Computing the cRMSD can be misleading, so this script will compute the #
# F-score of residues that are in contact.                               #
#                                                                        #
# The F-score is the harmonic mean between precision and recall, and can #
# be written as:                                                         #
#                                                                        #
#   F-score = 2 * (precision * recall) / (precision + recall)            #
#                                                                        #
# or, subsituting the values of precision and recall:                    #
#                                                                        #
#                 (#TP / #PredP) * (#TP / #ActP)                         #
#   F-score = 2 * ------------------------------                         #
#                 (#TP / #PredP) + (#TP / #ActP)                         #
#                                                                        #
# where:                                                                 #
#   #TP    : number of true positives                                    #
#   #PredP : number of predicted positives                               #
#   #ActP  : number of actual positives                                  #
#                                                                        #
##########################################################################

POST_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=$POST_SCRIPTS_DIR/../../
source $SCRIPTS_DIR/Makefile.def

BOUND_REC=$1
BOUND_LIG=$2
DOCKED_REC=$3
DOCKED_LIG=$4
F2DOCK_OUT=$5

# Comput
