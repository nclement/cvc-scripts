#! /bin/bash

input=$1
output=$2

grep -e "ATOM.*N   " -e "ATOM.*CA  " -e "ATOM.*C   " -e "ATOM.*O   " $input > $output
