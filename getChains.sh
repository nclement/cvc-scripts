#! /bin/bash

pdb=$1

echo -n $(grep "^ATOM" $pdb | cut -c 22 | sort | uniq) | sed 's/ //g'
