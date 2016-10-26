#! /bin/bash

FROM=$1
PATTERN=$2
export TO=$3

find $FROM -name "$PATTERN" -print0 | xargs -t -P 8 -0 -I file mv file $TO
