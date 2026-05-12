#!/bin/bash

SRR=$1
echo -e "$SRR\t"`curl -s "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=$SRR&accept-alternate-locations=yes" | tr ',' '\n'| perl -ne 'print "$1\n" if (m/\"(https.*SRR\d+)\"/)'`
