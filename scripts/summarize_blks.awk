#!/bin/sh
"exec" "gawk" -f "$0" "--" "$@" && 0 {}
#
# This script summarizes piawka results counted over several loci.
# It simply divides sum of numerators by sum of denominators.

BEGIN{ OFS="\t" }

NR==1 { 
  if ( NF < 13 ) {
    print "piawka output (13 columns) is required!" > "/dev/stderr"
    exit 1
  } 
}
$9 != "value" && $9==$9 { # exclude header and NaNs 
  idx=$1 SUBSEP $4 SUBSEP $5 SUBSEP $6 SUBSEP $8 SUBSEP
  if (!seen[idx]) {
    seen[idx]=1
    start[idx]=$2
    end[idx]=$3
  }
  if ($2 < start[idx]) { start[idx]=$2 }
  if ($3 > end[idx]) { end[idx]=$2 }
  allnUsed[idx]+=$7
  denominator[idx]+=$11
  numerator[idx]+=$10
  nGeno[idx]+=$12
  nMiss[idx]+=$13
}

END{
  printf "\033[2K\r" > "/dev/stderr" # empty stderr line
  for (i in seen) {
    split(i, pops, SUBSEP)
    finvalue[i]=numerator[i]/denominator[i]"\t"numerator[i]"\t"denominator[i]"\t"nGeno[i]"\t"nMiss[i]
    print pops[1], start[i], end[i], pops[2], pops[3], pops[4], allnUsed[i], pops[5], finvalue[i]
  }
}

