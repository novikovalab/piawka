#!/bin/sh
"exec" "gawk" -f "$0" "--" "$@" && 0 {}
#
# This script summarizes piawka results counted over several loci.
# It guesses if VERBOSE=1 was set by the number of fields.
# If not, it takes mean value weighted by nUsed across all lines with same locus, pop1, pop2 and metric.
# If yes, it divides sum of numerators by sum of denominators.

BEGIN{ OFS="\t" }

NR==1 { 
  if ( NF < 8 ) {
    print "piawka output (11 columns) is required!" > "/dev/stderr"
    exit 1
  } 
}
{ 
  nSites[$1,$3,$4,$6]+=$2
  allnUsed[$1,$3,$4,$6]+=$5
  # For theta Watterson, weigh denominator by nUsed = numerator = #segr sites
  # (need a smarter way to average harmonic numbers?)
  if ( $6=="thetaW" ) {
    nUsed[$1,$3,$4,$6]=$5" "nUsed[$1,$3,$4,$6]
    denominator[$1,$3,$4,$6]=$9" "denominator[$1,$3,$4,$6]
  } else {
    denominator[$1,$3,$4,$6]+=$9
  }
  numerator[$1,$3,$4,$6]+=$8
  nGeno[$1,$3,$4,$6]+=$10
  nMiss[$1,$3,$4,$6]+=$11
}

END{
  printf "\033[2K\r" > "/dev/stderr" # empty stderr line
  for (i in nSites) {
    split(i, pops, SUBSEP)
    if ( pops[4]=="thetaW" ) {
      split(nUsed[i], weights)
      split(denominator[i], values)
      for (j in weights) { finvalue[i] += ( values[j] * weights[j] ) }
      finvalue[i] = numerator[i]*allnUsed[i]/finvalue[i]"\t"numerator[i]"\t"finvalue[i]/allnUsed[i]"\t"nGeno[i]"\t"nMiss[i]
    } else {
      finvalue[i]=numerator[i]/denominator[i]"\t"numerator[i]"\t"denominator[i]"\t"nGeno[i]"\t"nMiss[i]
    }
    print pops[1], nSites[i], pops[2], pops[3], allnUsed[i], pops[4], finvalue[i]
  }
}

