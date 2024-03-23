#!/usr/bin/mawk -f
#
# This script summarizes piawka results counted over several loci.
# It guesses if VERBOSE=1 was set by the number of fields.
# If not, it takes mean value weighted by nUsed across all lines with same locus, pop1, pop2 and metric.
# If yes, it divides sum of numerators by sum of denominators.

BEGIN{ OFS="\t" }

{ 
  nSites[$1,$3,$4,$6]+=$2
  allnUsed[$1,$3,$4,$6]+=$5
  if ( NF>7 ) {
      numerator[$1,$3,$4,$6]+=$8
      denominator[$1,$3,$4,$6]+=$9
      nGeno[$1,$3,$4,$6]+=$10
      nMiss[$1,$3,$4,$6]+=$11
  } else { 
      nUsed[$1,$3,$4,$6]=$5" "nUsed[$1,$3,$4,$6]
      value[$1,$3,$4,$6]=$7" "value[$1,$3,$4,$6]
  }
}

END{
  for (i in nSites) {
    split(i, pops, SUBSEP)
    if ( denominator[i] ) {
        finvalue[i]=numerator[i]/denominator[i]"\t"numerator[i]"\t"denominator[i]"\t"nGeno[i]"\t"nMiss[i]
    } else {
        split(nUsed[i], weights, " ")
        split(value[i], values, " ")
        for (j in weights) { finvalue[i] += ( values[j] * weights[j] ) / allnUsed[i] }
    }
    print pops[1], nSites[i], pops[2], pops[3], allnUsed[i], pops[4], finvalue[i]
    }
}

