#!/usr/bin/mawk -f
#
# This script summarizes piawka results counted over several blocks.
# It takes mean value weighted by nUsed across all lines with same locus, pop1, pop2 and metric.

BEGIN{ OFS="\t" }

{ 
  nSites[$1,$3,$4,$6]+=$2
  nUsed[$1,$3,$4,$6]=$5" "nUsed[$1,$3,$4,$6]
  allnUsed[$1,$3,$4,$6]+=$5
  #metric[$1,$3,$4,$6]=$6
  value[$1,$3,$4,$6]=$7" "value[$1,$3,$4,$6]
}

END{
  for (i in nSites) {
    split(i, pops, SUBSEP)
    split(nUsed[i], weights, " ")
    split(value[i], values, " ")
    for (j in weights) { finvalue[i] += ( values[j] * weights[j] ) / allnUsed[i] }
    print pops[1], nSites[i], pops[2], pops[3], allnUsed[i], pops[4], finvalue[i]
    }
}

