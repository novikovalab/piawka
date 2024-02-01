#!/usr/bin/mawk -f
#
# This script summarizes piawka results counted over several blocks.
# It takes mean value weighted by nUsed across all lines with same locus, pop1 and pop2.

BEGIN{ OFS="\t" }

{ 
  nSites[$1,$3,$4]+=$2
  nUsed[$1,$3,$4]=$5" "nUsed[$3,$4]
  allnUsed[$1,$3,$4]+=$5
  metric[$1,$3,$4]=$6
  value[$1,$3,$4]=$7" "value[$3,$4]
}

END{
  for (i in nSites) {
    split(i, pops, SUBSEP)
    split(nUsed[i], weights, " ")
    split(value[i], values, " ")
    for (j in weights) { finvalue[i] += ( values[j] * weights[j] ) / allnUsed[i] }
    print pops[1], nSites[i], pops[2], pops[3], allnUsed[i], metric[i], finvalue[i]
    }
}

