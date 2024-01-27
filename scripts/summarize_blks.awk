#!/usr/bin/awk -f
#
# This script summarizes pi and dxy values calculated for VCF file blocks by `piawka_par_blk.sh`.

BEGIN{ OFS="\t" }

{ 
  if (!FLAG) {locus=$1; FLAG=1}
  nSites[$3,$4]+=$2
  nUsed[$3,$4]=$5" "nUsed[$3,$4]
  allnUsed[$3,$4]+=$5
  metric[$3,$4]=$6
  value[$3,$4]=$7" "value[$3,$4]
}

END{
  for (i in nSites) {
    split(i, pops, SUBSEP)
    split(nUsed[i], weights, " ")
    split(value[i], values, " ")
    for (j in weights) { finvalue[i] += ( values[j] * weights[j] ) / allnUsed[i] }
    print locus, nSites[i], pops[1], pops[2], allnUsed[i], metric[i], finvalue[i]
    }
}

