#!/bin/sh
"exec" "gawk" -f "$0" "--" "$@" && 0 {}

# Convert piawka Dxy values into NEXUS distance matrix.
# Input is piawka output.
# If loci are many, weighted average is taken.
# This version produces output compatible with:
#   R function phangorn::write.nexus.dist()
#   SplitsTree
#   spectre netmake
#
# By default, Dxy is used. Can be changed to other metric like
# piawka2nexus.awk METRIC=Fst_HUD file.tsv

BEGIN{
  OFS="\t"
  METRIC="Dxy"
}

$6 == METRIC {
  samples[$3]++
  samples[$4]++
  dist[$3,$4]+=( $9 ? $8 : ($7*$5) )
  scal[$3,$4]+=( $9 ? $9 : $5 )
}

END {
  print "#NEXUS"
  print ""
  print "BEGIN TAXA;"
  print "\tDIMENSIONS ntax="length(samples)";"
  printf "\tTAXLABELS";
  for (i in samples) { printf " %s", i }
  print " ;"
  print "END;"
  print ""
  print "BEGIN DISTANCES;"
  print "\tDIMENSIONS ntax="length(samples)";" # needed for spectre to work
  print "\tFORMAT TRIANGLE = BOTH LABELS = LEFT;" # labels=left needed for spectre to work
  print "\tMatrix"
  for(i in samples) {
    dist[i,i]=0
    scal[i,i]=1
    printf("\t%s", i)
    for(j in samples) { 
      if (i SUBSEP j in scal) {
        printf(" %.6f", dist[i,j] / scal[i,j] )
      } else if (j SUBSEP i in scal) {
        printf(" %.6f", dist[j,i] / scal[j,i] )
      } else {
        printf(" NaN")
      }
    }
  }
  printf "\n"
  print "\t;"
  print "END;"
}
