#!/bin/sh
"exec" "gawk" -f "$0" "--" "$@" && 0 {}

# Convert piawka Dxy values into NEXUS distance matrix.
# Input is piawka output (only lines with Dxy).
# If loci are many, weighted average is taken.
# This version produces output compatible with:
#   R function phangorn::write.nexus.dist()
#   SplitsTree
#   spectre netmake
#
# Assumptions about input data:
#  - only rows with Dxy present (either `piawka --nopi` or `grep Dxy piawka_output`)
#  - each comparison of groups has some data in the output

BEGIN{
  OFS="\t"
}

$2 != "nSites" {
  dist[$3,$4]+=( $9 ? $8 : ($7*$5) )
  scal[$3,$4]+=( $9 ? $9 : $5 )
}

END {
  for (i in dist) {
    s=index(i, SUBSEP)
    s1=substr(i,1,s-1)
    s2=substr(i,s+1)
    if (!(s1 in wasfirst)) { nsamples++; wasfirst[s1]=0 }
    if (!(s2 in wasfirst)) { nsamples++; wasfirst[s2]=0 }
    wasfirst[s1]++
  }
  for (i in wasfirst) { sortfirst[wasfirst[i]]=i }
  print "#NEXUS"
  print ""
  print "BEGIN TAXA;"
  print "\tDIMENSIONS ntax="nsamples";"
  printf "\tTAXLABELS";
  for (i=0;i<nsamples;i++) { printf " %s", sortfirst[i] }
  print " ;"
  print "END;"
  print ""
  print "BEGIN DISTANCES;"
  print "\tDIMENSIONS ntax="nsamples";" # needed for spectre to work
  print "\tFORMAT TRIANGLE = LOWER LABELS = LEFT;" # labels=left needed for spectre to work
  print "\tMatrix"
  for(i=0;i<nsamples;i++) {
    printf("\t%s", sortfirst[i])
    for(j=0;j<i;j++) { 
      printf(" %.6f", dist[sortfirst[i],sortfirst[i-j-1]] / scal[sortfirst[i],sortfirst[i-j-1]]) 
    }
    print " 0.000000"
  }
  print "\t;"
  print "END;"
}
