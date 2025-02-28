#!/bin/sh
"exec" "gawk" -f "$0" "--" "$@" && 0 {}

BEGIN{
help="\
    Convert piawka Dxy values into NEXUS distance matrix.\n\
    Input is piawka output.\n\
    If loci are many, weighted average is taken.\n\
    By default, Dxy is used. Can be changed to other metric like\n\
    EXAMPLE: \n\tpiawka2dist.awk [OPTIONS] file.tsv > distance_file\n\
    OPTIONS:\n\
    METRIC=Dxy   -- piawka metric to use as distance\n\
    OUT=nexus    -- output format (nexus or phylip)"
  if (ARGC == 1) {print help; quiet_exit=1; exit 1}
  OFS="\t"
  METRIC="Dxy"
  OUT="nexus"
}

$6 == METRIC {
  samples[$3]++
  samples[$4]++
  dist[$3,$4]+=( $9 ? $8 : ($7*$5) )
  scal[$3,$4]+=( $9 ? $9 : $5 )
}

END {
  if ( quiet_exit ) { exit }
    # Nexus output is compatible with:
    #   R function phangorn::write.nexus.dist()
    #   SplitsTree
    #   spectre netmake
  if (OUT=="nexus") {
    print "#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS ntax="length(samples)";\n\tTAXLABELS"
    for (i in samples) { printf " %s", i }
    print " ;\nEND;\n\nBEGIN DISTANCES;\n\tDIMENSIONS ntax="length(samples)";\n\tFORMAT TRIANGLE = BOTH LABELS = LEFT;\n\tMatrix"
  } else if (OUT=="phylip") { print length(samples) }
  for(i in samples) {
    dist[i,i]=0
    scal[i,i]=1
    if (OUT=="nexus") { printf "\t" }
    printf( "%s", i)
    for(j in samples) { 
      if (i SUBSEP j in scal) {
        printf(" %.6f", dist[i,j] / scal[i,j] )
      } else if (j SUBSEP i in scal) {
        printf(" %.6f", dist[j,i] / scal[j,i] )
      } else {
        printf(" NaN")
      }
    }
    printf "\n"
  }
  if (OUT=="nexus") { print "\n\t;\nEND;" }
}
