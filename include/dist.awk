@namespace "dist"

function run(){
  help="\
    Convert piawka Dxy values into a distance matrix.\n\
    Input is piawka output with some pairwise statistic in.\n\
    By default, Dxy is used. Can be changed to other metric (fst, rho etc.)\n\
    If loci are many, weighted average is taken.\n\
    EXAMPLE: \n\tpiawka2dist.awk [OPTIONS] file.tsv > distance_file"
  arg::add_argument("m", "metric", 0, "piawka pairwise metric to use as distance")
  arg::add_argument("n", "nexus", 1, "write matrix in NEXUS format instead of PHYLIP")
  arg::parse_args(2, help)

  if (arg::args["metric"]=="") {
    arg::args["metric"]="dxy"
  }
  exit calc2dist(ARGV[ARGC-1])
}

function calc2dist(f) {
  while (getline < f > 0) {
    if ( $8 == arg::args["metric"] ) {
      samples[$5]++
      samples[$6]++
      dist[$5,$6]+=$10
      scal[$5,$6]+=$11
    }
  }
  piawka::assert( awk::isarray(samples), "lines with metric "arg::args["metric"]" not found in file "f)
    # Nexus output is compatible with:
    #   R function phangorn::write.nexus.dist()
    #   SplitsTree
    #   spectre netmake
  if (arg::args["nexus"]) {
    print "#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS ntax="length(samples)";\n\tTAXLABELS"
    for (i in samples) { printf " %s", i }
    print " ;\nEND;\n\nBEGIN DISTANCES;\n\tDIMENSIONS ntax="length(samples)";\n\tFORMAT TRIANGLE = BOTH LABELS = LEFT;\n\tMatrix"
  } else { print length(samples) }
  for(i in samples) {
    dist[i,i]=0
    scal[i,i]=1
    if (arg::args["nexus"]) { printf "\t" }
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
  if (arg::args["nexus"]) { print "\n\t;\nEND;" }
  return 0
}
