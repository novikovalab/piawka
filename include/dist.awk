@namespace "dist"

function run(){
  help="\
    Convert `piawka calc` pairwise distance values as a PHYLIP/NEXUS matrix.\n\
    Input is piawka output with some pairwise statistics in.\n\
    By default, dxy is used. Can be changed to other statistic (fst, rho etc.)\n\
    If many values per sample pair are present, weighted average is taken.\n\
    EXAMPLE: \n\tpiawka dist [OPTIONS] file.bed > distmat.phy"
  arg::add_argument("s", "stat", 0, "piawka pairwise stat to use as distance")
  arg::add_argument("n", "nexus", 1, "write matrix in NEXUS format instead of PHYLIP")
  arg::parse_args(2, help, "no help if empty")
  narg=arg::parse_nonargs()
  if ( narg==0 ) {
    arg::nonargs[++narg]="/dev/stdin"
    calc::say("Warning: no input files given, reading from stdin!")
  }

  if (!( "stat" in arg::args )) {
    arg::args["stat"]="dxy"
  }
  exit calc2dist(arg::nonargs[1])
}

function calc2dist(f) {
  while (getline < f > 0) {
    if ( $7 == arg::args["stat"] ) {
      samples[$5]++
      samples[$6]++
      dist[$5,$6]+=$9
      scal[$5,$6]+=$10
    }
  }
  piawka::assert( awk::isarray(samples), "lines with stat "arg::args["stat"]" not found in file "f)
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
