@namespace "dist"

function run(){
  help="\
    Convert piawka Dxy values into a distance matrix.\n\
    Input is piawka output with some pairwise statistic in.\n\
    By default, Dxy is used. Can be changed to other metric (fst, rho etc.)\n\
    If loci are many, weighted average is taken.\n\
    EXAMPLE: \n\tpiawka2dist.awk [OPTIONS] file.tsv > distance_file"
  argparse::add_argument("h", "help", 1, "show this help message")
  argparse::add_argument("m", "metric", 0, "piawka pairwise metric to use as distance")
  argparse::add_argument("n", "nexus", 1, "write matrix in NEXUS format instead of PHYLIP")
  help = usage "\n" argparse::format_help()
  getopt::Optind = 2 # start with 2nd opt
  getopt::Opterr = 1 # print getopt errs
  if ( argparse::parse_args() != 0 ) { print help; exit 0 }
  for(i in argparse::args) { args[i]=argparse::args[i] } # shorten args array name
  if ( args["help"] || ARGC==2 ) { print help; exit 0 }
  OFS="\t"
  if (args["metric"]=="") {
    args["metric"]="dxy"
  }
  exit calc2dist(ARGV[ARGC-1])
}

function calc2dist(f) {
  while (getline < f > 0) {
    if ( $8 == args["metric"] ) {
      samples[$5]++
      samples[$6]++
      dist[$5,$6]+=$10
      scal[$5,$6]+=$11
    }
  }
  piawka::assert( awk::isarray(samples), "lines with metric "args["metric"]" not found in file "f)
    # Nexus output is compatible with:
    #   R function phangorn::write.nexus.dist()
    #   SplitsTree
    #   spectre netmake
  if (args["nexus"]) {
    print "#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS ntax="length(samples)";\n\tTAXLABELS"
    for (i in samples) { printf " %s", i }
    print " ;\nEND;\n\nBEGIN DISTANCES;\n\tDIMENSIONS ntax="length(samples)";\n\tFORMAT TRIANGLE = BOTH LABELS = LEFT;\n\tMatrix"
  } else { print length(samples) }
  for(i in samples) {
    dist[i,i]=0
    scal[i,i]=1
    if (args["nexus"]) { printf "\t" }
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
  if (args["nexus"]) { print "\n\t;\nEND;" }
  return 0
}
