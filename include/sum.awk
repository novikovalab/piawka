@namespace "sum"

function run() { 
  help="\
    This script summarizes piawka calc results counted over several loci. \n\
    It simply divides sum of numerators by sum of denominators. \n\
    The only arguments are the output file(-s) of piawka calc; stdin should be passed as `piawka sum -`."
  arg::parse_args(2, help)
  narg=arg::parse_nonargs()
  for (n=1; n<=narg; n++) {
    summarize_regions(arg::nonargs[n])
  }
  exit 0
}

function summarize_regions(f,    firstline) {
  firstline=1
  while (getline < f > 0) {
    if (firstline) { 
      firstline=0
      piawka::assert( NF == 13, "piawka output (13 columns) is required!" )
    }
    if ( $9 != "value" && $9==$9 ) { # exclude header and NaNs 
      idx=$1 SUBSEP $4 SUBSEP $5 SUBSEP $6 SUBSEP $8 SUBSEP
      if (!seen[idx]) {
        seen[idx]=1
        start[idx]=$2
        end[idx]=$3
      }
      if ($2 < start[idx]) { start[idx]=$2 }
      if ($3 > end[idx]) { end[idx]=$2 }
      allnUsed[idx]+=$7
      denominator[idx]+=$11
      numerator[idx]+=$10
      nGeno[idx]+=$12
      nMiss[idx]+=$13
    }
  }
  for (i in seen) {
    split(i, pops, SUBSEP)
    finvalue[i]=numerator[i]/denominator[i]"\t"numerator[i]"\t"denominator[i]"\t"nGeno[i]"\t"nMiss[i]
    print pops[1]"\t"start[i]"\t"end[i]"\t"pops[2]"\t"pops[3]"\t"pops[4]"\t"allnUsed[i]"\t"pops[5]"\t"finvalue[i]
  }
}

