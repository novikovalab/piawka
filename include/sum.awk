@namespace "sum"

# This script summarizes piawka results counted over several loci.
# It simply divides sum of numerators by sum of denominators.

function run() { 
  OFS="\t" 
  if (ARGC<=2) {
    summarize_regions("/dev/stdin")
  } else {
    for (a=2; a<ARGC; a++) {
      summarize_regions(ARGV[a])
    }
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
    print pops[1], start[i], end[i], pops[2], pops[3], pops[4], allnUsed[i], pops[5], finvalue[i]
  }
}

