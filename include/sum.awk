@namespace "sum"

function run() { 
  help="\
    This script summarizes piawka calc results counted over several loci. \n\
    It simply divides sum of numerators by sum of denominators. \n\
    The only arguments are the output file(-s) of piawka calc; stdin should be passed as `piawka sum -`."
  arg::add_argument("s", "stats", 0, "recalculate stats that cannot be summarized as sum(num)/sum(den) using dependencies, format as in piawka calc")
  arg::parse_args(2, help)
  narg=arg::parse_nonargs()

  if ( "stats" in arg::args ) {
    stats::parse_stats(arg::args["stats"])
  }
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
      piawka::assert( NF == 10, "piawka output (10 columns) is required!" )
    }
    if ( substr($0,1,1)!="#" && $8==$8 ) { # exclude header and NaNs 
      idx=$1 SUBSEP $4 SUBSEP $5 SUBSEP $6 SUBSEP $7 SUBSEP
      if (!seen[idx]) {
        seen[idx]=1
        start[idx]=$2
        end[idx]=$3
      }
      if ($2 < start[idx]) { start[idx]=$2 }
      if ($3 > end[idx]) { end[idx]=$2 }
      if (stats::summary_func[$7]=="sum") {
        numerator[idx]+=$9/$10
        denominator[idx]=1
      } else {
        denominator[idx]+=$10
        numerator[idx]+=$9
      }
    }
  }
  for (i in seen) {
    split(i, pops, SUBSEP)
    finvalue[i]=numerator[i]/denominator[i]"\t"numerator[i]"\t"denominator[i]
    print pops[1]"\t"start[i]"\t"end[i]"\t"pops[2]"\t"pops[3]"\t"pops[4]"\t"pops[5]"\t"finvalue[i]
    # populate num[locus][pop1,pop2][metric] array
    if ("stats" in arg::args) {
      allnum[ pops[1],pops[2] ][ pops[3] (pops[4]=="." ? "" : SUBSEP pops[4]) ][ pops[5] ]=numerator[i]
      allden[ pops[1],pops[2] ][ pops[3] (pops[4]=="." ? "" : SUBSEP pops[4]) ][ pops[5] ]=denominator[i]
      allstart[ pops[1],pops[2] ] = start[i]
      allend[ pops[1],pops[2] ] = end[i]
    }
  }
  if ("stats" in arg::args) {
    calc::tmpf="/dev/stdout"
    for (i in allden) {
      split(i, ii, SUBSEP)
      calc::chr=ii[1]
      calc::start=allstart[i]
      calc::end=allend[i]
      calc::locus=ii[2]
      copy_array(allnum[i],calc::num)
      copy_array(allden[i],calc::den)
      calc::yield_output()
    }
  }
}

# Array copy func by Ed Morton @ https://stackoverflow.com/a/62179751/14993291
function copy_array(orig, copy,      i) {
  for (i in orig) {
    if (awk::isarray(orig[i])) {
      copy_array(orig[i], copy[i])
    }
    else {
      copy[i] = orig[i]
    }
  }
}
