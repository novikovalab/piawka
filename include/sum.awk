@namespace "sum"

function run() { 
  help="\
    This script summarizes piawka calc results counted over several loci. \n\
    It simply divides sum of numerators by sum of denominators. \n\
    The only arguments are the output file(-s) of piawka calc; stdin should be passed as `piawka sum -`."
  arg::add_argument("s", "stats", 0, "recalculate stats that cannot be summarized as sum(num)/sum(den) using dependencies, format as in piawka calc")
  arg::add_argument("i", "ignore-chrs", 1, "summarize statistics across all chromosomes using only locus field to match")
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
      if ( "ignore-chrs" in arg::args ) {
        $1=$2=$3="_"
      }
      if ( !("stats" in arg::args) && !seenstat[$7] ) {
        seenstat[$7]=1
        stats=stats","$7
      }
      i1=$1 SUBSEP $4                      # locus-chr
      i2=$5 ( $6 == "." ? "" : SUBSEP $6 ) # pops
      i3= $7                               # metric
      if ($2 < start[ i1 ] || start[ i1 ]=="") { start[ i1 ]=$2 }
      if ($3 > end[ i1 ]) { end[ i1 ]=$3 }
      den[ i1 ][ i2 ][ i3 ]+=$10
      num[ i1 ][ i2 ][ i3 ]+=$9
    }
  }
  close(f)
  if (stats != "") {
    stats::parse_stats(substr(stats,2))
  }
  calc::tmpf="/dev/stdout"
  for (i1 in den) {
    split(i1, ii, SUBSEP)
    calc::chr=ii[1]
    calc::locus=ii[2]
    calc::start=start[i1]
    calc::end=end[i1]
    piawka::copy_array(num[i1],calc::num)
    piawka::copy_array(den[i1],calc::den)
    calc::yield_output() # if dependencies are not there, falls back to sum(num)/sum(den) because finalize_i3 returns nothing
  }
}

