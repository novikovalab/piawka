@namespace "sum"

function run() { 
  help="\
    Summarize `piawka calc` results counted over several loci. \n\
    If dependencies are not given (see `piawka list`), defaults to sum(numerator)/sum(denominator). \n\
    It only takes the output file(-s) of `piawka calc` passed over stdin. \n\
    EXAMPLE: \n\tpiawka sum [OPTIONS] file.bed > file_sum.bed"
  arg::add_argument("g", "groups", 0, "group file to average stats across individuals/subgroups")
  arg::add_argument("i", "ignore-chrs", 1, "summarize statistics across all chromosomes using only locus field to match")
  arg::add_argument("s", "stats", 0, "stats to be summarized, defaults to all stats found (see `piawka list`)")
  arg::parse_args(2, help, "no help if empty")
  narg=arg::parse_nonargs()
  if ( narg==0 ) {
    arg::nonargs[++narg]="/dev/stdin"
    calc::say("Warning: no input files given, reading from stdin!")
  }
  if ( "stats" in arg::args ) {
    stats::parse_stats(arg::args["stats"])
  }
  if ( "groups" in arg::args ) {
    calc::get_groups()
  }
  calc::print_header()
  calc::tmpf="/dev/stdout"
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

      # print if new locus/chr
      if ( locuschr != $1 SUBSEP $4 && locuschr != "" ) {
        print_sum()
      }

      if ( "ignore-chrs" in arg::args ) {
        $1=$2=$3="_"
      }
      if ( !("stats" in arg::args) && !seenstat[$7] ) {
        seenstat[$7]=1
        stats=stats","$7
        stats::parse_stats(substr(stats,2))
      }
      if ( "groups" in arg::args ) {
        if ( $5 in calc::groupmem && ($6=="." || $6 in calc::groupmem) ) {
          $5=calc::groupmem[$5]
          $6=( $6=="." ? "." : calc::groupmem[$6] )
        } else { 
          continue 
        }
      }
      locuschr=$1 SUBSEP $4                # locus-chr
      i2=$5 ( $6 == "." ? "" : SUBSEP $6 ) # pops
      i3= $7                               # metric
      if ($2 < calc::start || calc::start=="") { calc::start=$2 }
      if ($3 > calc::end) { calc::end=$3 }
      calc::den[ i2 ][ i3 ]+=$10
      calc::num[ i2 ][ i3 ]+=$9
    }
  }
  print_sum()
  close(f)
}

function print_sum() {
  split(locuschr, ii, SUBSEP)
  calc::chr=ii[1]
  calc::locus=ii[2]
  calc::yield_output() # if dependencies are not there, falls back to sum(num)/sum(den) because finalize_i3 returns nothing
  delete calc::num
  delete calc::den
  calc::start=calc::end=""
}
