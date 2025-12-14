@namespace "win"

function run() { 
  help="Usage:\nmake_windows.awk -v vcf_gz [OPTIONS]"
  arg::add_argument("n", "number", 1, "enumerate windows in 4th column")
  arg::add_argument("s", "skip", 0, "# lines to skip")
  arg::add_argument("T", "targets", 0, "BED file with small regions to aggregate into larger ones")
  arg::add_argument("v", "vcf", 0, "gzipped and tabixed VCF file")
  arg::add_argument("w", "windowsize", 0, "# VCF lines per window")
  arg::parse_args(2, help)
  DEFAULT_WINDOWSIZE=10000
  if ( arg::args["windowsize"] == "" ) { arg::args["windowsize"]=DEFAULT_WINDOWSIZE }
  piawka::assert( awk::xor(arg::args["vcf"]!="", arg::args["targets"]!=""), "provide either VCF or BED to generate windows from" )
  exit main()
}

function main() {
  if (arg::args["vcf"]!="") {
    piawka::check_file(arg::args["vcf"])
    make_windows(arg::args["windowsize"])
  } else if ( arg::args["targets"]!="" ) {
    piawka::check_file(arg::args["targets"])
    aggregate_regions()
  }
}

function aggregate_regions() {
# This script acts as `bedtools merge -d 1000`, joining BED records closer than 1kb to each other.
# This is done to reduce the number of tabix queries with numerous, small, adjacent BED regions.
# The aggregated file is then passed as -R argument and the small regions as -T argument.
  firstline=1
  while (getline < arg::args["targets"] > 0) {
    if (firstline) { 
      firstline=0
      printf $1"\t"$2"\t"
      lastchr=$1
      lastend=$3
      continue 
    }
    if ( $1!=lastchr || $2 - lastend > 1000 ) { 
      print lastend
      printf $1"\t"$2"\t" 
    }
    lastchr=$1; lastend=$3
  }
  print lastend
}

function make_windows(chunk_size,    
                      cmd, lines_seen, first_window,
                      last_chr, minpos, lastpos) {
  cmd = "bgzip -dc "arg::args["vcf"]
  if (arg::args["skip"]>0) {
    while( arg::args["skip"]-- > 0 ) { cmd | getline }
  }
  first_window=1
  while ( cmd | getline > 0 ) {
    if ( passed_header==0 && index($0,"#")==1 ) { continue }
    passed_header=1
    if ( first_window==1 ) {
      first_window=0
      last_chr=$1
      minpos=$2
      lastpos=$2
      lines_seen=1
    }
    if ( $2!=lastpos ) {
      lines_seen++
    }
    # If chromosome changes or #lines > size, print a window
    if ( $1 != last_chr || lines_seen == chunk_size ) {
      print last_chr"\t"minpos-1"\t"lastpos ( (arg::args["number"]==1) ? "\twin"++window_number : "" )
      fflush()
      last_chr=$1
      minpos=$2
      lines_seen=1
    }
    lastpos=$2
  }
  if (lines_seen > 0) { 
    print last_chr"\t"minpos-1"\t"lastpos ( (arg::args["number"]==1) ? "\twin"++window_number : "" )
  }
  close(cmd)
}
