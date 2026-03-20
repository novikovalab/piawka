@namespace "win"

function run() { 
  help="\
    Chunk a VCF file into windows for parallel processing or sub-chromosome resolution.\n\
    Default (slow): output windows containing --lines [100000] VCF lines each. \n\
    Alternatively, output windows of --window-size at --step \n\
    (requires --vcf with ##contig header lines or a --fai file). \n\
    If --targets are supplied instead of --vcf, produces bigger windows with multiple targets in. \n\
    EXAMPLE: \n\tpiawka win [OPTIONS] -v file.vcf.gz | -T targets.bed > win.bed"
  arg::add_argument("f", "fai", 0, "FAI index file to get chr sizes from (any TSV with chr & size as cols 1 & 2); prioritized over --vcf and --targets")
  arg::add_argument("l", "lines", 0, "# VCF lines per window")
  arg::add_argument("n", "number", 1, "enumerate windows in 4th column")
  arg::add_argument("s", "step", 0, "step in bp for sliding windows (requires --window-size)")
  arg::add_argument("T", "targets", 0, "BED file with small regions to aggregate into larger ones")
  arg::add_argument("v", "vcf", 0, "gzipped and tabixed VCF file; prioritized over --targets")
  arg::add_argument("w", "window-size", 0, "bp per window (requires ##contig lines in VCF header)")
  arg::parse_args(2, help)
  DEFAULT_WINDOWSIZE=100000
  piawka::assert( !( "lines" in arg::args && "window-size" in arg::args), "options --lines and --window-size are mutually exclusive!" )
  if ( "window-size" in arg::args && "step" in arg::args ) {
    calc::say("Warning: --step implies fixed-length windows with --window-size")
  }
  exit main()
}

function main() {
  if ( "fai" in arg::args ) {
    piawka::check_file(arg::args["fai"])
    make_windows_bp_fai( "window-size" in arg::args ? arg::args["window-size"] : DEFAULT_WINDOWSIZE )
  } else if ("vcf" in arg::args) {
    piawka::check_file(arg::args["vcf"])
    if ( !( "window-size" in arg::args || "step" in arg::args ) ) {
      make_windows( "lines" in arg::args ? arg::args["lines"] : DEFAULT_WINDOWSIZE )
    } else {
      make_windows_bp( "window-size" in arg::args ? arg::args["window-size"] : DEFAULT_WINDOWSIZE )
    }
  } else if ( "targets" in arg::args ) {
    piawka::check_file(arg::args["targets"])
    aggregate_regions()
  } else {
    calc::say("Error: no FAI, VCF, or BED file given!")
    return 1
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
  cmd = "tabix "arg::args["vcf"]" $( tabix -l "arg::args["vcf"]" )"
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

function make_windows_bp( window_size ) {
  cmd = "tabix -H "arg::args["vcf"]
  step = ( arg::args["step"] == "" ? window_size : arg::args["step"] )
  while ( cmd | getline > 0 ) {
    if ( /^##contig/ ) {
      if ( match($0,/ID=[^,>]+/) ) {
        id=substr($0, RSTART+3, RLENGTH-3)
      } else { continue }
      if ( match($0,/length=[0-9]+/) ) {
        len=int(substr($0, RSTART+7, RLENGTH-7))
      } else { continue }
      start=end=0
      while (end < len) {
        end = start+window_size
        print id"\t"start"\t"( end > len ? len : end ) ( (arg::args["number"]==1) ? "\twin"++window_number : "" )
        start += step 
      }
    }
  }
  close(cmd)
  if (start=="") {
    calc::say("Error: VCF header does not carry ##contig lines with ID and length values")
    exit 1
  }
}

function make_windows_bp_fai( window_size ) {
  step = ( arg::args["step"] == "" ? window_size : arg::args["step"] )
  FS="\t"
  while ( getline < arg::args["fai"] > 0 ) {
    id=$1
    len=int($2)
    start=end=0
    while (end < len) {
      end = start+window_size
      print id"\t"start"\t"( end > len ? len : end ) ( (arg::args["number"]==1) ? "\twin"++window_number : "" )
      start += step 
    }
  }
  close(arg::args["fai"])
  exit 0
}
