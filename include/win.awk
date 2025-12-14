@namespace "win"

function run() { 
  OFS="\t"
  usage="Usage:\nmake_windows.awk -v vcf_gz [OPTIONS]"
  DEFAULT_WINDOWSIZE=10000
  exit main()
}

function main() {
  # add_argument args: shortopt, longopt, is_flag, description
  argparse::add_argument("h", "help", 1, "show this help message")
  argparse::add_argument("n", "number", 1, "enumerate windows in 4th column")
  argparse::add_argument("s", "skip", 0, "# lines to skip")
  argparse::add_argument("T", "targets", 0, "BED file with small regions to aggregate into larger ones")
  argparse::add_argument("v", "vcf", 0, "gzipped and tabixed VCF file")
  argparse::add_argument("w", "windowsize", 0, "# VCF lines per window")
  help = usage "\n" argparse::format_help()
  getopt::Optind = 2 # start with 1st opt
  getopt::Opterr = 1 # print getopt errs
  if ( argparse::parse_args() != 0 ) { print help; return 0 }
  for(i in argparse::args) { args[i]=argparse::args[i] } # shorten args array name
  if (args["help"] || ARGC==1 ) { print help; return 0 }
  if ( args["windowsize"] == "" ) { args["windowsize"]=DEFAULT_WINDOWSIZE }
  piawka::assert( awk::xor(args["vcf"]!="", args["targets"]!=""), "provide either VCF or BED to generate windows from" )
  if (args["vcf"]!="") {
    check_file(args["vcf"])
    make_windows(args["windowsize"])
  } else if ( args["targets"]!="" ) {
    check_file(args["targets"])
    aggregate_regions()
  }
}

function check_file(file,   cmd) {
   if (awk::stat(file, _) < 0) {
        print "Error: could not read file "file > "/dev/stderr"
        exit 1
    }
}

function aggregate_regions() {
# This script acts as `bedtools merge -d 1000`, joining BED records closer than 1kb to each other.
# This is done to reduce the number of tabix queries with numerous, small, adjacent BED regions.
# The aggregated file is then passed as -R argument and the small regions as -T argument.
  firstline=1
  while (getline < args["targets"] > 0) {
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
  cmd = "bgzip -dc "args["vcf"]
  if (args["skip"]>0) {
    while( args["skip"]-- > 0 ) { cmd | getline }
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
      print last_chr, minpos-1, lastpos,
             (args["number"]==1) ? "win"++window_number : ""
      fflush()
      last_chr=$1
      minpos=$2
      lines_seen=1
    }
    lastpos=$2
  }
  if (lines_seen > 0) { 
    print last_chr, minpos-1, lastpos,
           (args["number"]==1) ? "win"++window_number : ""
  }
  close(cmd)
}
