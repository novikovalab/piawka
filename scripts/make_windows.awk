#!/bin/sh
"export" "AWKPATH=$( dirname $( which piawka ) ):$AWKPATH"
"exec" "gawk" -f "$0" "--" "$@" && 0 {}
#
# Convert VCF to BED, each region has a certain #lines 

@load "filefuncs"

@include "getopt.awk"
@include "argparse.awk"

BEGIN{ 
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
  argparse::add_argument("v", "vcf", 0, "gzipped and tabixed VCF file")
  argparse::add_argument("w", "windowsize", 0, "# VCF lines per window")
  help = usage "\n" argparse::format_help()
  getopt::Optind = 1 # start with 1st opt
  getopt::Opterr = 1 # print getopt errs
  if ( argparse::parse_args() != 0 ) { print help; return 0 }
  for(i in argparse::args) { args[i]=argparse::args[i] } # shorten args array name
  if (args["help"] || ARGC==1 ) { print help; return 0 }
  check_file(args["vcf"])
  if ( args["windowsize"] == "" ) { args["windowsize"]=DEFAULT_WINDOWSIZE }
  make_windows(args["windowsize"])
  exit 0
}

function check_file(file,   cmd) {
   if (stat(file, _) < 0) {
        print "Error: could not read file "file > "/dev/stderr"
        exit 1
    }
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
    }
    # If chromosome changes or #lines > size, print a window
    if ( $1 != last_chr || ++lines_seen == chunk_size ) {
      print last_chr, minpos-1, lastpos,
             (args["number"]==1) ? "win"++window_number : args["vcf"]
      fflush()
      last_chr=$1
      minpos=$2
      lines_seen=0
    }
    lastpos=$2
  }
  if (lines_seen > 0) { 
    print last_chr, minpos-1, lastpos,
           (args["number"]==1) ? "win"++window_number : args["vcf"]
  }
  close(cmd)
}
