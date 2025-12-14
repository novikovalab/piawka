# An interface to getopt.awk to simplify adding args and writing help messages.
# taprs, 2025-11-24

@include "getopt.awk"

@namespace "arg"

function add_argument(short, long, is_flag, desc) {
  narg++
  args_short[narg] = short
  args_long[narg] = long
  args_isflag[narg] = is_flag
  args_desc[narg] = desc
  shortopts = shortopts short (is_flag ? "" : ":")
  longopts = longopts long (is_flag ? "" : ":") ","
}

function format_help(   help, help_col1, total_width) {
  total_width=80
  for (i=1;i<=narg;i++) {
    help_col1[i] = "-"args_short[i]", --"args_long[i] (args_isflag[i] ? "" : " <arg>")
    if ( (l=length(help_col1[i])) > help_col1_width ) { help_col1_width=l }
  }
  help_col1_width++
  remain_width = total_width - help_col1_width
  help="Options:"
  for (i=1;i<=narg;i++) {
    nlines=1
    if (length(args_desc[i]) > remain_width ) {
      for (j=remain_width+index(args_desc[i], "\n"); j<length(args_desc[i]); j+=remain_width+1) {
        for (k=0; k<remain_width/2; k++) {
          if (substr(args_desc[i], j-k, 1) ~ /[ ,)]/) { break }
        }
        if (k < remain_width/2-1) { j-=k }
        args_desc[i] = substr(args_desc[i], 1, j) "\n" substr(args_desc[i],j+1)
      }
      nlines=split(args_desc[i], desclines, "\n")
    }
    help = help "\n" sprintf("%-"help_col1_width"s", help_col1[i]) ( nlines > 1 ? desclines[1] : args_desc[i] )
    for (j=2;j<=nlines;j++) {
      help = help "\n" sprintf("%-"help_col1_width"s", " ") desclines[j] 
    }
  }
  return help
}

function parse_args(firstopt, helpmsg, helpltr) {
  if (firstopt=="") { firstopt=1 }
  getopt::Optind = firstopt # start with 1st opt
  getopt::Opterr = 1 # print getopt errs

  if (helpltr=="") { helpltr="h" }
  add_argument(helpltr, "help", 1, "show this help message and exit")

  while ((c=getopt::getopt(ARGC,ARGV,shortopts,longopts)) != -1) {
    if (c == "?") { 
      print helpmsg
      exit 1
    }
    for (i=1;i<=narg;i++) {
      if (getopt::Optopt == args_short[i] || getopt::Optopt == args_long[i]) {
        args[args_long[i]]=(args_isflag[i] ? 1 : getopt::Optarg)
        break
      }
    }
  }
  help = helpmsg "\n" format_help()
  if ( args["help"] || ARGC <= firstopt ) { 
    print help
    exit 0 
  }
  return 0
}

function parse_nonargs(    counter) {
  for (; getopt::Optind < ARGC; getopt::Optind++) {
    nonargs[++counter]=ARGV[getopt::Optind]
  }
  return counter
}
