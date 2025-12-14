# getopt.awk --- Do C library getopt(3) function in awk
#                Also supports long options.
#
# From: https://www.gnu.org/software/gawk/manual/html_node/Getopt-Function.html

# External variables:
#    Optind -- index in ARGV of first nonoption argument
#    Optarg -- string value of argument to current option
#    Opterr -- if nonzero, print our own diagnostic
#    Optopt -- current option letter

# Returns:
#    -1     at end of options
#    "?"    for unrecognized option
#    <s>    a string representing the current option

# Private Data:
#    _opti  -- index in multiflag option, e.g., -abc

@include "getopt.awk"

@namespace "argparse"

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

function parse_args() {
  while ((c=getopt::getopt(ARGC,ARGV,shortopts,longopts)) != -1) {
    if (c == "?") { return 1 }
    for (i=1;i<=narg;i++) {
      if (getopt::Optopt == args_short[i] || getopt::Optopt == args_long[i]) {
        args[args_long[i]]=(args_isflag[i] ? 1 : getopt::Optarg)
        break
        }
      }
  }
  return 0
}

