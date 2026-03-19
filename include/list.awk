@namespace "list"

function run() { 
  help="\
    List all statistics available for calculation with `piawka calc`. Takes no input. \n\
    EXAMPLE: \n\tpiawka list [OPTIONS]"
  arg::parse_args(1, help)

  arg::args["dependencies"]=1
  print stats::format_stats()
  exit 0
}

