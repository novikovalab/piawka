@namespace "list"

function run() { 
  help="\
    List all statistics available for calculation with `piawka calc`. Takes no input."
  arg::parse_args(1, help)

  arg::args["dependencies"]=1
  print stats::format_stats()
  exit 0
}

