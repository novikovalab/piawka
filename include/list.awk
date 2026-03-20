@namespace "list"

function run() { 
  help="\
    List all statistics available for calculation with `piawka calc`. Takes no input. \n\
    EXAMPLE: \n\tpiawka list [OPTIONS]"
  arg::add_argument("d", "dependencies", 1, "output dependencies stats as well (best for piping to `piawka sum`)")
  arg::parse_args(2, help, "no help if empty")

  print stats::format_stats()
  exit 0
}

