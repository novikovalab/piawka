@namespace "stats"

function parse_stats(s) {
  split( tolower(s), statargs, "," )
  for (i in statargs) { 
    piawka::assert( statargs[i] in statslist, "statistic not found: "statargs[i]" (check options with piawka -l)" )
    nr=statslist[statargs[i]]
    is_between=( stat_isbetween[nr] ? "between" : "within" )
    stats_print[is_between][statargs[i]]=1 
    stats[is_between][statargs[i]]=1 
    if (stat_dep[nr] != "") {
      ndeps=split(stat_dep[nr], stat_deps, ",")
      for (j=1;j<=ndeps;j++) {
        piawka::assert( stat_deps[j] in statslist, "dependent statistic not found: "statargs[i]" (check options with piawka -l)" )
        nr=statslist[stat_deps[j]]
        is_between=( stat_isbetween[nr] ? "between" : "within" )
        stats[is_between][stat_deps[j]]=1
      }
    }
  }
}

function add_stat(name, desc, is_between, dep) {
  nstat++
  statslist[name]=nstat
  stat_name[nstat] = name
  stat_desc[nstat] = desc (is_between ? " (within population)" : " (between populations)")
  stat_isbetween[nstat] = is_between
  stat_dep[nstat] = dep
}

function format_stats(   help, help_col1, total_width) {
  total_width=80
  for (i=1;i<=nstat;i++) {
    if (stat_desc[i] ~ /^helper:/) { continue }
    help_col1[i] = stat_name[i] 
    if ( (l=length(help_col1[i])) > help_col1_width ) { help_col1_width=l }
  }
  help_col1_width++
  remain_width = total_width - help_col1_width
  help="Available statistics:"
  for (i=1;i<=nstat;i++) {
    if (stat_desc[i] ~ /^helper:/) { continue }
    nlines=1
    if (length(stat_desc[i]) > remain_width ) {
      for (j=remain_width+index(stat_desc[i], "\n"); j<length(stat_desc[i]); j+=remain_width+1) {
        for (k=0; k<remain_width/2; k++) {
          if (substr(stat_desc[i], j-k, 1) ~ /[ ,)]/) { break }
        }
        if (k < remain_width/2-1) { j-=k }
        stat_desc[i] = substr(stat_desc[i], 1, j) "\n" substr(stat_desc[i],j+1)
      }
      nlines=split(stat_desc[i], desclines, "\n")
    }
    help = help "\n" sprintf("%-"help_col1_width"s", help_col1[i]) ( nlines > 1 ? desclines[1] : stat_desc[i] )
    for (j=2;j<=nlines;j++) {
      help = help "\n" sprintf("%-"help_col1_width"s", " ") desclines[j] 
    }
  }
  return help
}
