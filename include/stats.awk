@namespace "stats"

function parse_stats(commalist, dependency,    nstats,args,idx,i,j,s) {
  delete between
  delete within
  delete used
  delete printed
  nstats=split( tolower(commalist), args, "," )
  for (i=1;i<=nstats;i++) { 
    s=args[i]
    piawka::assert( s in list, (dependency ? "dependency" : "statistic" ) " not found: "s" (check options with `piawka list`)" )
    idx=list[s]
    used[s]=idx
    if (!dependency || "dependencies" in arg::args) { 
      printed[s]=idx
    }
    if ( isbetween[idx] == "within" ) {
      within[s]=idx
    } else if ( isbetween[idx] == "between" ) {
      between[s]=idx
    }
    parse_stats(dep[idx], "dependency")
  }
  if ( !dependency ) {
    n_between = awk::asort(between)
    n_within  = awk::asort(within)
    n_used = awk::asort(used)
    for (j in between) { between[j]=name[between[j]] }
    for (j in within) { within[j]=name[within[j]] }
    for (j in used) { used[j]=name[used[j]] }
  }
}

function add_stat(thisname, thisdesc, is_between, thisdep) {
  nstat++
  list[thisname]=nstat
  name[nstat] = thisname
  isbetween[nstat] = ( is_between!=0 ? "between" : "within" )
  dep[nstat] = thisdep
  desc[nstat] = thisdesc "\n\t(" isbetween[nstat] " pops, dependencies: " ( thisdep == "" ? "none" : thisdep ) ")"
}

function format_stats(   help, help_col1, total_width) {
  total_width=80
  for (i=1;i<=nstat;i++) {
    if (desc[i] ~ /^helper:/ && !("dependencies" in arg::args) ) { continue }
    help_col1[i] = name[i] 
    if ( (l=length(help_col1[i])) > help_col1_width ) { help_col1_width=l }
  }
  help_col1_width++
  remain_width = total_width - help_col1_width
  help="Available statistics:"
  for (i=1;i<=nstat;i++) {
    if (desc[i] ~ /^helper:/ && !("dependencies" in arg::args)) { continue }
    nlines=1
    if (length(desc[i]) > remain_width || desc[i]~/\n/) {
      for (j=remain_width+index(desc[i], "\n"); j<length(desc[i]); j+=remain_width+1) {
        for (k=0; k<remain_width/2; k++) {
          if (substr(desc[i], j-k, 1) ~ /[ ,)]/) { break }
        }
        if (k < remain_width/2-1) { j-=k }
        desc[i] = substr(desc[i], 1, j) "\n" substr(desc[i],j+1)
      }
      nlines=split(desc[i], desclines, "\n")
    }
    help = help "\n" sprintf("%-"help_col1_width"s", help_col1[i]) ( nlines > 1 ? desclines[1] : desc[i] )
    for (j=2;j<=nlines;j++) {
      help = help "\n" sprintf("%-"help_col1_width"s", " ") desclines[j] 
    }
  }
  return help
}
