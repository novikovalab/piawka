@namespace "calc"

BEGIN{
  stats::add_stat("diff", "average allele frequency difference", "between")
}

function initiate_diff(){
  if ( arg::args["mult"] == 1 ) {
    say("Warning: -s diff is not reliable in multiallelic sites")
  }
}

function increment_diff(i,j,    thisdiff){
  for ( x in a_ij ) { 
    thisdiff=( a[i][x]*n[j] - a[j][x]*n[i] ) / 2
    thisnum[i,j]["diff"]+=( thisdiff > 0 ? thisdiff : -thisdiff )
  }
  thisden[i,j]["diff"]=n[i]*n[j]
}

