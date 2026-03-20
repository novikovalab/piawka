@namespace "calc"

BEGIN{
  stats::add_stat("afd", "average allele frequency difference", "between")
}

function initiate_afd(){
  if ( arg::args["mult"] == 1 ) {
    say("Warning: -s afd is not reliable in multiallelic sites")
  }
}

function increment_afd(i,j,    thisdiff){
  for ( x in a_ij ) { 
    thisdiff=( a[i][x]*n[j] - a[j][x]*n[i] ) / 2
    thisnum[i,j]["afd"]+=( thisdiff > 0 ? thisdiff : -thisdiff )
  }
  thisden[i,j]["afd"]=n[i]*n[j]
}

