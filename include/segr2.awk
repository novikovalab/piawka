@namespace "calc"

BEGIN{
  stats::add_stat("segr2", "# sites segregating in pop1 + pop2 pool", "between", "")
}

function increment_segr2(i,j){
  if ( length(a_ij)>1 ) {
    thisnum[i,j]["segr2"]=1
  }
}

function finalize_segr2(i,j){
  den[i,j]["segr2"]=1
}
