@namespace "calc"

BEGIN{
  stats::add_stat("segr", "helper: count of segregating sites", 0, "")
}

function increment_segr(i){
  if ( length(a[i])>1 ) { 
    thisnum[i]["segr"]=1
  }
}

function finalize_segr(i){
  den[i]["segr"]=1
}
