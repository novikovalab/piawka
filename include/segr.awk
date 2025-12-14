@namespace "calc"

BEGIN{
  stats::add_stat("segr", "helper: count of segregating sites", 0, "", "sum")
}

function increment_segr(i){
  if ( length(a[i])>1 ) { 
    thisnum[i]["segr"]++
    thisden[i]["segr"]=1 
  }
}
