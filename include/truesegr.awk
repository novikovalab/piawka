@namespace "calc"

BEGIN{
  stats::add_stat("truesegr", "helper: count of segregating sites", 0)
}

function increment_truesegr(i){
  if ( length(a[i])==2 ) { 
    thisnum[i]["truesegr"]++
    thisden[i]["truesegr"]=1 
  }
}
