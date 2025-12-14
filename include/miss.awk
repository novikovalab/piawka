@namespace "calc"

BEGIN{
  stats::add_stat("miss", "share of missing genotype calls", 0, "", "sum")
}

function increment_miss(i){
  thisnum[i]["miss"]+=miss[i]
  thisden[i]["miss"]+=n[i]+miss[i]
}
