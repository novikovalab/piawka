@namespace "calc"

BEGIN{
  stats::add_stat("theta_w", "Watterson's theta", 0, "a1,segr,lines")
}

function initiate_theta_w(){
  piawka::assert( arg::args["persite"]!=1, "-s theta_w at a single site is meaningless" )
}

function finalize_theta_w(i){
  if ( num[i]["segr"]==0 ) { return 0 }
  num[i]["theta_w"]=num[i]["segr"]/num[i]["lines"]
  den[i]["theta_w"]=num[i]["a1"]/num[i]["segr"]
}
