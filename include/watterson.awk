@namespace "calc"

BEGIN{
  stats::add_stat("watterson", "Watterson's theta", 0, "a1,truesegr,nlines")
}

function initiate_watterson(){
  piawka::assert( arg::args["persite"]!=1, "-s watterson at a single site is meaningless" )
}

function finalize_watterson(i){
  if ( num[i]["truesegr"]==0 ) { return 0 }
  num[i]["watterson"]=num[i]["truesegr"]/num[i]["nlines"]
  den[i]["watterson"]=num[i]["a1"]/num[i]["truesegr"]
}
