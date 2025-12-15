@namespace "calc"

BEGIN{
  stats::add_stat("nei", "Nei's D standard genetic distance", "between", "pi,dxy")
}

function finalize_nei(i,j,    jx, jy, jxy){
  if (den[i]["pi"]==0) { return 0 }
  if (den[j]["pi"]==0) { return 0 }
  if (den[i,j]["dxy"]==0) { return 0 }
  if (num[i]["pi"]==den[i]["pi"]) { return 0 }
  if (num[j]["pi"]==den[j]["pi"]) { return 0 }
  jx=1-num[i]["pi"]/den[i]["pi"]
  jy=1-num[j]["pi"]/den[j]["pi"]
  jxy=1-num[i,j]["dxy"]/den[i,j]["dxy"]
  num[i,j]["nei"]=-1 * log( jxy / sqrt( jx * jy ) )
  den[i,j]["nei"]=1
}
