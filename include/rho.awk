@namespace "calc"

BEGIN{
  stats::add_stat("rho", "Ronfort's rho", "between", "pi")
}

function initiate_rho(){
  if ( arg::args["mult"] == 1 ) {
    say("Warning: -s rho is not reliable in multiallelic sites")
  }
  for (i in groups) {
    if ( ploidy[i] % 1 != 0 ) {
      say("Warning: mixed ploidy population "g" will give unreliable results with -s rho" )
    }
  }
}

function increment_rho(i,j){
  # Here Hs = average of pi values of two populations,
  #      Ht = pi of two populations pooled,
  #      Hsp = Hs corrected for ploidy
  thisnum[i,j]["Hs"] = ( ( thisnum[i]["pi"]*thisden[j]["pi"] ) + ( thisnum[j]["pi"] * thisden[i]["pi"] ) )
  thisden[i,j]["Hs"] = 2 * thisden[i]["pi"] * thisden[j]["pi"] # same as thisden[i,j]["Hsp"]
  thisnum[i,j]["Hsp"] = ( ( thisnum[i]["pi"] * thisden[j]["pi"] * ( ploidy[i] - 1 ) / ploidy[i] ) + ( thisnum[j]["pi"] * thisden[i]["pi"] * (ploidy[j]-1) / ploidy[j] ) )
  thisden[i,j]["Hsp"] = thisden[i,j]["Hs"]
  thisnum[i,j]["Ht"]=(n[i] + n[j])^2
  for ( x in a_ij ) { 
    thisnum[i,j]["Ht"]-=(a[i][x]+a[j][x])^2 
  }
  thisden[i,j]["Ht"]=(n[i] + n[j]) * (n[i] + n[j] - 1)

  num[i,j]["Hs"]+=thisnum[i,j]["Hs"]
  den[i,j]["Hs"]+=thisden[i,j]["Hs"]
  num[i,j]["Hsp"]+=thisnum[i,j]["Hsp"]
  den[i,j]["Hsp"]+=thisden[i,j]["Hsp"]
  num[i,j]["Ht"]+=thisnum[i,j]["Ht"]
  den[i,j]["Ht"]+=thisden[i,j]["Ht"]
}

function finalize_rho(i,j,    Hs, Ht, Hsp, Hpt) {
  if ( den[i,j]["Hs"]==0 ) { return 0 }
  if ( den[i,j]["Hsp"]==0 ) { return 0 }
  if ( den[i,j]["Ht"]==0 ) { return 0 }
  Hs = num[i,j]["Hs"]/den[i,j]["Hs"]
  Ht = num[i,j]["Ht"]/den[i,j]["Ht"]
  Hsp= num[i,j]["Hsp"]/den[i,j]["Hsp"]
  Hpt= Hs + 2 * (Ht - Hs)
  num[i,j]["rho"]=Hpt-Hs
  den[i,j]["rho"]=Hpt-Hsp
}
