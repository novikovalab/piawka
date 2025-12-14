@namespace "calc"

BEGIN{
  stats::add_stat("tajimalike", "Tajima's D-like statistic", 0, "a1,a2,segrcorr,pi,nlines")
}

function initiate_tajimalike(){
  piawka::assert( arg::args["persite"]!=1, "-s tajimalike cannot be calculated for a single site" )
  if ( arg::args["jobs"] > 1 ) {
    say("Warning: -s tajimalike is a bit less precise in multithreaded mode due to averaging across windows")
  }
  for (i in groups) {
    if ( ploidy[i] % 1 != 0 ) {
      say("Warning: mixed ploidy population "g" will give unreliable results with -s tajimalike" )
    }
  }
}

function finalize_tajimalike(i,     tP, a1, a2){
  if ( num[i]["segr"]==0 ) { return 0 }
  tP=num[i]["pi"]/den[i]["pi"]*num[i]["lines"]
  a1=num[i]["a1"]/num[i]["segr"]
  a2=num[i]["a2"]/num[i]["segr"]
  num[i]["tajimalike"]= tP - num[i]["segrcorr"]/a1
  den[i]["tajimalike"]=calcTajimaVar( num[i]["segrcorr"], a1, a2, n[i]+miss[i] )
}

# calcTajimaVar function defined in tajima.awk:
# # Tajima's D-like statistic calculation
# # The formula below is classic Tajima's D but with missing data we calculate a1 and a2 in a slightly different way
# function calcTajimaVar( S, a1, a2, n ) {
#   if (n==3) { return 0 }
#   if (n==1) { return 0 }
#   if (n==0) { return 0 }
#   if (a1==0) { return 0 }
#   if (S==0) {return 0}
#   b1=(n+1)/(3*n-3)
#   b2=(2*(n^2+n+3))/(9*n*(n-1))
#   c1=b1-1/a1
#   c2=b2-(n+2)/(a1*n)+a2/(a1^2)
#   e1=c1/a1
#   e2=c2/(a1^2+a2)
#   return sqrt(e1*S + e2*S*(S-1))
# }
