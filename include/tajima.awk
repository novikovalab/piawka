@namespace "calc"

BEGIN{
  stats::add_stat("tajima", "Tajima's D", 0, "a1,a2,segr,pi,lines,miss")
}

function initiate_tajima(){
  piawka::assert( arg::args["persite"]!=1, "-s tajima cannot be calculated for a single site" )
  if ( arg::args["mult"]==1 ) {
    piawka::say("Warning: -s tajima is unreliable when calculated for multiallelic sites")
  }
  for (i in groups) {
    if ( ploidy[i] % 1 != 0 ) {
      say("Warning: mixed ploidy population "i" will give unreliable results with -s tajima" )
    }
  }
}

function finalize_tajima(i,     tP, a1, a2, N){
  if ( num[i]["segr"]==0 ) { return 0 }
  tP=num[i]["pi"]/den[i]["pi"]*num[i]["lines"]
  a1=num[i]["a1"]/num[i]["segr"]
  a2=num[i]["a2"]/num[i]["segr"]
  N=den[i]["miss"]/num[i]["lines"]
  num[i]["tajima"]= tP - num[i]["segr"]/a1
  den[i]["tajima"]= calcTajimaVar( num[i]["segr"], a1, a2, N )
}

# Tajima's D-like statistic calculation
# The formula below is classic Tajima's D but with missing data we calculate a1 and a2 in a slightly different way
function calcTajimaVar( S, a1, a2, n ) {
  if (n==3) { return 0 }
  if (n==1) { return 0 }
  if (n==0) { return 0 }
  if (a1==0) { return 0 }
  if (S==0) {return 0}
  b1=(n+1)/(3*n-3)
  b2=(2*(n^2+n+3))/(9*n*(n-1))
  c1=b1-1/a1
  c2=b2-(n+2)/(a1*n)+a2/(a1^2)
  e1=c1/a1
  e2=c2/(a1^2+a2)
  return sqrt(e1*S + e2*S*(S-1))
}
