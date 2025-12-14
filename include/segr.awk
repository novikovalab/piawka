@namespace "calc"

BEGIN{
  stats::add_stat("segr", "helper: count of segregating sites adjusted for missingness", 0)
}

function increment_segr(i){
  if ( length(a[i])==2 ) { 
    thisnum[i]["segr"]+=recalcS_expected( 1, n[i], n[i]+miss[i] )
    thisden[i]["segr"]=1 
  }
}

function recalcS_expected(S, n1, n2,    coef1, coef2, i ) {
  for (i=1; i<n1; i++) { coef1+=(1/i+1/(n1-i)) }
  for (i=1; i<n2; i++) { coef2+=(1/i+1/(n2-i)) }
  return S*(coef2/coef1)
}
