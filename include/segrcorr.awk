@namespace "calc"

BEGIN{
  stats::add_stat("segrcorr", "helper: count of segregating sites corrected for missingness", 0, "", "sum")
}

function increment_segrcorr(i){
  if ( length(a[i])>1 ) { 
    thisnum[i]["segrcorr"]=recalcS_expected( 1, n[i], n[i]+miss[i] )
    thisden[i]["segrcorr"]=1 
  }
}

function recalcS_expected(S, n1, n2,    coef1, coef2, i ) {
  for (i=1; i<n1; i++) { coef1+=(1/i+1/(n1-i)) }
  for (i=1; i<n2; i++) { coef2+=(1/i+1/(n2-i)) }
  return S*(coef2/coef1)
}

function finalize_segrcorr(i){
  den[i]["segrcorr"]=1
}
