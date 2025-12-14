@namespace "calc"

BEGIN{
  stats::add_stat("a2", "helper: 2nd harmonic number", 0)
}

function increment_a2(i){
  if ( length(a[i])==2 ) { 
    thisnum[i]["a2"]+=harm2(n[i]-1)
    thisden[i]["a2"]=1 
  }
}

function harm2( n ) {
 # Return Nth second harmonic number
 if ( n == 0 ) { return 0 }
 if ( !(n in _harm2) ) { _harm2[n]= 1/n^2 + harm2(n-1) }
 return _harm2[n]
}
