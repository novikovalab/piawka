@namespace "calc"

BEGIN{
  stats::add_stat("a1", "helper: 1st harmonic number", 0)
}

function increment_a1(i){
  if ( length(a[i])>1 ) { 
    thisnum[i]["a1"]+=harm(n[i]-1)
    thisden[i]["a1"]=1 
  }
}

function harm( n ) {
 # Return Nth harmonic number
 if ( n == 0 ) { return 0 }
 if ( !(n in _harm) ) { _harm[n]= 1/n + harm(n-1) }
 return _harm[n]
}
