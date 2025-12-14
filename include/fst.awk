@namespace "calc"

BEGIN{
  stats::add_stat("fst", "fixation index, Hudson's estimator", "between")
}

function increment_fst(i,j,    a1, a2, n1, n2, hw, hb){
  for ( x in a_ij ) { 
    a1=a[i][x]
    a2=a[j][x]
    n1=n[i]
    n2=n[j]
    pi1 = a1 * (n1 - a1) / (n1*(n1-1))
    pi2 = a2 * (n2 - a2) / (n2*(n2-1))
    hw = pi1 + pi2
    hb = ( a1 * (n2-a2) + a2*(n1-a1) ) / (n1*n2)
    thisnum[i,j]["fst"]+=hb-hw
    thisden[i,j]["fst"]+=hb
    if ( arg::args["mult"] != 1 ) { break } # no need to increment twice for biallelic comparisons
  }
} 
