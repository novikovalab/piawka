@namespace "calc"

BEGIN{
  stats::add_stat("fstwc", "fixation index, Weir & Cockerham's estimator", "between")
}

function increment_fstwc(i,j,    a1, a2, n1, n2, sizes, frac, mism, den){
  for ( x in a_ij ) { 
    a1=a[i][x]
    a2=a[j][x]
    n1=n[i]
    n2=n[j]
    # Formula from Bhatia et al. 2013, eq. (6)
    sizes = n1 * n2 / ( n1 + n2 )
    frac = 1 / ( n1 + n2 - 2 )
    mism = a1 * ( 1 - a1 / n1 ) + a2 * ( 1 - a2 / n2 )

    den = sizes * ( a1 / n1 - a2 / n2 )^2 + ( 2 * sizes - 1 ) * frac * mism
    thisnum[i,j]["fstwc"] += den - 2 * sizes * frac * mism
    thisden[i,j]["fstwc"] += den
    if ( arg::args["mult"] != 1 ) { break } # no need to increment twice for biallelic comparisons
  }
}
