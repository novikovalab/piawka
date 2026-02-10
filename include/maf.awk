@namespace "calc"

BEGIN{
  stats::add_stat("maf", "minor allele frequency")
}

function increment_maf(i,    maf){
  maf=n[i]
  for ( x in a[i] ) {
    if ( a[i][x] < maf ) {
      maf=a[i][x]
    }
  }
  thisnum[i]["maf"]=maf
  thisden[i]["maf"]=n[i]
}


