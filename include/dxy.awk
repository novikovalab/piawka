@namespace "calc"

BEGIN{
  stats::add_stat("dxy", "absolute nucleotide divergence", "between")
}

function increment_dxy(i,j){
  # Add to dxy: probability that two a picked from two groups differ
  # subtraction rather than addition inspired by https://pubmed.ncbi.nlm.nih.gov/36031871/
  thisnum[i,j]["dxy"]=n[i]*n[j]
  thisden[i,j]["dxy"]=thisnum[i,j]["dxy"]
  for ( x in a_ij ) { 
    thisnum[i,j]["dxy"]-=a[i][x]*a[j][x]
  }
}
