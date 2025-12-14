@namespace "calc"

BEGIN{
  stats::add_stat("pi", "expected heterozygosity = nucleotide diversity")
}

function increment_pi(i){
  # Add to pi: probability that two randomly picked a differ
  # New formula below from https://pubmed.ncbi.nlm.nih.gov/36031871/
  thisnum[i]["pi"]=n[i]^2
  thisden[i]["pi"]=n[i]*(n[i]-1)
  for ( x in a[i] ) { thisnum[i]["pi"]-=a[i][x]^2 }
}
