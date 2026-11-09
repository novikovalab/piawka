@namespace "calc"

# helper: sum of the harmonic number a1(n) over all sites with a genotype count n>1, variable or not.
# E[segregating sites at a site] = theta * a1(n), so S / sum(a1) is the unbiased, minimum-variance estimate of theta per site
# when the sample size n differs between sites (missing data).
BEGIN{
  stats::add_stat("a1all", "helper: sum of 1st harmonic number over all sites", 0)
}

function increment_a1all(i){
  if ( n[i]>1 ) { thisnum[i]["a1all"]=harm(n[i]-1) }
}

function finalize_a1all(i){
  den[i]["a1all"]=1
}

