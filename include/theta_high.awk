@namespace "calc"

BEGIN{
  stats::add_stat("theta_high", "Theta estimator based on sites with 0.33<=allele_freq<0.66", 0, "lines")
}

function initiate_theta_high(){
  say("Warning: -s theta_high concerns derived allele frequencies, make sure the VCF is polarized correctly!")
  piawka::assert( arg::args["persite"]!=1, "-s theta_high at a single site is meaningless" )
  piawka::assert( arg::args["mult"]!=1, "-s theta_high is not defined in multiallelic comparisons" )
}

function increment_theta_high(i,    alt, weight){
  if ( length(a[i]) > 1 && n[i] > 3 && 0 in a[i] ) {
    if ( a[i][0] <= n[i]/3 ) {
      weight=1/int(n[i]/3)
      alt=n[i]-a[i][0]
      thisnum[i]["theta_high"]= alt * weight
    }
  }
}

function finalize_theta_high(i) {
  den[i]["theta_high"]=num[i]["lines"]
}


