@namespace "calc"

BEGIN{
  stats::add_stat("theta_low", "Theta estimator based on sites with 0<allele_freq<0.33", 0, "lines")
}

function initiate_theta_low(){
  say("Warning: -s theta_low concerns derived allele frequencies, make sure the VCF is polarized correctly!")
  piawka::assert( arg::args["persite"]!=1, "-s theta_low at a single site is meaningless" )
  piawka::assert( arg::args["mult"]!=1, "-s theta_low is not defined in multiallelic comparisons" )
}

function increment_theta_low(i,    alt, weight){
  if ( length(a[i]) > 1 && n[i] > 3 && 0 in a[i] ) {
    if (a[i][0] > (2/3)*n[i]) {
      weight=1/int(n[i]/3)
      alt=n[i]-a[i][0]
      thisnum[i]["theta_low"]= alt * weight
    }
  }
}

function finalize_theta_low(i) {
  den[i]["theta_low"]=num[i]["lines"]
}
