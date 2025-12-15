@namespace "calc"

BEGIN{
  stats::add_stat("theta_mid", "Theta estimator based on sites with 0.33<=allele_freq<0.66", 0, "lines")
}

function initiate_theta_mid(){
  say("Warning: -s theta_mid concerns derived allele frequencies, make sure the VCF is polarized correctly!")
  piawka::assert( arg::args["persite"]!=1, "-s theta_mid at a single site is meaningless" )
  piawka::assert( arg::args["mult"]!=1, "-s theta_mid is not defined in multiallelic comparisons" )
}

function increment_theta_mid(i,    alt, weight){
  if ( length(a[i]) > 1 && n[i] > 3 && 0 in a[i] ) {
    if ( a[i][0] > n[i]/3 && a[i][0] <= (2/3)*n[i] ) {
      weight=1/int(n[i]/3)
      alt=n[i]-a[i][0]
      thisnum[i]["theta_mid"]= alt * weight
    }
  }
}

function finalize_theta_mid(i) {
  den[i]["theta_mid"]=num[i]["lines"]
}

