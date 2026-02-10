@namespace "calc"

BEGIN{
  stats::add_stat("daf", "derived allele frequency")
}

function initiate_daf(){
  if ( arg::args["mult"] == 1 ) {
    say("Warning: -s daf relies on REF alleles, make sure VCF is polarized correctly!")
  }
}

function increment_daf(i){
  thisnum[i]["daf"]=n[i]-a[i]["0"]
  thisden[i]["daf"]=n[i]
}

