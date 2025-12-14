@namespace "calc"

function run(){ 
  check_gawk_version()
  check_htslib()
  help="\
    Calculate statistics within & between groups of samples in a VCF file.\n\
    Mandatory arguments are VCF file (--vcf) and grouping (--groups).\n\
    Default --stats are pi,dxy; other can be provided as a spaceless comma-seprarted list.\n\
    EXAMPLE: \n\tpiawka calc [OPTIONS] -v file.tsv -g  > distance_file"
  # add_argument args: shortopt, longopt, is_flag, description
  arg::add_argument("1", "persite", 1, "output values for each site")
  arg::add_argument("b", "bed", 0, "BED file with regions to be analyzed")
  arg::add_argument("B", "targets", 0, "BED file with targets (faster for numerous small regions)")
  arg::add_argument("g", "groups", 0, "either 2-columns sample / group table or \nkeywords \"unite\" (all samples in one group) or \"divide\" (each sample is a separate group)")
  arg::add_argument("j", "jobs", 0, "number of parallel jobs to run")
  arg::add_argument("l", "list", 1, "list all available statistics and exit")
  arg::add_argument("m", "mult", 1, "use populations with multiple alleles at a site")
  arg::add_argument("M", "miss", 0, "max share of missing GT per group at site, 0.0-1.0")
  arg::add_argument("q", "quiet", 1, "do not output progress and warning messages")
  arg::add_argument("R", "rand", 0, "randomly use this share of sites, 0.0-1.0")
  arg::add_argument("s", "stats", 0, "stats to calculate, comma-separated, e.g. \"pi,dxy,fst\"; full list under `piawka calc -l`")
  arg::add_argument("v", "vcf", 0, "gzipped and tabixed VCF file")
  arg::parse_args(2, help)

  SIGNAL_END_OF_BUFFER=SUBSEP SUBSEP SUBSEP
  piawka=ENVIRON["AWKPATH"]
  sub(/include:.*$/,"piawka",piawka)
  parse_arguments()
  print_header()
  configure_run()
  exit main()
}

function check_htslib(    bgzip_status, tabix_status) {
  bgzip_status = system("bgzip --version > /dev/null")
  tabix_status = system("tabix --version > /dev/null")
  if ( bgzip_status || tabix_status ) {
    print "Error: could not access dependencies: " ( bgzip_status ? "bgzip" : "" ) " " ( tabix_status ? "tabix" : "" )> "/dev/stderr"
    exit 1
  }
}
function check_gawk_version(    gawk_version) {
  gawk_version = substr(PROCINFO["version"],1,index(PROCINFO["version"],".")-1)
  if ( gawk_version < 5 ) {
    print "Error: GNU AWK v5.0.0 or above is needed to run piawka" > "/dev/stderr"
    exit 1
  }
}

function parse_arguments() {

  stats::add_stat("diff", "average allele frequency difference", "between")
  stats::add_stat("dxy", "absolute nucleotide divergence", "between")
  stats::add_stat("fst", "fixation index, Hudson's estimator", "between")
  stats::add_stat("fstwc", "fixation index, Weir & Cockerham's estimator", "between")
  stats::add_stat("pi", "expected heterozygosity = nucleotide diversity")
  stats::add_stat("tajima", "Tajima's D", "a1,a2,truesegr,pi")
  stats::add_stat("tajimalike", "Tajima's D-like statistic", 0, "a1,a2,segr,pi")
  stats::add_stat("rho", "Ronfort's rho", "between", "pi")
  stats::add_stat("watterson", "Watterson's theta", 0, "a1,truesegr")
  stats::add_stat("lines", "number of lines used in calculation", 0)
  stats::add_stat("miss", "share of missing genotype calls", 0)
  # Helper statistics: accumulators for other stats, are not shown with piawka -l
  stats::add_stat("a1", "helper: 1st harmonic number", 0)
  stats::add_stat("a2", "helper: 2nd harmonic number", 0)
  stats::add_stat("truesegr", "helper: count of segregating sites", 0)
  stats::add_stat("segr", "helper: count of segregating sites adjusted for missingness", 0)

  if ( arg::args["list"] ) {
    print stats::format_stats()
    exit 0
  }
  if ( arg::args["stats"]=="" ) {
    arg::args["stats"]="pi,dxy"
  }

  stats::parse_stats(arg::args["stats"])
  any_within="within" in stats::stats
  any_between="between" in stats::stats

  # Some arg checks
  piawka::assert( arg::args["vcf"] != "", "required argument: -v <file.vcf.gz>" )
  piawka::assert( arg::args["groups"] != "", "required argument: -g <groups.tsv>" )
  if ( arg::args["bed"] != "" ) { piawka::check_file( arg::args["bed"] ) }
  if ( arg::args["targets"] != "" ) { piawka::check_file( arg::args["targets"] ) }
}

function configure_run() {
# Set default variable values, prepare index arrays and pipes

  if ( arg::args["miss"] == "" ) { arg::args["miss"]=1 }

  if ( ( arg::args["groups"] != "unite" ) && ( arg::args["groups"] != "divide" ) ) {
    get_groups()
  }

  get_header()

  for (i in stats::statargs) {
    s=stats::statargs[i]
    init="calc::initiate_"s
    if ( init in FUNCTAB ) {
      @init()
    }
  }

  tmpcmd="mktemp -d piawkatmp.XXXXXXXX"
  if ( tmpcmd | getline tmpdir <= 0 ) {
    say("Error: failed to create temporary directory!")
    exit 1
  }
  close(tmpcmd)
  # We can't trap an exec'ed script but we can employ a process to clean the tmpdir after it's done:
  system( "{ while kill -0 "PROCINFO["pid"]" 2>/dev/null; do sleep 1; done; rm -r "tmpdir"; } &" )
  if ( arg::args["jobs"]==0 ) {
    arg::args["jobs"]=1 
  }
}

function main() {

  # Children: listen to query regions dispenser
  for ( jobnum=0; jobnum < arg::args["jobs"]; jobnum++ ) { 
    if ( arg::args["rand"] ) { srand( awk::xor( systime(), PROCINFO["pid"] ) ) } # processes have different random seeds
    buffer = "cat #"jobnum
    printf "" |& buffer # initialize pipe
    pid=awk::fork()
    if ( pid>0 ) { continue } 
    close(buffer, "to") # children don't print to buffer
    bufl=0
    while ( $0 != SIGNAL_END_OF_BUFFER ) {
      curr= jobnum + arg::args["jobs"] * (++bufl - 1)
      tmpf=tmpdir "/" curr ".tmp"
      getl=(buffer |& getline)
      if (getl <= 0) {
        print "Error: job "jobnum" did not reach end of buffer:" > "/dev/stderr"
        print $0 > "/dev/stderr"
        return getl
      }
      chr=$1
      start=$2
      end=$3
      locus=$4
      process_sites()
      # in case tmpf is empty, add a few bytes
      print SIGNAL_END_OF_BUFFER > tmpf
      close(tmpf)
    }
    close(buffer)
    return 0
  }
  
  # Parent process: dispense query regions or just scan the whole file
  bedcmd = get_bedcmd()
  bufl=0
  while ( bedcmd | getline > 0 ) {
    bufl++
    jobnum=bedline++ % arg::args["jobs"]
    buffer="cat #"jobnum # number of pipes with regions == # jobs
    print $0 |& buffer
  }
  for ( jobnum=0; jobnum < arg::args["jobs"]; jobnum++ ) { 
    buffer="cat #"jobnum 
    print SIGNAL_END_OF_BUFFER |& buffer
    close(i, "to")
  }
  close( bedcmd )
  return 0
}

function get_bedcmd() {
  if ( arg::args["bed"] != "" ) {
    return "cat "arg::args["bed"]
  } else if (arg::args["targets"] != "") {
    return piawka" win -T "arg::args["targets"] 
  } else {
    return piawka" win -v "arg::args["vcf"]" -s "header_length
  }
}

# Groups file (first in command line): store lists of group members in `groups` array
function get_groups() { 
  piawka::check_file( arg::args["groups"] )
  while ( getline < arg::args["groups"] > 0 ) {
    if (!checked_groups) { piawka::assert(NF==2, "the groups file must contain two columns"); checked_groups=1 }
      groupmem[$1]=$2
  }
  close( arg::args["groups"] )
}

# VCF header: assign samples to groups
function get_header() {
  cmd="tabix -H " arg::args["vcf"]
  while ( cmd | getline > 0 ) { 
    header_length++
    if ($0 ~ /^#CHROM/ ) {

      # Assign sample positions to groups
      for (i=10; i<=NF; i++) {
        if ( arg::args["groups"]=="unite" ) {
          groupindex[i]="all_samples"
          groups["all_samples"]++
        } else if ( arg::args["groups"]=="divide" ) {
          groupindex[i]=$i
          groups[$i]++
        } else if ( groupmem[$i] != "" ) { 
          groupindex[i]=groupmem[$i]
          groups[groupmem[$i]]++
        }
      }
      # Assign ploidy to groups using variant calls from next line, detect mixed-ploidy groups
      cmd | getline
      # Count alleles for each group, then divide by number of samples
      for (i in groupindex) {
        ploidy[groupindex[i]] += ( ( p = index($i, ":") ) > 0 ? p : length($i)+1 )
      }
      for (g in groups) {
        ploidy[g] /= 2*groups[g]
      }
      break
    }
  }
  close(cmd)
}

function initiate_diff(){
  if ( arg::args["mult"] == 1 ) {
    say("Warning: -s diff is not reliable in multiallelic sites")
  }
}

function initiate_tajima(){
  piawka::assert( arg::args["persite"]!=1, "-s tajima cannot be calculated for a single site" )
  if ( arg::args["jobs"] > 1 ) {
    say("Warning: -s tajima is a bit less precise in multithreaded mode due to averaging across windows")
  }
  for (g in groups) {
    if ( ploidy[g] % 1 != 0 ) {
      say("Warning: mixed ploidy population "g" will give unreliable results with -s tajima" )
    }
  }
}

function initiate_watterson(){
  piawka::assert( arg::args["persite"]!=1, "-s watterson at a single site is meaningless" )
}

function initiate_tajimalike(){
  piawka::assert( arg::args["persite"]!=1, "-s tajimalike cannot be calculated for a single site" )
  if ( arg::args["jobs"] > 1 ) {
    say("Warning: -s tajimalike is a bit less precise in multithreaded mode due to averaging across windows")
  }
  for (g in groups) {
    if ( ploidy[g] % 1 != 0 ) {
      say("Warning: mixed ploidy population "g" will give unreliable results with -s tajimalike" )
    }
  }
}

function initiate_rho(){
  if ( arg::args["mult"] == 1 ) {
    say("Warning: -s rho is not reliable in multiallelic sites")
  }
  for (g in groups) {
    if ( ploidy[g] % 1 != 0 ) {
      say("Warning: mixed ploidy population "g" will give unreliable results with -s rho" )
    }
  }
}

# Process VCF lines
function process_sites() {
  if ( locus == "" ) { locus="_" } # empty locus breaks piawka sum
  cmd = "tabix " arg::args["vcf"] ( arg::args["targets"]!="" ? " -T " arg::args["targets"]:"" ) " "chr":"start+1"-"end  
  while ( cmd | getline > 0 ) {

    if ( arg::args["rand"]!="" && ( rand() > arg::args["rand"] ) ) { continue }

    # Process only SNPs (possibly monomorphic or multiallelic)
    # To obtain results identical to ksamuk/pixy, set $4 !~ /^[ACGT]$/ && $5 !~ /\*|,|[ACGT][ACGT]/
    if (length($4)>1) { continue } # if REF is multi-nucleotide, skip
    delete ALT
    ALT[a=0]=1

    # only consider single-char ALT that are not *
    while (x=index($5,",")) {
      a++
      if (x == 2 && substr($5,1,1)!="*") {
        ALT[a]=1
      }
      $5=substr($5,x+1)
    }
    if ( length($5) == 1 && $5 != "*" ) {
      ALT[++a]=1
    }

    if ( arg::args["persite"] == 1 ) { 
      start=$2-1
      end=$2
    }
  
    # Reset site-specific parameters
    for (g in groups) { miss[g]=0; misind[g]=0; nalleles[g]=0; exclude[g]=0 }
    delete alleles
    delete thisnum
    delete thisden

    # Pool GT values for groups, count each state and missing data
    for (i in groupindex) {
      g=groupindex[i]
      if ( exclude[g] ) { continue }
      gtend=index( $i, ":" )
      if ( gtend==0 ) { gtend=length($i)+1 }
      for (c=1; c<gtend; c+=2) {
        al=substr($i,c,1)
        if ( al == "." ) {
          miss[g]+=gtend/2
          break
        } else {
          if ( !(al in ALT) ) { 
            exclude[g]=1
            break
          }
          alleles[g][al]++
          nalleles[g]++
        }
      }
      if ( arg::args["groups"] == "divide" && g in alleles ) {
        calculate_within(g)
        for (g2 in alleles) {
          if ( g2 > g ) {
            calculate_between(g,g2)
          } else if ( g2 < g ) {
            calculate_between(g2,g)
          }
        }
        if ( arg::args["persite"] ) { yield_output() }
      }
    }

    if ( arg::args["groups"] == "divide" ) { continue }
    for ( g in alleles ) {
      if ( exclude[g] || ( arg::args["miss"] < 1 && miss[g]/(miss[g]+nalleles[g]) > arg::args["miss"] ) || ( arg::args["mult"] != 1 && length(alleles[g]) > 2 ) ) { 
        delete alleles[g]
        continue
      }
      calculate_within( g )
    }
    for ( g in alleles ) {
      for ( g2 in alleles ) {
        # Only making comparisons one way
        if (g2>g) { 
          calculate_between(g,g2)
        }
      }
    }
    if ( arg::args["persite"] ) { yield_output() }
  }
  if ( !arg::args["persite"] ) { yield_output() }
  close(cmd)
}

function yield_output() {
  if ( pid>0 ) { return 0 }
  for (i in den) {
    if ( split(i, ii, SUBSEP) == 1 ) {
      for ( s in stats::stats["within"] ) {
        fin="calc::finalize_"s
        if ( fin in FUNCTAB ) {
          @fin(i)
        }
        if (s in stats::stats_print["within"]) {
          printOutput( i, "", s ) 
        }
      }
    } else {
      for ( s in stats::stats["between"] ) {
        fin="calc::finalize_"s
        if ( fin in FUNCTAB ) {
          @fin(ii[1],ii[2])
        }
        if (s in stats::stats_print["between"]) {
          printOutput( ii[1], ii[2], s )
        }
      }
    }
  }
  # Reset locus-specific parameters
  delete num
  delete den
}

function finalize_watterson(g){
  if ( num[g]["truesegr"]==0 ) { return 0 }
  num[g]["watterson"]=num[g]["truesegr"]/nUsed[g]
  den[g]["watterson"]=num[g]["a1"]/num[g]["truesegr"]
}

function finalize_tajima(g,     tP, a1, a2){
  if ( num[g]["truesegr"]==0 ) { return 0 }
  tP=num[g]["pi"]/den[g]["pi"]*nUsed[g]
  a1=num[g]["a1"]/num[g]["truesegr"]
  a2=num[g]["a2"]/num[g]["truesegr"]
  num[g]["tajima"]= tP - num[g]["truesegr"]/a1
  den[g]["tajima"]= calcTajimaVar( num[g]["truesegr"], a1, a2, nalleles[g]+miss[g] )
}

function finalize_tajimalike(g,     tP, a1, a2){
  if ( num[g]["truesegr"]==0 ) { return 0 }
  tP=num[g]["pi"]/den[g]["pi"]*nUsed[g]
  a1=num[g]["a1"]/num[g]["truesegr"]
  a2=num[g]["a2"]/num[g]["truesegr"]
  num[g]["tajimalike"]= tP - num[g]["segr"]/a1
  den[g]["tajimalike"]=calcTajimaVar( num[g]["segr"], a1, a2, nalleles[g]+miss[g] )
}

END {
  if (pid>0) {
    say("Total jobs to run: "bufl + (bufl % arg::args["jobs"] > 0 ? arg::args["jobs"] : 0)-1"\n", 1 )
    for (f=0; f<bufl; f++) {
      tmpf=tmpdir"/"f".tmp"
      # to avoid race condition given persite, wait until the child writes to next region
      tmpfn=tmpdir"/"f+arg::args["jobs"]".tmp"
      while ( awk::stat(tmpfn, _) < 0 || _["size"]==0 ) {
        system("sleep 1")
      }
      say("Finalizing job " f+arg::args["jobs"], 1)
      while ( getline < tmpf > 0 ) {
        if ( $0 == SIGNAL_END_OF_BUFFER ) {
          break
        }
        if ( arg::args["bed"] == "" && arg::args["persite"] != 1 ) {
          print | piawka" sum -"
        } else {
          print
        }
      }
      close(tmpf)
      system("rm -f "tmpf)
    }
    # no need to wait since all children confirmedly written output?..
    close(piawka" sum -")
  }
}

function increment_pi(g){
  # Add to pi: probability that two randomly picked alleles differ
  # New formula below from https://pubmed.ncbi.nlm.nih.gov/36031871/
  thisnum[g]["pi"]=nalleles[g]^2
  thisden[g]["pi"]=nalleles[g]*(nalleles[g]-1)
  for ( x in alleles[g] ) { thisnum[g]["pi"]-=alleles[g][x]^2 }
}

function increment_truesegr(g){
  if ( length(alleles[g])==2 ) { 
    thisnum[g]["truesegr"]++
    thisden[g]["truesegr"]=1 
  }
}

function increment_segr(g){
  if ( length(alleles[g])==2 ) { 
    thisnum[g]["segr"]+=recalcS_expected( 1, nalleles[g], nalleles[g]+miss[g] )
    thisden[g]["segr"]=1 
  }
}

function increment_a1(g){
  if ( length(alleles[g])==2 ) { 
    thisnum[g]["a1"]+=harm(nalleles[g]-1)
    thisden[g]["a1"]=1 
  }
}
   
function increment_a2(g){
  if ( length(alleles[g])==2 ) { 
    thisnum[g]["a2"]+=harm2(nalleles[g]-1)
    thisden[g]["a2"]=1 
  }
}

function increment_lines(g){
  thisnum[g]["lines"]++
}

function finalize_lines(g){
  den[g]["lines"]=1
}

function increment_miss(g){
  thisnum[g]["miss"]+=miss[g]
  thisden[g]["miss"]+=nalleles[g]+miss[g]
}

function calculate_within(g) {
  if (!any_within) { return 0 }
  for ( s in stats::stats["within"] ) {
    incr="calc::increment_"s
    if ( incr in FUNCTAB ) {
      @incr(g)
    }
    num[g][s]+=thisnum[g][s]
    den[g][s]+=thisden[g][s] 
  }
}

function increment_dxy(g,g2){
  # Add to dxy: probability that two alleles picked from two groups differ
  # subtraction rather than addition inspired by https://pubmed.ncbi.nlm.nih.gov/36031871/
  thisnum[g,g2]["dxy"]=nalleles[g]*nalleles[g2]
  thisden[g,g2]["dxy"]=thisnum[g,g2]["dxy"]
  for ( x in bothalleles ) { 
    thisnum[g,g2]["dxy"]-=alleles[g][x]*alleles[g2][x]
  }
}

function increment_diff(g,g2,    thisdiff){
  for ( x in bothalleles ) { 
    thisdiff=( alleles[g][x]*nalleles[g2] - alleles[g2][x]*nalleles[g] ) / 2
    thisnum[g,g2]["diff"]+=( thisdiff > 0 ? thisdiff : -thisdiff )
  }
  thisden[g,g2]["diff"]=nalleles[g]*nalleles[g2]
}

function increment_fst(g,g2,    a1, a2, n1, n2, hw, hb){
  for ( x in bothalleles ) { 
    a1=alleles[g][x]
    a2=alleles[g2][x]
    n1=nalleles[g]
    n2=nalleles[g2]
    pi1 = a1 * (n1 - a1) / (n1*(n1-1))
    pi2 = a2 * (n2 - a2) / (n2*(n2-1))
    hw = pi1 + pi2
    hb = ( a1 * (n2-a2) + a2*(n1-a1) ) / (n1*n2)
    thisnum[g,g2]["fst"]+=hb-hw
    thisden[g,g2]["fst"]+=hb
    if ( arg::args["mult"] != 1 ) { break } # no need to increment twice for biallelic comparisons
  }
} 

function increment_fstwc(g,g2,    a1, a2, n1, n2, sizes, frac, mism, den){
  for ( x in bothalleles ) { 
    a1=alleles[g][x]
    a2=alleles[g2][x]
    n1=nalleles[g]
    n2=nalleles[g2]
    # Formula from Bhatia et al. 2013, eq. (6)
    sizes = n1 * n2 / ( n1 + n2 )
    frac = 1 / ( n1 + n2 - 2 )
    mism = a1 * ( 1 - a1 / n1 ) + a2 * ( 1 - a2 / n2 )

    den = sizes * ( a1 / n1 - a2 / n2 )^2 + ( 2 * sizes - 1 ) * frac * mism
    thisnum[g,g2]["fstwc"] += den - 2 * sizes * frac * mism
    thisden[g,g2]["fstwc"] += den
    if ( arg::args["mult"] != 1 ) { break } # no need to increment twice for biallelic comparisons
  }
}

function increment_rho(g,g2){
  # Here Hs = average of pi values of two populations,
  #      Ht = pi of two populations pooled,
  #      Hsp = Hs corrected for ploidy
  thisnum[g,g2]["Hs"] = ( ( thisnum[g]["pi"]*thisden[g2]["pi"] ) + ( thisnum[g2]["pi"] * thisden[g]["pi"] ) )
  thisden[g,g2]["Hs"] = 2 * thisden[g]["pi"] * thisden[g2]["pi"] # same as thisden[g,g2]["Hsp"]
  thisnum[g,g2]["Hsp"] = ( ( thisnum[g]["pi"] * thisden[g2]["pi"] * ( ploidy[g] - 1 ) / ploidy[g] ) + ( thisnum[g2]["pi"] * thisden[g]["pi"] * (ploidy[g2]-1) / ploidy[g2] ) )
  thisden[g,g2]["Hsp"] = thisden[g,g2]["Hs"]
  thisnum[g,g2]["Ht"]=(nalleles[g] + nalleles[g2])^2
  for ( x in bothalleles ) { 
    thisnum[g,g2]["Ht"]-=(alleles[g][x]+alleles[g2][x])^2 
  }
  thisden[g,g2]["Ht"]=(nalleles[g] + nalleles[g2]) * (nalleles[g] + nalleles[g2] - 1)

  num[g,g2]["Hs"]+=thisnum[g,g2]["Hs"]
  den[g,g2]["Hs"]+=thisden[g,g2]["Hs"]
  num[g,g2]["Hsp"]+=thisnum[g,g2]["Hsp"]
  den[g,g2]["Hsp"]+=thisden[g,g2]["Hsp"]
  num[g,g2]["Ht"]+=thisnum[g,g2]["Ht"]
  den[g,g2]["Ht"]+=thisden[g,g2]["Ht"]
}

function calculate_between(g,g2) {
  if (!any_between) { return 0 }
  # Extract alleles of the group
  delete bothalleles
  for (xx in alleles[g]) {
    bothalleles[xx]++
    # gawk breaks if you don't set uninitialized array elems to 0
    if ( !(xx in alleles[g2]) ) { alleles[g2][xx]=0 }
  }
  for (yy in alleles[g2]) {
    if ( !(yy in bothalleles) ) { 
      bothalleles[yy]++ # so far keeps alleles from comparisons of g with previous groups
      alleles[g][yy]=0
    }
  }
  # If not arg::args["mult"] == 1, proceed only if common allele pool has <=2 alleles
  if ( arg::args["mult"] != 1 && length(bothalleles) > 2 ) { return 1 }
  for ( s in stats::stats["between"] ) {
    incr="calc::increment_"s
    if ( incr in FUNCTAB ) {
      @incr(g,g2)
    }
    num[g,g2][s]+=thisnum[g,g2][s]
    den[g,g2][s]+=thisden[g,g2][s] 
  }
}

function finalize_rho(g,g2,    Hs, Ht, Hsp, Hpt) {
  if ( den[g,g2]["Hs"]==0 ) { return 0 }
  if ( den[g,g2]["Hsp"]==0 ) { return 0 }
  if ( den[g,g2]["Ht"]==0 ) { return 0 }
  Hs = num[g,g2]["Hs"]/den[g,g2]["Hs"]
  Ht = num[g,g2]["Ht"]/den[g,g2]["Ht"]
  Hsp= num[g,g2]["Hsp"]/den[g,g2]["Hsp"]
  Hpt= Hs + 2 * (Ht - Hs)
  num[g,g2]["rho"]=Hpt-Hs
  den[g,g2]["rho"]=Hpt-Hsp
}

function printOutput( pop1, pop2, metric,    idx ) {
  if (pop2=="") {
    idx=pop1
    pop2="."
  } else {
    idx=pop1 SUBSEP pop2 
  }
  if ( den[idx][metric]==0 ) { return 0 }
  out=chr"\t"start"\t"end"\t"locus"\t"pop1"\t"pop2"\t"metric"\t"num[idx][metric]/den[idx][metric]"\t"num[idx][metric]"\t"den[idx][metric]
  print out > tmpf
}

function harm( n ) {
 # Return Nth harmonic number
 if ( n == 0 ) { return 0 }
 if ( !(n in _harm) ) { _harm[n]= 1/n + harm(n-1) }
 return _harm[n]
}

function harm2( n ) {
 # Return Nth second harmonic number
 if ( n == 0 ) { return 0 }
 if ( !(n in _harm2) ) { _harm2[n]= 1/n^2 + harm2(n-1) }
 return _harm2[n]
}

# Tajima's D-like statistic calculation
# The formula below is classic Tajima's D but with missing data we calculate a1 and a2 in a slightly different way
function calcTajimaVar( S, a1, a2, n ) {
  if (n==3) { return 0 }
  if (n==1) { return 0 }
  if (n==0) { return 0 }
  if (a1==0) { return 0 }
  if (S==0) {return 0}
  b1=(n+1)/(3*n-3)
  b2=(2*(n^2+n+3))/(9*n*(n-1))
  c1=b1-1/a1
  c2=b2-(n+2)/(a1*n)+a2/(a1^2)
  e1=c1/a1
  e2=c2/(a1^2+a2)
  return sqrt(e1*S + e2*S*(S-1))
}

function recalcS_expected(S, n1, n2,    coef1, coef2 ) {
  for (i=1; i<n1; i++) { coef1+=(1/i+1/(n1-i)) }
  for (i=1; i<n2; i++) { coef2+=(1/i+1/(n2-i)) }
  return S*(coef2/coef1)
}

function print_header() {
    print "#chr", "start", "end", "locus", "pop1", "pop2", "metric", "value", "numerator", "denominator"
}

function say(string, same_line) {
  if (arg::args["quiet"]) { return }
  if (same_line) {
    printf "\033[2K\r" string > "/dev/stderr"
  } else {
    print string > "/dev/stderr"
  }
}
