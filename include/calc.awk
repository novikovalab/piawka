@namespace "calc"

function run(){ 
  time_start=awk::systime()
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
  arg::add_argument("m", "mult", 1, "use populations with multiple a at a site")
  arg::add_argument("M", "miss", 0, "max share of missing GT per group at site, 0.0-1.0")
  arg::add_argument("q", "quiet", 1, "do not output progress and warning messages")
  arg::add_argument("R", "rand", 0, "randomly use this share of sites, 0.0-1.0")
  arg::add_argument("s", "stats", 0, "stats to calculate, comma-separated, e.g. \"pi,dxy,fst\"; full list under `piawka calc -l`")
  arg::add_argument("v", "vcf", 0, "gzipped and tabixed VCF file")
  arg::parse_args(2, help)

  SIGNAL_END_OF_BUFFER=SUBSEP SUBSEP SUBSEP
  piawka=ENVIRON["AWKPATH"]
  sub(/include:.*$/,"piawka",piawka)
  check_arguments()
  print_header()
  make_tmpdir()
  exit main()
}

function print_header() {
  # has to be before make_tmpdir otherwise children print own headers
  print "#chr\tstart\tend\tlocus\tpop1\tpop2\tmetric\tvalue\tnumerator\tdenominator"
}

function check_htslib(    bgzip_status, tabix_status) {
  bgzip_status = system("bgzip --version > /dev/null")
  tabix_status = system("tabix --version > /dev/null")
  if ( bgzip_status || tabix_status ) {
    say("Error: could not access dependencies: " ( bgzip_status ? "bgzip" : "" ) " " ( tabix_status ? "tabix" : "" ) )
    exit 1
  }
}
function check_gawk_version(    gawk_version) {
  gawk_version = substr(PROCINFO["version"],1,index(PROCINFO["version"],".")-1)
  if ( gawk_version < 5 ) {
    say("Error: GNU AWK v5.0.0 or above is needed to run piawka")
    exit 1
  }
}

function check_arguments() {
  if ( arg::args["list"] ) {
    print stats::format_stats()
    exit 0
  }
  if ( arg::args["stats"]=="" ) { arg::args["stats"]="pi,dxy" }
  if ( arg::args["miss"] == "" ) { arg::args["miss"]=1 }
  if ( arg::args["jobs"]==0 ) { arg::args["jobs"]=1 }
  stats::parse_stats(arg::args["stats"])
  any_within="within" in stats::stats
  any_between="between" in stats::stats
  piping_to_sum=( arg::args["bed"] == "" && arg::args["persite"] != 1  )
  # if piping to sum, dependencies should be also passed over
  if ( piping_to_sum ) {
    piawka::copy_array(stats::stats, stats::stats_print)
  }
  # Some arg checks
  piawka::assert( arg::args["vcf"] != "", "required argument: -v <file.vcf.gz>" )
  piawka::assert( arg::args["groups"] != "", "required argument: -g <groups.tsv>" )
  if ( arg::args["bed"] != "" ) { piawka::check_file( arg::args["bed"] ) }
  if ( arg::args["targets"] != "" ) { piawka::check_file( arg::args["targets"] ) }

  divide=arg::args["groups"]=="divide"
  if ( ( arg::args["groups"] != "unite" ) && !divide ) {
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
}

# Groups file (first in command line): store lists of group members in `groups` array
function get_groups() { 
  piawka::check_file( arg::args["groups"] )
  while ( getline < arg::args["groups"] > 0 ) {
    if (!seen_groupsfile) { 
      seen_groupsfile=1 
      piawka::assert(NF==2, "the groups file must contain two columns")
    }
    groupmem[$1]=$2
  }
  close( arg::args["groups"] )
}

# VCF header: assign samples to groups
function get_header() {
  cmd="tabix -H " arg::args["vcf"]
  while ( cmd | getline > 0 ) { 
    if ($0 ~ /^#CHROM/ ) {

      # Assign sample positions to groups
      for (col=10; col<=NF; col++) {
        if ( arg::args["groups"]=="unite" ) {
          groupindex[col]="all_samples"
          groups["all_samples"]++
        } else if ( divide ) {
          groupindex[col]=$col
          groups[$col]++
          for ( j in groups ) {
            if ( j != $col ) { pairs[$col][j]=1 }
          }
        } else if ( $col in groupmem ) { 
          groupindex[col]=groupmem[$col]
          groups[groupmem[$col]]++
        }
      }
      if ( seen_groupsfile ) {
        divide=1
        for ( i in groups ) {
          if ( divide && groups[i]>1 ) { divide=0 }
          for ( j in groups ) {
            if ( j > i ) {
              pairs[i][j]=1
            }
          }
        }
      }
      # Assign ploidy to groups using variant calls from next line, detect mixed-ploidy groups
      cmd | getline
      # Count a for each group, then divide by number of samples
      for (col in groupindex) {
        ploidy[groupindex[col]] += ( ( p = index($col, ":") ) > 0 ? p : length($col)+1 )
      }
      for (i in groups) {
        ploidy[i] /= 2*groups[i]
      }
      break
    }
  }
  close(cmd)
}

function make_tmpdir(){
  tmpcmd="mktemp -d piawkatmp.XXXXXXXX"
  if ( tmpcmd | getline tmpdir <= 0 ) {
    say("Error: failed to create temporary directory!")
    exit 1
  }
  close(tmpcmd)
  # We can't trap an exec'ed script but we can employ a process to clean the tmpdir after it's done:
  system( "{ while kill -0 "PROCINFO["pid"]" 2>/dev/null; do sleep 1; done; rm -r "tmpdir"; } &" )
}

function main() {
 
  # Children: listen to query regions dispenser
  for ( jobnum=0; jobnum < arg::args["jobs"]; jobnum++ ) { 
    if ( arg::args["rand"] ) { srand( awk::xor( awk::systime(), PROCINFO["pid"] ) ) } # processes have different random seeds
    buffer=tmpdir"/buffer_"jobnum".tmp" # number of pipes with regions == # jobs
    pid=awk::fork()
    if ( pid>0 ) { continue } 
    bufl=0
    bytes=0
    while ( $0 != SIGNAL_END_OF_BUFFER ) {
      while ( _["size"] <= bytes ) { awk::stat(buffer, _) } #wait till file gets bigger
      curr= jobnum + arg::args["jobs"] * (++bufl - 1)
      tmpf=tmpdir "/" curr ".tmp"
      piawka::assert((getl=getline < buffer)>0, "could not reach end of buffer "buffer": "bufl": "ERRNO)
      bytes+=length($0)+1
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
  say("Total jobs to run: counting...", 1)
  while ( bedcmd | getline > 0 ) {
    bufl++
    if ( bufl % 1000 == 0 ) {
      say("Total jobs to run: "bufl" and counting...", 1)
    }
    jobnum=bedline++ % arg::args["jobs"]
    buffer=tmpdir"/buffer_"jobnum".tmp" # number of pipes with regions == # jobs
    print $0 > buffer
    fflush()
  }
  for ( jobnum=0; jobnum < arg::args["jobs"]; jobnum++ ) { 
    buffer=tmpdir"/buffer_"jobnum".tmp" # number of pipes with regions == # jobs
    print SIGNAL_END_OF_BUFFER > buffer
    close(buffer)
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
    return piawka" win -v "arg::args["vcf"]
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
    ALT[nalt=0]=1

    # only consider single-char ALT that are not *
    while (x=index($5,",")) {
      nalt++
      if (x == 2 && substr($5,1,1)!="*") {
        ALT[nalt]=1
      }
      $5=substr($5,x+1)
    }
    if ( length($5) == 1 && $5 != "*" ) {
      ALT[++nalt]=1
    }

    if ( arg::args["persite"] == 1 ) { 
      start=$2-1
      end=$2
    }
  
    # Reset site-specific parameters
    for (i in groups) { miss[i]=0; misind[i]=0; n[i]=0; exclude[i]=0 }
    delete a
    delete thisnum
    delete thisden

    # Pool GT values for groups, count each state and missing data
    for (col in groupindex) {
      i=groupindex[col]
      if ( exclude[i] ) { continue }
      gtend=index( $col, ":" )
      if ( gtend==0 ) { gtend=length($col)+1 }
      for (c=1; c<gtend; c+=2) {
        al=substr($col,c,1)
        if ( al == "." ) {
          miss[i]+=gtend/2
          break
        } else {
          if ( !(al in ALT) ) { 
            exclude[i]=1
            break
          }
          a[i][al]++
          n[i]++
        }
      }
      if ( divide && i in a ) {
        if ( exclude[i] || ( arg::args["mult"] != 1 && length(a[i]) > 2 ) ) { 
          delete a[i]
          continue
        }
        calculate_within(i)
        calculate_between(i)
        if ( arg::args["persite"] ) { yield_output() }
      }
    }

    if ( divide ) { continue }
    for ( i in a ) {
      if ( exclude[i] || ( arg::args["miss"] < 1 && miss[i]/(miss[i]+n[i]) > arg::args["miss"] ) || ( arg::args["mult"] != 1 && length(a[i]) > 2 ) ) { 
        delete a[i]
        continue
      }
      calculate_within( i )
    }
    # Between is in second pass to make sure dependencies are there
    for ( i in a ) {
      calculate_between( i )
    }
    if ( arg::args["persite"] ) { yield_output() }
  }
  if ( !arg::args["persite"] ) { yield_output() }
  close(cmd)
}

function calculate_within(i) {
  if (!any_within) { return 0 }
  for ( s in stats::stats["within"] ) {
    incr="calc::increment_"s
    if ( incr in FUNCTAB ) {
      @incr(i)
    }
    num[i][s]+=thisnum[i][s]
    den[i][s]+=thisden[i][s] 
  }
}

function calculate_between(i) {
  if (!any_between) { return 0 }
  for ( j in a ) {
    if ( divide || j in pairs[i] ) { 
      # Populate a_ij array -- the union of a[i] and a[j] arrays
      # Also set a missing from a[i] or a[j] to zero
      delete a_ij
      for (x in a[i]) {
        if ( !(x in a[j]) ) { 
          a[j][x]=0 
        }
        a_ij[x]=a[i][x]
      }
      for (x in a[j]) {
        if ( !(x in a[i]) ) { 
          a[i][x]=0
        }
        a_ij[x]+=a[j][x] 
      }
      # If not arg::args["mult"] == 1, proceed only if common allele pool has <=2 alleles
      if ( arg::args["mult"] != 1 && length(a_ij) > 2 ) { return 0 }
      for ( s in stats::stats["between"] ) {
        incr="calc::increment_"s
        if ( incr in FUNCTAB ) {
          @incr(i,j)
        }
        num[i,j][s]+=thisnum[i,j][s]
        den[i,j][s]+=thisden[i,j][s] 
      }
    }
  }
}

function yield_output() {
  if ( pid>0 ) { return 0 }
  for (_ij in den) {
    if ( split(_ij, ij, SUBSEP) == 1 ) {
      i=ij[1]
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
          @fin(ij[1],ij[2])
        }
        if (s in stats::stats_print["between"]) {
          printOutput( ij[1], ij[2], s )
        }
      }
    }
  }
  # Reset locus-specific parameters
  delete num
  delete den
}

function printOutput( i, j, metric,    idx, ij ) {
  if (j=="") {
    idx=i
    j="."
  } else {
    idx=i SUBSEP j 
  }
  if ( den[idx][metric]==0 ) { return 0 }
  # reverse pop1 and pop2 to keep pop1 alphabetically smaller
  if ( divide && i > j ) {
    ij=j"\t"i
  } else {
    ij=i"\t"j
  }
}
  out=chr"\t"start"\t"end"\t"locus"\t"ij"\t"metric"\t"num[idx][metric]/den[idx][metric]"\t"num[idx][metric]"\t"den[idx][metric]
  print out > tmpf
}

END {
  if (pid>0) {
    total_jobs=bufl + (bufl % arg::args["jobs"] > 0 ? arg::args["jobs"]: 0 ) - 1
    jobs_width=int(log(total_jobs)/log(10))+1
    say(sprintf("Finishing job %*d of %d, seconds elapsed: %d", 
                jobs_width, 0, total_jobs, awk::systime()-time_start),
        1)
    for (f=0; f<bufl; f++) {
      tmpf=tmpdir"/"f".tmp"
      # to avoid race condition given persite, wait until the child writes to next region
      tmpfn=tmpdir"/"f+arg::args["jobs"]".tmp"
      while ( awk::stat(tmpfn, _) < 0 || _["size"]==0 ) {
        continue
      }
      say(sprintf("Finishing job %*d of %d, seconds elapsed: %d", 
                  jobs_width, f+arg::args["jobs"], total_jobs, awk::systime()-time_start),
          1)
      while ( getline < tmpf > 0 ) {
        if ( $0 == SIGNAL_END_OF_BUFFER ) {
          break
        }
        if ( piping_to_sum ) {
          sumcmd = piawka" sum -s "arg::args["stats"]" -"
          print | sumcmd
        } else {
          print
        }
      }
      close(tmpf)
      system("rm -f "tmpf)
    }
    # no need to wait since all children confirmedly written output?..
    say("") # clear stdout
    close(sumcmd)
  }
}

function say(string, same_line) {
  if (arg::args["quiet"] && string !~ /^Error/) { return }
  if (same_line) {
    printf "\033[2K\r" string > "/dev/stderr"
  } else {
    print string > "/dev/stderr"
  }
}
