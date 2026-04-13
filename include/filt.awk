@namespace "filt"

function run() { 
  help="\
    Filter `piawka` output by values of calculated statistics. \n\
    BED entries for groups that satisfy the expression fully are kept. \n\
    Expressions can include arithmetic, comparison and logical operators, \n\
    numbers, statistics names and some field names (chr,start,end,locus,pop1,pop2). \n\
    Takes piawka output file(-s) on stdin and an AWK-compatibe expression. \n\
    EXAMPLE: \n\tpiawka filt -e 'pi > 0.01 && dxy > 0.01' file.bed > file_filt.bed"
  arg::add_argument("e", "expr", 0, "expression to be evaluated")
  arg::parse_args(2, help, "no help if empty")
  narg=arg::parse_nonargs()
  if ( !("expr" in arg::args) ) {
    arg::args["expr"]=1
  }
  if ( narg==0 ) {
    arg::nonargs[++narg]="/dev/stdin"
    calc::say("Warning: no input files given, reading from stdin!")
  }

  calc::print_header()
  parse_expr()
  prepare_awkscript()
  for (n=1; n<=narg; n++) {
    filter_regions(arg::nonargs[n])
  }
  close(awkscript)
  exit 0
}

function parse_expr() {
  split("chr start end locus pop1 pop2", fieldsi)
  for (i in fieldsi) {
    fields[fields[i]]=1
  }
  split(arg::args["expr"], ex, / +|\^|\*\*|\*|\/|%|\+|-|\(|\)|!=|!|&&|\|\||>=|<=|>|<|==|=/)
  for (i in ex) {
    if ( ex[i] ~ /[0-9]+(\.[0-9]+)?/ ) {
      if ( ex[i] == 0 && !made_zero_warning++ ) {
        calc::say("Warning: checks against zero should be quoted (e.g. pi==\"0\" and not pi==0)! Otherwise empty string matches the criterion as well") 
      }
    } else if (ex[i] in stats::list) {
      seenstat[ex[i]]=1
      stats=stats","ex[i]
    } else if ( !(ex[i] == "" || ex[i] in fields || ex[i] in FUNCTAB) ) {
      calc::say("Warning: not a recognized statistic/field name or value:" ex[i])
    }
  }
  stats::parse_stats(substr(stats,2))
}

function prepare_awkscript(    vars, s) {
  close(awkscript)
  vars=""
  header="chr\tstart\tend\tlocus\tpop1\tpop2"
  for ( si in stats::used ) {
    s = stats::used[si]
    vars = vars" -v "s"="
    header = header "\t" s
  }
  awkscript=" gawk -v chr= -v start= -v end= -v locus= -v pop1= -v pop2= " vars " 'NR==1{split($0,s)} \n\
              NR>1{for(i=1;i<=NF;i++){split($i,f,SUBSEP);SYMTAB[s[i]]=f[1]} \n\
                if(!(" arg::args["expr"] ")){next} \n\
                if($0 ~ SUBSEP SUBSEP SUBSEP){for(i=1;i<=NF;i++){if($i!~SUBSEP SUBSEP SUBSEP){continue};split($i,f,SUBSEP);SYMTAB[s[i]]=f[6]} \n\
                  if(!(" arg::args["expr"] ")){next}} \n\
                for(i=7;i<=NF;i++){if($i==\"\"){continue};if($i!~SUBSEP SUBSEP SUBSEP){gsub(SUBSEP,\"\\t\",$i);print $1,$2,$3,$4,$5,$6,s[i],$i}} \n\
              }' FS=\\\\t OFS=\\\\t -"
  print header | awkscript
}

function filter_regions(f,    firstline) {
  firstline=1
  while (getline < f > 0) {
    if (firstline) { 
      firstline=0
      piawka::assert( NF == 10, "piawka output (10 columns) is required!" )
    }
    if ( substr($0,1,1)!="#" ) { # exclude header 

      # print if new locus/chr
      if ( locus != $1 "\t" $2 "\t" $3 "\t" $4 && locus != "" ) {
        print_wide()
      }

      if ( !seenstat[$7] ) {
        seenstat[$7]=1
        stats=stats","$7
        stats::parse_stats(substr(stats,2))
        prepare_awkscript()
      }

      locus=$1 "\t" $2 "\t" $3 "\t" $4  # locus-chr
      stat[$5,$6][$7]=$8 SUBSEP $9 SUBSEP $10
    }
  }
  print_wide()
  close(f)
}

function print_wide() {
  for (pop in stat) {
    split(pop, p1p2, SUBSEP)
    out = locus "\t" p1p2[1] "\t" p1p2[2]
    for ( si in stats::used ) {
      s = stats::used[si]
      if (p1p2[2]!="." && si in stats::within ) {
        out = out"\t"stat[p1p2[1] SUBSEP "."][s] SUBSEP SUBSEP SUBSEP stat[p1p2[2] SUBSEP "."][s]
      } else {
        out = out"\t"stat[pop][s]
      }
    }
    print out | awkscript
  }
  delete stat
  locus=""
}
