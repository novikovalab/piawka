#!/bin/bash

print_help() {
echo \
'Calculate popgen statistics in a VCF file using `piawka` in parallel processes.

Usage: piawka_par.sh -a parallel_options -b bed_file -g groups_file -v vcf_gz -p piawka_options

Options:
-g grp_file: the groups file for piawka (see piawka docs).
-v vcf_gz: the compressed VCF file.

-a parallel_options: a string of space-separated options for GNU parallel (e.g. -a "--block 10M")
-b bed_file: the BED file with regions to analyze in parallel jobs.
             If it contains 4+ columns, the 4th is passed as the locus name (LOCUS) to piawka_options
             ( but can be overridden by LOCUS passed with -p . )
             Thus, to suppress locus name with PERSITE=1, add -p "LOCUS=''" .
-p piawka_options: a string of space-separated options for piawka (e.g. -p "PIXY=1 MULT=1").

-h : print this help message'
exit 0
}

while getopts ":a:b:g:hp:v:" opt; do
  case $opt in
    a) paropts="$OPTARG"
    ;;
    b) bed="$OPTARG"
    ;;
    g) grp="$OPTARG"
    ;;
    p) piopts="$OPTARG"
    ;;
    v) vcf="$OPTARG"
    ;;
    h | *) print_help
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# Help if no options provided
[ "$OPTIND" -eq 1 ] && print_help

if [ -z "$grp" ] || [ -z "$vcf" ]; then
  echo 'Mandatory arguments: -g groups_file -v vcf_gz' >&2
  exit 1
fi

# default parameter values
if [[ $piopts != *LOCUS=* ]] && [ -z "$bed" ]; then 
  piopts="LOCUS=$( basename $vcf .vcf.gz ) "$piopts
fi

if [ -z "$bed" ]; then
  zcat $vcf | grep -v '^##' |
  parallel $paropts --pipe --header : --block 10M --halt now,fail=1 \
  piawka VERBOSE=1 $piopts $grp - |
  { if [[ $piopts == *PERSITE=1* ]]; then cat -; else summarize_blks.awk -; fi; }
else
# BED start field is 0-based and bcftools index -r option is not, so increment second field first
  awk -v OFS="\t" '{$2++}1' $bed |
  parallel --colsep '\t' $paropts --halt now,fail=1 \
  tabix -h $vcf {1}:{2}-{3} \| \
  piawka LOCUS={4} NSITES='$(( {3} - {2} + 1 ))' $piopts $grp -
fi

