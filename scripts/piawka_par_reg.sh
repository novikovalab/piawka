#!/bin/bash

# Subset a VCF by chunks of a BED file (genes/windows)
# and calculate pi within and dxy between groups specified by groups file.

# Usage: piawka_par_reg.sh -a parallel_options -b bed_file -g groups_file -v vcf_gz -p piawka_options

# Options:
#
# -a parallel_options: a string of space-separated options for GNU parallel (e.g. -a "-j20 --block 10M")
#
# -b bed_file: the BED file with regions to analyze in parallel jobs.
#              If it contains 4+ columns, the 4th is passed as the locus name (LOCUS) to piawka_options
#              ( but can be overridden by LOCUS passed with -p . )
#              Thus, to suppress locus name with PERSITE=1, add -p "LOCUS=''" .
#
# -g grp_file: the groups file for piawka (see piawka docs).
#
# -p piawka_options: a string of space-separated options for piawka (e.g. -p "PIXY=1 MULT=1").
#
# -v vcf_gz: the compressed VCF file.

# Parse arguments
while getopts ":a:b:g:p:v:" opt; do
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
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# BED start field is 0-based and bcftools index -r option is not, so increment second field first
awk -v OFS="\t" '{$2++}1' $bed |
  parallel --colsep '\t' $paropts \
  bcftools view -r {1}:{2}-{3} $vcf \| \
  piawka LOCUS={4} NSITES='$(( {3} - {2} + 1 ))' $piopts $grp -


