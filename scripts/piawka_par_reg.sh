#!/bin/bash

# Subset a VCF by chunks of a BED file (genes/windows)
# and calculate pi within and dxy between groups specified by groups file.

# Usage: piawka_par_reg.sh -b bed_file -g groups_file -v vcf_file -p piawka_options

# Options:
#
# -a parallel_options: a string of space-separated options for GNU parallel (e.g. -a "-j20 --block 10M")
#
# -b bed_file: the BED file with regions to analyze in parallel jobs.
#              If it contains 4+ columns, the 4th is passed as the locus name (LOCUS) to piawka.
#
# -g grp_file: the groups file for piawka (see piawka docs).
#
# -p piawka_options: a string of space-separated options for piawka (e.g. -p "PIXY=1 MULT=1").
#                    Note that LOCUS value, if provided, will be overridden.
#
# -v vcf_file: the VCF file for piawka (see piawka docs).

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

# If gene names contain slashes, replace {//} with {= s@/.*@@ =} and {/} with {= s@^[^/]/+@@ =}
awk '{if ($4) {$4="/"$4}; print $1 ":" $2+1 "-" $3 $4}' $bed |
  parallel $paropts \
  bcftools view -r {//} $vcf \| \
  piawka $piopts LOCUS={/} $grp -


