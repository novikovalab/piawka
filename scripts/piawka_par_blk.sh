#!/bin/bash

# Run piawka on parallel for chunks of the VCF file and aggregate the result.

# Usage: piawka_par_blk.sh -g groups_file -v vcf_file -p piawka_options

# Options:
#
# -a parallel_options: a string of space-separated options for GNU parallel (e.g. -a "-j20 --block 10M")
#
# -g grp_file: the groups file for piawka (see piawka docs).
#
# -p piawka_options: a string of space-separated options for piawka (e.g. -p "PIXY=1 LOCUS=mylocus").
#
# -v vcf_file: the VCF file for piawka (see piawka docs).

# Parse arguments

while getopts ":a:b:g:p:v:" opt; do
  case $opt in
    a) paropts="${OPTARG}"
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

# default `-a` value
if [ -z ${paropts+x} ]; then paropts="--block 10M"; fi

# If gene names contain slashes, replace {//} with {= s@/.*@@ =} and {/} with {= s@^[^/]/+@@ =}
zcat $vcf | grep -v '^##' |
  parallel $paropts --pipe --header : \
  piawka $piopts $grp - |
  summarize_blks.awk

