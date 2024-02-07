#!/bin/bash

# Run piawka on parallel for chunks of the VCF file and aggregate the result.

# Usage: piawka_par_blk.sh -a parallel_options -g groups_file -v vcf_gz -p piawka_options

# Options:
#
# -a parallel_options: a string of space-separated options for GNU parallel (e.g. -a "-j20 --block 100M")
#                      default: "--block 10M" (parallel default 1M is too small for a genomic VCF)
#
# -g grp_file: the groups file for piawka (see piawka docs).
#
# -p piawka_options: a string of space-separated options for piawka (e.g. -p "PIXY=1 LOCUS=mylocus").
#                    default: "LOCUS=[vcf basename]"
#
# -v vcf_gz: the compressed VCF file.

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

# default parameter values
if [ -z ${paropts+x} ]; then paropts="--block 10M"; fi
if [[ $piopts != *LOCUS=* ]]; then piopts="LOCUS=$( basename $vcf .vcf.gz ) "$piopts; fi

zcat $vcf | grep -v '^##' |
  parallel $paropts --pipe --header : \
  piawka $piopts $grp - |
  { if [[ $piopts == *PERSITE=1* ]]; then cat -; else summarize_blks.awk -; fi; }

