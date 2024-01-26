#!/bin/bash

# Subset a big VCF (including non-variant sites!) by chunks of a BED file (genes/windows)
# and calculate pi within and dxy between groups specified by groups file.

# The bed file should have at least four columns (chr, start, end, locus name). gff2bed or bedtools makewindows -i winnum output suffices.
# Group file consists of sample and group name columns (see scripts/table_processing/create_groups.sh).

# This version produces long-format output suitable for plotting without much preprocessing.
# WARNING: Not yet tested in the submission system!

bed="/netscratch/dep_mercier/grp_novikova/Potamogeton/assembly/coloratus/final/annotation/windows_100k.bed"
vcf="/netscratch/dep_mercier/grp_novikova/Potamogeton/vcf/coloratus_for-pixy.vcf.gz"
groups="/netscratch/dep_mercier/grp_novikova/nikita/tables/groups/groups_phylo-species.tsv"
outfile=${vcf/.vcf.gz/_pixy_win100k_phylo-species_bial_testnew.tsv}

piscript="/netscratch/dep_mercier/grp_novikova/nikita/scripts/vcf_analysis/pi_final.awk"

echo $'locus\tnSites\tpop1\tpop2\tnUsed\tmetric\tvalue' > $outfile

# If gene names contain slashes, replace {//} with {= s@/.*@@ =} and {/} with {= s@^[^/]/+@@ =}
awk '{print $1 ":" $2+1 "-" $3 "/" $4}' $bed |
  parallel -j20 --bar \
    bcftools view -r {//} $vcf \| \
    $piscript -v LOCUS={/} $groups - >> $outfile

# Possible arguments and flags (passed before input files):
# -v LOCUS="locus_name"
# -v PERSITE=1 for per-site pi output
# -v MULT=1 for inclusion of multiallelic loci

