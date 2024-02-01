``piawka`` <img src="logo/logo.svg" align="right" width="25%">
==========

The powerful `awk` script to calculate π and Dxy (or πxy, or Nei's D) in VCF files.

Largely inspired by [`pixy`](https://github.com/ksamuk/pixy), it builds upon it in a few aspects:

 - **supports arbitrary ploidy level**, which can be different between samples and/or regions in the VCF
 - supports both average weighted π and `pixy`-like π calculation
 - can use multiallelic SNPs, in biallelic mode also uses multiallelic SNPs that have two alleles in the analyzed groups
 - **lightweight and portable**, runs wherever vanilla AWK can run (Windows, macOS, Linux...) and requires no installation
 - **faster on a single core** (and can be parallelized with shell tools, e.g. GNU `parallel` -- see [Usage](#usage))

By default, it reports `pixy`-like π and Dxy:

$$ π_{pixy}, Dxy_{pixy} = { \sum^n N_{diff} \over \sum^n N_{comp} } $$

Where $N_{diff}$ and $N_{comp}$ denote numbers of differences versus comparisons (within-group for π, between groups for Dxy) and $n$ stands for the number of sites used for calculation. This means that only one division per VCF file is performed after numerators and denominators from all sites are summarized. This metric gives lower weight to sites with fewer genotyped alleles (i.e. fewer possible comparisons) and should be more robust against missing data.

With option `PIXY=0` `piawka` will calculate **average weighted** π and Dxy like this:

$$ π_{w}, Dxy_{w} = { { \sum^n { N_{diff} \over N_{comp} } } \over n } $$

This metric might give unpredictable values at sites with lots of missing data, so we deliberately chose to only use sites with >50% alleles genotyped in the current group for weighted π and Dxy calculation.

## Running `piawka`

### Installation

Just clone the repo (or even simply download the `scripts` folder if you don't need 80Mb example VCF file) and you are good to go!

```
git clone https://github.com/novikovalab/piawka.git
cd piawka
```

We recommend the `mawk` AWK interpreter for `piawka` (it's much faster!). If you don't have it, just change the shebangs to `awk` like

```
mawk || sed -i '1s/mawk/awk/' ./scripts/piawka ./scripts/summarize_blks.awk
```

It might be useful to add `piawka` location to `PATH` environmental variable to run it from anywhere by either executing the following code or adding it to your `.bashrc` file:

```
echo 'export PATH="/path/to/piawka/scripts:$PATH"' >> ~/.bashrc
```

### Usage

To get help, run `piawka` without any arguments:

```
piawka
```

`piawka` works with decompressed VCF files. The most convenient way to use it is streaming the VCF file via stdin:

```
zcat file.vcf.gz | piawka [OPTIONS] groups_file - > piawka_pi-dxy.tsv
```
If you want to parallelize the counting and have GNU parallel installed, try our wrapper scripts:

```
# Parallelize VCF reading and count summary statistics for entire file
piawka_par_blk.sh -a parallel_options -g groups_file -p piawka_options -v vcf_gz

# Split VCF by BED regions and count stats for each region in parallel
piawka_par_reg.sh -a parallel_options -b bed_file -g groups_file -p piawka_options -v vcf_gz
```

See [Options](#options) and [Examples](#example-data) for further details.

### Input files

 - a **VCF file**: it is most sensible to include invariant sites for an unbiased estimation of π and Dxy. Consult the well-written [guide](https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html) by `pixy` authors. `piawka` only looks at what looks like discrete genotype calls, so any cell can have any imaginable number of genotypes. Regardless of the calculation method, `piawka` does not make assumptions about missing genotype calls and does not include them in the calculation (see why this is good [here](https://pixy.readthedocs.io/en/latest/about.html)). `piawka` parses multiallelic sites but does not parse sites with indels.
 - a **groups file**: this is a 2-column TSV file with no header, first column being the sample IDs from the VCF and the second being the group ID. `piawka` can handle arbitrary number of groups, which means one can also use it to calculate ploidy- and missing data-aware heterozygosity if all groups contain a single sample. If a sample from the VCF file is missing in the groups file, it will not be used for calculation. Groups file can also contain non-existent sample IDs, they will not be considered.

### Options

Options are provided as KEY=value pairs (no spaces around the `=` sign!) before input files. Flags can be set to 1 (true) or 0 (false). Following options exist:

 - `MULT=1` : counts pi and Dxy including multiallelic sites. Default is biallelic sites only. Higher values, lower comparability with other tools, but maybe more honest?
 - `PIXY=1` (default) : use the missingness-based site weighting as in [`pixy`](https://github.com/ksamuk/pixy). Might be better for groups with lots (>10%) of missing data according to [this paper](https://doi.org/10.1111/1755-0998.13707). `piawka` results might be slightly different (and more precise) because here we also make use of sites marked as multiallelic if they have two alleles in a given group. Thus, **`piawka` might be more suitable for multi-species VCF files** with higher share of multiallelic SNPs. Full convergence with `pixy` can be enforced by filtering out multiallelic sites before running `piawka` (e.g. `bcftools view -M2 file.vcf.gz`), or I can make it an option if there is demand for that.
 - `PERSITE=1` : returns per-site estimates instead of default VCF-wide average. Note that adding `PIXY=1` will not make any difference in this case.
 - `LOCUS="locus_name"` : the name of the locus in the output. Meaningless with `PERSITE=1`. Default is "chr\_start\_end" (first chromosome encountered in the file is taken).
 - `DXY=1` (default): output Dxy values along with pi values

Helper `parallel` scripts (`piawka_par_reg.sh` and `piawka_par_blk.sh`) accept following options:

- `-a parallel_options` : a string of space-separated options for GNU parallel (e.g. `-a "-j20"`)
- `-b bed_file` : the BED file with regions to analyze in parallel jobs.
              If it contains 4+ columns, the 4th is passed as the locus name (LOCUS) to piawka.
 - `-g grp_file` : the groups file for piawka (see piawka docs).
 - `-p piawka_options` : a string of space-separated options for piawka (e.g. -p "PIXY=1 MULT=1"). *Note that with `piawka_par_reg.sh` the LOCUS value, if provided, will be overridden.*
 - `-v vcf_gz` : the compressed VCF file for piawka (see piawka docs).

Here are examples of useful `parallel` options to be passed as `-a parallel_options`:

 - `-j 20` : number of parallel jobs (defaults to available CPUs)
 - `--block 10M` : for `piawka_par_blk.sh`, the size of the VCF block for 1 `piawka` job.
 - `--keep-order` : for `piawka_par_blk.sh`, make the order of lines in the output same as in non-parallel `piawka` run.
 - `--bar` : for `piawka_par_reg.sh`, display progress bar (the share of processed regions) in `stderr`.

Check the [`parallel` tutorial](https://www.gnu.org/software/parallel/parallel_tutorial.html) for more details on GNU `parallel`.

### Output

`piawka` outputs a long-format table with no header and following columns:

 - **locus** : either genomic position of analyzed locus or custom `LOCUS` value.
 - **nSites** : number of "potentially useful" lines in the VCF file (SNPs or invariant sites before filtering for number of alleles (and genotyping rate for weighted pi or dxy))
 - **pop1** : analyzed group 1
 - **pop2** : "." for pi values or group 2 for Dxy values
 - **nUsed** : number of sites used for pi calculation (i.e. SNPs and invariant sites; for weighted pi or dxy these should also pass the 50% genotyping rate threshold)
 - **metric** : pi or dxy, average weighted or `pixy`-like
 - **value** of the metric

## Example data

You can try `piawka` with the (part of) genomic variant data we made for Siberian *Arabidopsis lyrata* populations. `alyrata_scaff_1_10000k-10500k.vcf.gz` contains data for diploids and polyploids with various amounts of missing data split into several admixture groups defined in `groups.tsv` file. There are ~0.5M sites and 4 groups.

The following lazy time tests were run on a Lenovo laptop with a 12th gen Core i7.

Software versions: `bcftools==1.19`, `parallel==20230822`, `mawk==1.3.4` (but `piawka` should not be too version-dependent).

### Single-threaded execution (pure AWK)

In the most simple case, one would calculate summary pi and Dxy for each group and combination of groups for the entire region covered by the VCF file using site-weighted pi like this:

```
cd ./examples
vcf=alyrata_scaff_1_10000k-10500k.vcf.gz
grp=groups.tsv
out=piawka.tsv

zcat $vcf | piawka PIXY=0 $grp - > $out
head $out | column -t
```

The output is:

```

real    0m18.834s
user    0m18.816s
sys     0m0.887s

scaffold_1_9999942_10500000  310091  UKScandinavia_2n  .             215090  pi_w   0.00988236
scaffold_1_9999942_10500000  310091  PUWS_4n           .             219820  pi_w   0.010075
scaffold_1_9999942_10500000  310091  CESiberia_2n      .             215849  pi_w   0.00867701
scaffold_1_9999942_10500000  310091  LE_2n             .             218437  pi_w   0.00870315
scaffold_1_9999942_10500000  310091  PUWS_4n           CESiberia_2n  210149  dxy_w  0.0108729
scaffold_1_9999942_10500000  310091  UKScandinavia_2n  CESiberia_2n  202752  dxy_w  0.0149389
scaffold_1_9999942_10500000  310091  UKScandinavia_2n  PUWS_4n       205623  dxy_w  0.014567
scaffold_1_9999942_10500000  310091  LE_2n             CESiberia_2n  211095  dxy_w  0.0101943
scaffold_1_9999942_10500000  310091  UKScandinavia_2n  LE_2n         205216  dxy_w  0.0147753
scaffold_1_9999942_10500000  310091  PUWS_4n           LE_2n         213396  dxy_w  0.0106753
```

### Parallel execution (with GNU `parallel`)

Now, let's test the parallel VCF reading:

```
vcf=alyrata_scaff_1_10000k-10500k.vcf.gz
grp=groups.tsv
out2=piawka_blks.tsv
time piawka_par_blk.sh -a "-j20 --block 10M" -p "PIXY=0" -g $grp -v $vcf > $out2

#real    0m4.261s
#user    0m30.856s
#sys     0m3.198s
```

The output is identical to the first one except that the lines might be shuffled (you can change that by adding `--keep-order` `parallel` option.)

A common usecase is running `piawka` for a set of genomic regions (genes or windows). Having a BED file with these regions (like example `genes.bed`) at hand and helper tools like `bcftools` GNU `parallel` installed, this can be parallelized like:

```
vcf=alyrata_scaff_1_10000k-10500k.vcf.gz
grp=groups.tsv
bed=genes.bed
out3=piawka_genes.tsv

time piawka_par_reg.sh -a "-j20" -p "PIXY=0" -b $bed -g $grp -v $vcf > $out3

head -5 $out3 | column -t
```

```
real    0m3.900s
user    0m25.943s
sys     0m17.486s

AL5G20950  116  UKScandinavia_2n  .             85   pi_w   0.00630252
AL5G20950  116  PUWS_4n           .             103  pi_w   0.00341604
AL5G20950  116  CESiberia_2n      .             102  pi_w   0.0103426
AL5G20950  116  LE_2n             .             103  pi_w   0.0123319
AL5G20950  116  PUWS_4n           CESiberia_2n  102  dxy_w  0.00665266
```

-- and we have the stats for many genes after a single pass through the file in 4 seconds!

### Advanced usage example: genewise 4-fold and 0-fold sites' pi and Dxy

One can limit calculations to synonymous/non-synonymous sites inferred using an external tool. Example below was made with [`degenotate`](https://github.com/harvardinformatics/degenotate). At preparation step, `degeneracy-all-sites.bed` is made using `degenotate` with annotation file and reference genome sequence[^1]. Then following steps are needed to extract 0-folds and 4-folds and calculate genewise pi and dxy from them:

```
cd examples/degenotate
vcf=../alyrata_scaff_1_10000k-10500k.vcf.gz
grp=../groups.tsv

# Extract 0fold and 4fold BED lines
awk -v OFS="\t" '$5==0 { print $0 }' degeneracy-all-sites.bed > zerofolds.bed
awk -v OFS="\t" '$5==4 { print $0 }' degeneracy-all-sites.bed > fourfolds.bed

# Need unique gene names, so strip position info from NAME field of BED
zerogenes=( $( cut -f4 zerofolds.bed | sed 's/:[0-9]+//' | sort | uniq ) )
fourgenes=( $( cut -f4 fourfolds.bed | sed 's/:[0-9]+//' | sort | uniq ) )

# For each gene, extract part of VCF by BED and feed to piawka, name loci as genes
parallel -j20 \
  bcftools view -R \<( grep -w {} zerofolds.bed ) $vcf \| \
  piawka LOCUS={} $grp - > piawka_zerofolds.tsv ::: ${zerogenes[@]}

parallel -j20 \
  bcftools view -R \<( grep -w {} fourfolds.bed ) $vcf \| \
  piawka LOCUS={} $grp - > piawka_fourfolds.tsv ::: ${fourgenes[@]}
```

Another possibility is treating each line of `zerofolds.bed` and `fourfolds.bed` as a region for `piawka_par_reg.sh` and then summarize the result by genes using `summarize_blks.awk` helper script, but this way the result will be unweighted and might be less precise.

[^1] I was using [`liftoff`](https://github.com/agshumate/Liftoff) annotation that lacks phase info so I had to sanitize it first with `agat_sp_fix_cds_phase.pl` from [AGAT toolkit](https://github.com/NBISweden/AGAT). Details [here](https://github.com/harvardinformatics/degenotate/issues/32).

## Alternatives

This script shows all strong and weak points of being written in a simple text-processing language of `awk`. While combining high speed with not too complicated code, it lacks good error handling. It means that with corrupted inputs the script will do its best to produce *some* result silently. If you need a more fool-proof solution, consider some better-developed alternatives.

For diploid VCFs, one can use [`pixy`](https://github.com/ksamuk/pixy). To make it work with polyploids, one would need to randomly sample two GT values from each cell with >2 genotypes (it should not affect diversity metrics much at the genomic scale).

## References

The `pixy` method of dealing with missing data for π calculation was introduced in [this paper](https://doi.org/10.1111/1755-0998.13326). You can find discussion on the applicability of this metric [here](https://doi.org/10.1111/1755-0998.13707) and [here](https://doi.org/10.1111/1755-0998.13738). We suggest to parallelize `piawka` using GNU `parallel` which is introduced [here](https://doi.org/10.5281/zenodo.1146014).

## Citing `piawka`

If you want to express your gratitude for `piawka`, please cite our [Siberian *Arabidopsis* paper]() (coming soon!) where we have introduced and first used it.

