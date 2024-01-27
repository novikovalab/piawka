``piawka`` <img src="logo/logo.png" align="right" width="20%">
==========

The powerful `awk` script to calculate π and Dxy (or πxy, or Nei's D) in VCF files.

Largely inspired by [`pixy`](https://github.com/ksamuk/pixy), it builds upon it in a few aspects:

 - **supports arbitrary ploidy level**, which can be different between samples and/or regions in the VCF
 - supports both average weighted π and `pixy`-like π calculation
 - can use multiallelic SNPs, in biallelic mode also uses multiallelic SNPs that have two alleles in the analyzed groups
 - **lightweight and portable**, runs wherever vanilla AWK can run and requires no installation
 - **faster on a single core** (and can be parallelized to some extent with shell tools, e.g. GNU `parallel`)

By default, it reports **average weighted** π and Dxy that are calculated like this:

$$ π_{w}, Dxy_{w} = { { \sum^n { N_{diff} \over N_{comp} } } \over n } $$

Where $N_{diff}$ and $N_{comp}$ denote numbers of differences versus comparisons (within-group for π, between groups for Dxy) and $n$ stands for the number of sites used for calculation. This metric might give unpredictable values at sites with lots of missing data, so we deliberately chose to only use sites with >50% alleles genotyped in the current group for π and Dxy calculation.

With option `PIXY=1` `piawka` will calculate `pixy`-like π and Dxy:

$$ π_{pixy}, Dxy_{pixy} = { \sum^n N_{diff} \over \sum^n N_{comp} } $$

— this means that only one division per VCF file is performed once numerators and denominators from all sites are summarized. This metric gives lower weight to sites with fewer genotyped alleles (i.e. fewer possible comparisons) and should be more robust against missing data.

## Running `piawka`

### Installation

Just clone the repo (or even simply download `piawka` if you don't need 80Mb example VCF file) and you are good to go!

```
git clone https://github.com/novikovalab/piawka.git
cd piawka
```

We recommend the `mawk` AWK interpreter for `piawka` (it's much faster!). If you don't have it, just change the shebang to `awk` like

```
mawk || sed -i '1s/mawk/awk/' ./piawka
```

It might be useful to add `piawka` location to `PATH` environmental variable to run it from anywhere by either executing the following code or adding it to your `.bashrc` file:

```
echo 'export PATH="/path/to/piawka/scripts:$PATH"' >> ~/.bashrc
```

### Usage

`piawka` works with decompressed VCF files. The most convenient way to use it is streaming the VCF file via stdin:

```
zcat file.vcf.gz | piawka [OPTIONS] groups_file - > piawka_pi-dxy.tsv
```
If you want to parallelize the counting and have GNU parallel installed, try our wrapper scripts:

```
# Parallelize VCF reading and count summary statistics for entire file
piawka_par_blk.sh -g groups_file -p piawka_options -v vcf_gz

# Split VCF by BED regions and count stats for each region in parallel
piawka_par_reg.sh -b bed_file -g groups_file -p piawka_options -v vcf_gz
```

See [Options](#options) and [Examples](#example-data) for further details.

### Input files

 - a **VCF file**: it is most sensible to include invariant sites for an unbiased estimation of π and Dxy. Consult the well-written [guide](https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html) by `pixy` authors. `piawka` only looks at what looks like discrete genotype calls, so any cell can have any imaginable number of genotypes. Regardless of the calculation method, `piawka` does not make assumptions about missing genotype calls and does not include them in the calculation (see why this is good [here](https://pixy.readthedocs.io/en/latest/about.html)).
 - a **groups file**: this is a 2-column TSV file with no header, first column being the sample IDs from the VCF and the second being the group ID. `piawka` can handle arbitrary number of groups, which means one can also use it to calculate missing data-aware heterozygosity if `groups_file` has two identical columns with unique sample names. If a VCF file is missing in the groups file, it will not be used for calculation. Groups file can contain non-existent sample IDs, they will not be considered.

### Options

Options are provided as KEY=value pairs before input files. Following options exist:

 - `MULT=1` : counts pi and Dxy including multiallelic sites. Default is biallelic sites only. Higher values, lower comparability with other tools, but more honest and insightful.
 - `PIXY=1` : use the missingness-based site weighting as in [`pixy`](https://github.com/ksamuk/pixy). Might be better for groups with lots (>10%) of missing data according to [this paper](https://doi.org/10.1111/1755-0998.13707). `piawka` results might be slightly different (and more precise) because here we also make use of sites marked as multiallelic if they have two alleles in a given group. Thus, **`piawka` might be more suitable for multi-species VCF files** with higher share of multiallelic SNPs. Full convergence with `pixy` can be enforced by changing `$5 !~ /\*|[ACGT]{2}/` to `$5 !~ /\*|,|[ACGT]{2}` in the script (I can make it an option if there is demand for that).
 - `PERSITE=1` : returns per-site estimates instead of default VCF-wide average. Note that adding `PIXY=1` will not make any difference in this case.
 - `LOCUS="locus_name"` : the name of the locus in the output. Meaningless with `PERSITE=1`. Default is "chr\_start\_end" (first chromosome encountered in the file is taken).

Helper scripts (`piawka_par_*`) accept following options:

- `-a parallel_options` : a string of space-separated options for GNU parallel (e.g. `-a "-j20"`)
- `-b bed_file` : the BED file with regions to analyze in parallel jobs.
              If it contains 4+ columns, the 4th is passed as the locus name (LOCUS) to piawka.
 - `-g grp_file` : the groups file for piawka (see piawka docs).
 - `-p piawka_options` : a string of space-separated options for piawka (e.g. -p "PIXY=1 MULT=1"). *Note that with `piawka_par_reg.sh` the LOCUS value, if provided, will be overridden.*
 - `-v vcf_file` : the VCF file for piawka (see piawka docs).

Here is a subset of `parallel` options to be passed as `-a parallel_options`:

 - `-j [num]` : number of parallel jobs (defaults to available CPUs)
 - `--block 10M` : for `piawka_par_blk.sh`, the size of the VCF block for 1 `piawka` job. Bigger blocks => higher RAM usage.
  - `--keep-order` : for `piawka_par_blk.sh`, make the order of lines in the output same as in non-parallel `piawka` run.

### Output

`piawka` outputs a long-format table with no header and following columns:

 - **locus** : either genomic position of analyzed locus or custom `LOCUS` value.
 - **nSites** : number of lines in the VCF file (not all of which might be useful, e.g. in case of indels, for default non-weighted pi measure sites with <50% genotyping rate for a given group are also excluded)
 - **pop1** : analyzed group 1
 - **pop2** : "." for pi values or group 2 for Dxy values
 - **nUsed** : number of sites used for pi calculation (i.e. SNPs and invariant sites)
 - **metric** : pi or dxy, average weighted or `pixy`-like
 - **value** of the metric

## Example data

You can try `piawka` with the (part of) genomic variant data we made for Siberian *Arabidopsis lyrata* populations. `alyrata_scaff_1_10000k-10500k.vcf.gz` contains data for diploids and polyploids with various amounts of missing data split into several admixture groups defined in `groups.tsv` file. There are ~0.5M sites and 4 groups.

The following lazy time tests were run on a Lenovo laptop with a 12th gen Core i7.

In the most simple case, one would calculate summary pi and Dxy for each group and combination of groups using site-weighted pi like this:

```
cd ./examples
vcf=alyrata_scaff_1_10000k-10500k.vcf.gz
grp=groups.tsv
out=piawka.tsv

zcat $vcf | piawka $grp - > $out
head $out
```

The output is:

```

real    0m18.834s
user    0m18.816s
sys     0m0.887s

scaffold_1_9999942_10500000     310091  UKScandinavia_2n        .       215090  pi_w    0.00988236
scaffold_1_9999942_10500000     310091  PUWS_4n .       219820  pi_w    0.010075
scaffold_1_9999942_10500000     310091  CESiberia_2n    .       215849  pi_w    0.00867701
scaffold_1_9999942_10500000     310091  LE_2n   .       218437  pi_w    0.00870315
scaffold_1_9999942_10500000     310091  PUWS_4n CESiberia_2n    210149  dxy_w   0.0108729
scaffold_1_9999942_10500000     310091  UKScandinavia_2n        CESiberia_2n    202752  dxy_w    0.0149389
scaffold_1_9999942_10500000     310091  UKScandinavia_2n        PUWS_4n 205623  dxy_w   0.014567
scaffold_1_9999942_10500000     310091  LE_2n   CESiberia_2n    211095  dxy_w   0.0101943
scaffold_1_9999942_10500000     310091  UKScandinavia_2n        LE_2n   205216  dxy_w   0.0147753
scaffold_1_9999942_10500000     310091  PUWS_4n LE_2n   213396  dxy_w   0.0106753
```

Now, let's test the parallel VCF reading:

```
out2=piawka_blks.tsv
time piawka_par_blk.sh -a "-j20 --block 10M" -g $grp -v $vcf > $out2

#real    0m4.261s
#user    0m30.856s
#sys     0m3.198s
```

The output is identical to the first one except that the lines might be shuffled (you can change that by adding `--keep-order` `parallel` option.)
A common usecase is running `piawka` for a set of genomic regions (genes or windows). Having a BED file with these regions (like example `genes.bed`) at hand and helper tools like `bcftools` GNU `parallel` installed, this can be parallelized like:

```
bed=genes.bed
out3=piawka_genes.tsv

time piawka_par_reg.sh -a "-j20" -b $bed -g $grp -v $vcf > $out3

head -5 $out3
```


```
real    0m5.570s
user    0m34.291s
sys     0m6.471s

y_num_ID=AL5G20950.v2.1_0       116     UKScandinavia_2n        .       85      pi_w    0.00630252
y_num_ID=AL5G20950.v2.1_0       116     PUWS_4n .       103     pi_w    0.00341604
y_num_ID=AL5G20950.v2.1_0       116     CESiberia_2n    .       102     pi_w    0.0103426
y_num_ID=AL5G20950.v2.1_0       116     LE_2n   .       103     pi_w    0.0123319
y_num_ID=AL5G20950.v2.1_0       116     PUWS_4n CESiberia_2n    102     dxy_w   0.00665266
```

-- and we have the stats for many genes after a single pass through the file!

## Alternatives

This script shows all strong and weak points of being written in a simple text-processing language of `awk`. While combining high speed with not too complicated code, it lacks good error handling. It means that with corrupted inputs the script will do its best to produce *some* result silently. If you need a more fool-proof solution, consider some better-developed alternatives.

For diploid VCFs, one can use [`pixy`](https://github.com/ksamuk/pixy). To make it work with polyploids, one would need to randomly sample two GT values from each cell with >2 genotypes (it should not affect diversity metrics much at the genomic scale).

## References

The `pixy` method of dealing with missing data for π calculation was introduced in [this paper](https://doi.org/10.1111/1755-0998.13326). You can find discussion on the applicability of this metric [here](https://doi.org/10.1111/1755-0998.13707) and [here](https://doi.org/10.1111/1755-0998.13738). We suggest to parallelize `piawka` using GNU `parallel` which is introduced [here](https://doi.org/10.5281/zenodo.1146014).

## Citing `piawka`

If you want to express your gratitude for `piawka`, please cite our [Siberian *Arabidopsis* paper]() where we have introduced and first used it.

