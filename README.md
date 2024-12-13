``piawka`` <img src="logo/logo.svg" align="right" width="25%">
==========

[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/piawka?style=for-the-badge&label=bioconda%20downloads&color=FA9E89)](http://bioconda.github.io/recipes/piawka/README.html)

The powerful `awk` script to calculate π, Dxy (or πxy, or Nei's D) and some more simple stats (Fst, Tajima's D, Ronfort's rho) in VCF files in the command line. Developed to analyze arbitrary-ploidy groups with substantial amounts of missing data.

> :warning: `piawka` is under development. If something does not seem to work well, check for updates and do not hesitate to file an issue!

[wiki](https://github.com/novikovalab/piawka/wiki)

# Install it!

## Quickly: `conda`

```bash
conda install -c bioconda piawka
```

## Quickly but slower

Make the following programs available in the command line (install and add to `PATH`):

 - `gawk` v5.0.0 and above 
 - `tabix`

Then, get `piawka` by cloning the repo and add the scripts to `PATH`:

```bash
git clone https://github.com/novikovalab/piawka.git
export PATH="$( realpath ./piawka/scripts ):${PATH}"
```

# Use it!

```console
$ piawka
piawka v0.8.10
Usage:
piawka -g groups_tsv -v vcf_gz [OPTIONS]
Options:
-1, --persite       output values for each site
-a, --all           output more cols: numerator, denominator, nGeno, nMiss
-b, --bed <arg>     BED file with regions to be analyzed
-B, --targets <arg> BED file with targets (faster for numerous small regions)
-D, --nodxy         do not output Dxy
-f, --fst           output Hudson Fst
-F, --fstwc         output Weir and Cockerham Fst instead
-g, --groups <arg>  either 2-columns sample / group table or 
                    keywords "unite" (1 group) or "divide" (n_samples groups)
-h, --help          show this help message
-H, --het           output only per-sample pi = heterozygosity
-j, --jobs <arg>    number of parallel jobs to run
-m, --mult          use multiallelic sites
-M, --miss <arg>    max share of missing GT per group at site, 0.0-1.0
-P, --nopi          do not output pi
-q, --quiet         do not output progress and warning messages
-r, --rho           output Ronfort's rho
-t, --tajimalike    output TajimaD-like stat (manages missing data but untested)
-T, --tajima        output classic TajimaD instead (affected by missing data)
-v, --vcf <arg>     gzipped and tabixed VCF file
```

See the [wiki](https://github.com/novikovalab/piawka/wiki) for further details.

# Cite it!

If you want to express your gratitude for having `piawka`, please cite our [Siberian *Arabidopsis* paper](https://www.biorxiv.org/content/10.1101/2024.08.27.609292) where we have introduced and first used it.

