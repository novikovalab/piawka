``piawka`` <img src="logo/logo.svg" align="right" width="25%">
==========

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/piawka/README.html)

The powerful `awk` script to calculate π, Dxy (or πxy, or Nei's D) and some more simple stats (Fst, Tajima's D, Ronfort's rho) in VCF files in the command line. Developed to analyze arbitrary-ploidy groups with substantial amounts of missing data.

[wiki](https://github.com/novikovalab/piawka/wiki)

# Install it!

## Quickly: `conda`

```bash
conda install -c bioconda piawka
```

## Quickly but slower

Make the following programs available in the command line (install and add to `PATH`):

 - `awk` (we recommend `gawk` >= 1.3.4)
 - `tabix` (optional, for region-based analyses)
 - `parallel` (optional, for the parallel wrapper `piawka_par.sh`)

Then, get `piawka` by cloning the repo and add the scripts to `PATH`:

```bash
git clone https://github.com/novikovalab/piawka.git
export PATH="$( realpath ./piawka/scripts ):${PATH}"
```

**Note** that the shebang is set to `gawk` in all AWK scripts. If you have another AWK or if your `/usr/bin/env` does not support the `-S` option, change the shebangs.

# Use it!

`piawka` takes following input data:

 - **VCF file**
 - **groups file** (two columns: sample ID and group ID, tab or space as separator)
 - (optional) **BED file** to restrict analysis to certain regions (windows/genes etc.)

Most usecases are covered by the wrapper script that runs `piawka` in parallel:

```bash
piawka_par.sh -a parallel_options -g groups_file -p piawka_options -v vcf_gz [ -b bed_file ]
```

`parallel_options` include all arguments passed to the `parallel` command, e.g. `-a "-j 16"` to run `piawka` in 16 threads.

`piawka_options` include analysis parameters such as calculated statistics. Run `piawka` with no arguments to get the list of options:

```bash
piawka
```

See the [wiki](https://github.com/novikovalab/piawka/wiki) for further details.

# Cite it!

If you want to express your gratitude for having `piawka`, please cite our [Siberian *Arabidopsis* paper](https://www.biorxiv.org/content/10.1101/2024.08.27.609292) where we have introduced and first used it.

