---
layout: default
title: Overview
nav_order: 1
description: "piawka – calculate pi and other population statistics in VCF files"
permalink: /
---

# piawka <img src="logo/logo.svg" align="right" width="180">

**piawka** (pronounced pee-af-kah = *leech* in Russian) calculates SNP-based population statistics for arbitrary groups of samples directly from bgzipped, tabixed VCF files. It outputs a [BED](https://en.wikipedia.org/wiki/BED_(file_format)) file that can be intersected with any genomic annotation and piped into the companion subcommands for filtering, summarizing, and distance matrix creation.

[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/piawka?style=flat-square&label=bioconda%20downloads&color=FA9E89)](http://bioconda.github.io/recipes/piawka/README.html)
[![GitHub](https://img.shields.io/badge/source-GitHub-181717?style=flat-square&logo=github)](https://github.com/novikovalab/piawka)

---

## What piawka does

`piawka calc` reads a VCF file and, for each genomic window or locus supplied in a BED file, counts alleles per group and per site, then accumulates per-site numerators and denominators for each requested statistic. At the end of each window the accumulated sums are emitted as a single BED row per group (or group pair) per statistic. The companion subcommands transform or aggregate this output for downstream analysis.

### Example output

```console
$ cd piawka/examples
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz -b genes.bed -g groups.tsv -s pi,dxy \
  | head -5 | column -t
#chr        start     end       locus      pop1          pop2      stat  value       numerator  denominator
scaffold_1  10035093  10035276  AL5G20950  CESiberia_2n  LE_2n     dxy   0.0071137   460        64664
scaffold_1  10035093  10035276  AL5G20950  PUWS_4n       .         pi    0.00588993  640        108660
scaffold_1  10035093  10035276  AL5G20950  LE_2n         PUWS_4n   dxy   0.00881262  1102       125048
scaffold_1  10035093  10035276  AL5G20950  LE_2n         .         pi    0.00772461  1078       139554
```

Each output row carries a **numerator** and **denominator** in addition to the final value. This is essential for correctly averaging statistics across windows (see [Technical details](technical.html)).

---

## Advantages

### Indexable BED output

The output is sorted, tab-delimited BED with a standard header. It can be `bgzip`-compressed & indexed with `tabix` and fed to other genomics programs with no or little restructuring. The `locus` field allows labeling results by gene, exon, or any other genomic feature.

### Correct handling of missing data

All statistics are computed as sums of per-site numerators and denominators, with additional summarization on top if needed. A site with many missing genotypes contributes a smaller numerator-denominator pair than a site with full coverage. When the output is later averaged over windows or genes with `piawka sum`, the averaging is automatically missing-data aware — no bias is introduced.

### Support for polyploid variant calls

`piawka` counts individual alleles rather than diploid genotypes. Mixed-ploidy cohorts are handled natively; sample ploidy is detected automatically from the first VCF data line.

### Higher data yield: ALT-agnostic SNP retrieval

Many tools check if a site's REF and ALT alleles are both single nucleotides to identify it as a biallelic SNP. When analysing a subset of samples from a large VCF, many sites biallelic in the subset may be multi-allelic in the full dataset. `piawka` checks only that the REF allele is a single nucleotide; it then counts only the alleles actually present in each group. This recovers substantially more SNPs for small subsets of large VCF files. See [Technical details](technical.html#snp-retrieval) for a worked example.

### Optional support for multiallelic sites 

The flexible retrieval of genotypes allows `piawka` to process sites with more than one alternative allele. Many statistics in `piawka` are defined in a way that naturally extends to multiallelic sites; this allows to explore genetic patters less accessible to classic tools.

### Extensible statistics

`piawka` ships with {{ site.data.stats_count | default: "18" }} statistics covering nucleotide diversity, divergence, differentiation, allele frequency spectra and standard neutrality tests. Any user-defined statistic can be added as a plain AWK file following a simple interface. See [Technical details](technical.html#adding-modules) and `include/example.awk` for further details.

### GNU AWK: no installation, fast, low memory

The entire tool is a single AWK/shell script with no compiled dependencies beyond `gawk ≥ 5.2` with default extensions, `tabix`, and `bgzip`. It runs at competitive speed and uses modest amounts of RAM. Multiprocessing is built in via `gawk`'s `fork()` extension.

---

## Installation

```bash
# Recommended: conda
conda install -c bioconda piawka

# Alternative: clone and add to PATH
git clone https://github.com/novikovalab/piawka.git
export PATH="$( realpath ./piawka ):${PATH}"
# Dependencies: gawk>=5.2.0, tabix, bgzip
```

---

## Input and output

**Required** (for `piawka calc`):

| Input | Description |
|-------|-------------|
| VCF file | bgzipped (`.vcf.gz`) and tabix-indexed |

**Optional**:

| Input | Description |
|-------|-------------|
| Groups file | 2-column TSV: sample ID → group ID OR keywords "unite" (all samples in one group) or "divide" (each sample is a separate group) |
| BED file (`-b`) | Regions to analyse; optionally with locus name in column 4 |
| Targets BED (`-T`) | Small regions for targeted extraction (faster than `-b` for many small regions); optionally with locus name in column 4 |

**Output** columns:

| Column | Description |
|--------|-------------|
| `chr`, `start`, `end` | BED coordinates |
| `locus` | Region name (from BED col. 4, or `.` if absent) |
| `pop1` | Group name (within-group stats) or first group (pairwise) |
| `pop2` | `.` for within-group stats; second group for pairwise |
| `stat` | Statistic name |
| `value` | Statistic value |
| `numerator` | Raw numerator (for re-averaging with `piawka sum`) |
| `denominator` | Raw denominator (for re-averaging with `piawka sum`) |

---

## Subcommands

| Subcommand | Description |
|------------|-------------|
| [`calc`](subcommands.html#calc) | Calculate within- and between-group statistics from a VCF file |
| [`dist`](subcommands.html#dist) | Convert pairwise `calc` output to a PHYLIP or NEXUS distance matrix |
| [`filt`](subcommands.html#filt) | Filter `calc` output rows using AWK expressions on statistic values |
| [`list`](subcommands.html#list) | List all statistics available for `piawka calc` |
| [`sum`](subcommands.html#sum) | Re-summarize `calc` output across custom regions or super-groups |
| [`win`](subcommands.html#win) | Partition a VCF or genome into windows for `calc` |

---

## Citation

If you use `piawka`, please cite:

> Scott et al. (2025) *Mol. Biol. Evol.* [https://doi.org/10.1093/molbev/msaf153](https://doi.org/10.1093/molbev/msaf153)
