---
layout: default
title: Subcommands
nav_order: 2
---

# Subcommands reference
{: .no_toc }

<details open markdown="block">
  <summary>Contents</summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---

## calc

Calculate within- and between-group population statistics from a VCF file.

```
Usage: piawka calc [OPTIONS] -v file.vcf.gz [-g groups.tsv] [-b regions.bed]

    Calculate statistics within & between groups of samples in a VCF file.
    Mandatory arguments are VCF file (--vcf) and grouping (--groups). --bed is recommended.
    Default --stats are pi,dxy; other can be provided as a spaceless comma-separated list.
    EXAMPLE:
        piawka calc [OPTIONS] -v file.vcf.gz -g ( groups.tsv | unite | divide ) > file.bed

Options:
  -1, --persite          output values for each site
  -b, --bed <file>       BED file with regions to be analysed
  -d, --dependencies     output dependency stats as well (best for piping to piawka sum)
  -g, --groups <file>    2-column sample/group table, or keywords "unite" / "divide"
  -j, --jobs <n>         number of parallel jobs to run
  -m, --mult             use populations with multiple alleles at a site
  -q, --quiet            do not output progress and warning messages
  -R, --rand <frac>      randomly use this share of sites, 0.0–1.0
  -s, --stats <list>     comma-separated list of stats to calculate (see piawka list)
  -T, --targets <file>   BED file with targets (faster for many small regions)
  -v, --vcf <file>       bgzipped and tabix-indexed VCF file
  -h, --help             show this help message and exit
```

### Description

`piawka calc` is the core subcommand. It streams the VCF through `tabix` region queries and, for every site within each requested region, counts alleles per group. Per-site numerators and denominators are accumulated for each requested statistic, then a single summary row is printed per region per group (pair).

**Groups** can be specified three ways:

| Value | Behaviour |
|-------|-----------|
| `groups.tsv` | 2-column TSV mapping sample IDs to group names |
| `unite` | All VCF samples form one group named `all_samples` |
| `divide` | Each sample is its own group; all pairwise statistics are computed |

When `--groups` is omitted, `unite` is assumed.

**Regions** can be provided via:
- `-b / --bed`: a BED file with one row per region; column 4 sets the locus name in the output.
- `-T / --targets`: a BED file with small targets (e.g., CDS coordinates). `piawka` merges nearby targets into larger tabix queries automatically.
- Neither: the whole VCF is processed as a single window per chromosome.

**Parallelism**: with `-j N`, `piawka calc` forks `N` child processes. Each child reads from a shared FIFO buffer of region lines; results are written to temporary files and assembled in order by the parent before printing. See [Technical details](technical.html#parallelism).

**Output**: BED format with 10 columns:

```
#chr  start  end  locus  pop1  pop2  stat  value  numerator  denominator
```

`pop2` is `.` for within-group statistics.

### Examples

**Basic run — pi and dxy per gene:**

```console
$ cd piawka/examples
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -g groups.tsv -s pi,dxy \
  | head -6 | column -t
#chr        start     end       locus      pop1          pop2      stat  value       numerator  denominator
scaffold_1  10035093  10035276  AL5G20950  CESiberia_2n  LE_2n     dxy   0.0071137   460        64664
scaffold_1  10035093  10035276  AL5G20950  PUWS_4n       .         pi    0.00588993  640        108660
scaffold_1  10035093  10035276  AL5G20950  LE_2n         PUWS_4n   dxy   0.00881262  1102       125048
scaffold_1  10035093  10035276  AL5G20950  LE_2n         .         pi    0.00772461  1078       139554
scaffold_1  10035093  10035276  AL5G20950  CESiberia_4n  LE_2n     dxy   0.00721949  548        75888
```

**Per-site output for a single region:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g groups.tsv -s pi --persite \
              -b minitest.bed | head -4 | column -t
#chr        start     end       locus  pop1          pop2  stat  value    numerator  denominator
scaffold_1  10000112  10000113  sas    CESiberia_2n  .     pi    0.28571  8          28
scaffold_1  10000130  10000131  sas    LE_2n         .     pi    0.13333  2          15
scaffold_1  10000130  10000131  sas    PUWS_4n       .     pi    0.05882  2          34
```

**All samples as one group (no groups file):**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -s pi | head -3 | column -t
#chr        start     end       locus      pop1         pop2  stat  value       numerator  denominator
scaffold_1  10035093  10035276  AL5G20950  all_samples  .     pi    0.00735569  5088        691588
scaffold_1  10035511  10036000  AL2G11890  all_samples  .     pi    0.00831956  21240       2552940
```

**Parallel jobs with Fst and multiple stats:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -g groups.tsv -s pi,dxy,fst -j 4 -q \
  | wc -l
```

---

## dist

Convert pairwise `piawka calc` output to a PHYLIP or NEXUS distance matrix.

```
Usage: piawka dist [OPTIONS] file.bed

    Convert `piawka calc` pairwise distance values as a PHYLIP/NEXUS matrix.
    Input is piawka output with some pairwise statistics in.
    By default, dxy is used. Can be changed to other statistic (fst, rho etc.)
    If many values per sample pair are present, weighted average is taken.
    EXAMPLE:
        piawka dist [OPTIONS] file.bed > distmat.phy

Options:
  -s, --stat <name>   piawka pairwise stat to use as distance (default: dxy)
  -n, --nexus         write matrix in NEXUS format instead of PHYLIP
  -h, --help          show this help message and exit
```

### Description

`piawka dist` reads the output of `piawka calc` and builds a pairwise distance matrix. If multiple rows exist for the same population pair (e.g., multiple windows), a weighted average is computed using the denominators. The matrix can be output in:

- **PHYLIP** format (default): compatible with PHYLIP, RAxML, FastTree, and R's `ape` package.
- **NEXUS** format (`--nexus`): compatible with SplitsTree, SpecTRE `netmake`, and R's `phangorn`.

Only the pairwise statistic specified with `-s` is used; within-group rows are ignored. Diagonal entries are set to 0.

### Example

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -g groups.tsv -s dxy \
  | piawka dist -
5
CESiberia_2n  0.000000  0.009476  0.010348  0.009248  0.015265
CESiberia_4n  0.009476  0.000000  0.009521  0.006877  0.014528
LE_2n         0.010348  0.009521  0.000000  0.009551  0.015448
PUWS_4n       0.009248  0.006877  0.009551  0.000000  0.014237
SibA_2n       0.015265  0.014528  0.015448  0.014237  0.000000
```

NEXUS output (compatible with SplitsTree):

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -g groups.tsv -s dxy \
  | piawka dist --nexus -
#NEXUS

BEGIN TAXA;
    DIMENSIONS ntax=5;
    TAXLABELS CESiberia_2n CESiberia_4n LE_2n PUWS_4n SibA_2n ;
END;
...
```

---

## filt

Filter `piawka calc` output using AWK expressions on statistic values and field names.

```
Usage: piawka filt [OPTIONS] [file.bed ...]

    Filter `piawka` output by values of calculated statistics.
    BED entries for groups that satisfy the expression fully are kept.
    Expressions can include arithmetic, comparison and logical operators,
    numbers, statistics names and some field names (chr,start,end,locus,pop1,pop2).
    Takes piawka output file(-s) on stdin and an AWK-compatible expression.
    EXAMPLE:
        piawka filt -e 'pi > 0.01 && dxy > 0.01' file.bed > file_filt.bed

Options:
  -e, --expr <expr>   AWK expression to evaluate (default: 1, passes all lines)
  -h, --help          show this help message and exit
```

### Description

`piawka filt` pivots the long-format `piawka calc` output to a wide row per (region, group pair) and evaluates the expression. Rows where the expression is true are emitted back in long format. Available variables in the expression:

- **Field names**: `chr`, `start`, `end`, `locus`, `pop1`, `pop2`
- **Statistic names**: any stat name present in the input (e.g., `pi`, `dxy`, `fst`, `miss`)
- **Standard AWK operators**: `>`, `<`, `==`, `!=`, `&&`, `||`, `!`

When a statistic is not present for a particular row, its value is empty (evaluates as 0 in numeric context).

### Examples

**Keep only regions where pi > 0.005 in both groups of a pair:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -g groups.tsv -s pi,dxy \
  | piawka filt -e 'pi > 0.005' \
  | head -4 | column -t
#chr        start     end       locus      pop1          pop2  stat  value       numerator  denominator
scaffold_1  10035093  10035276  AL5G20950  CESiberia_2n  .     pi    0.00717748  652        90837
scaffold_1  10035093  10035276  AL5G20950  LE_2n         .     pi    0.00772461  1078       139554
scaffold_1  10035093  10035276  AL5G20950  PUWS_4n       .     pi    0.00588993  640        108660
```

**Keep sites with no missing data:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g groups.tsv -s pi,miss --persite \
  | piawka filt -e 'miss == 0' \
  | head -4 | column -t
```

**Keep sites that are fixed differences between two populations (privately differentially fixed):**

```console
$ piawka calc ... -s pi,dxy --persite \
  | piawka filt -e 'pi == 0 && dxy > 0' | head -4 | column -t
```

---

## list

List all statistics available for `piawka calc`.

```
Usage: piawka list [OPTIONS]

    List all statistics available for calculation with `piawka calc`.
    Takes no input.
    EXAMPLE:
        piawka list [OPTIONS]

Options:
  -d, --dependencies   output helper/dependency stats as well
  -h, --help           show this help message and exit
```

### Description

`piawka list` prints all available statistics with their descriptions, scope (within/between groups), and any dependencies on other statistics. By default, helper statistics (prefixed with `helper:`) are hidden; use `-d` to show them.

### Example

```console
$ piawka list
Available statistics:
lines      number of lines used in calculation
           	(within pops, dependencies: none)
miss       share of missing genotype calls
           	(within pops, dependencies: none)
pi         expected heterozygosity = nucleotide diversity
           	(within pops, dependencies: none)
maf        minor allele frequency
           	(within pops, dependencies: none)
daf        derived allele frequency
           	(within pops, dependencies: none)
theta_w    Watterson's theta
           	(within pops, dependencies: a1,segr,lines)
theta_low  Theta estimator based on sites with 0<allele_freq<0.33
           	(within pops, dependencies: lines)
theta_mid  Theta estimator based on sites with 0.33<=allele_freq<0.66
           	(within pops, dependencies: lines)
theta_high Theta estimator based on sites with allele_freq>=0.66
           	(within pops, dependencies: lines)
tajima     Tajima's D
           	(within pops, dependencies: a1,a2,segr,pi,lines,miss)
tajimalike Tajima's D interpolated for missing SNPs (experimental)
           	(within pops, dependencies: a1,a2,segrcorr,pi,lines)
afd        average allele frequency difference
           	(between pops, dependencies: none)
dxy        absolute nucleotide divergence
           	(between pops, dependencies: none)
fst        fixation index, Hudson's estimator
           	(between pops, dependencies: none)
fstwc      fixation index, Weir & Cockerham's estimator
           	(between pops, dependencies: none)
rho        Ronfort's rho
           	(between pops, dependencies: pi)
nei        Nei's D standard genetic distance
           	(between pops, dependencies: pi,dxy)
```

With `-d`, helper statistics (`a1`, `a2`, `segr`, `segrcorr`) are also shown.

---

## sum

Re-summarize `piawka calc` output across windows, genes, or any other grouping.

```
Usage: piawka sum [OPTIONS] [file.bed ...]

    Summarize `piawka calc` results counted over several loci.
    If dependencies are not given (see `piawka list`), defaults to sum(numerator)/sum(denominator).
    It only takes the output file(-s) of `piawka calc` passed over stdin.
    EXAMPLE:
        piawka sum [OPTIONS] file.bed > file_sum.bed

Options:
  -b, --bed <file>      summarize stats by regions in a BED file
  -g, --groups <file>   group file to average stats across individuals/subgroups
  -i, --ignore-chrs     summarize across all chromosomes using only the locus field to match
  -I, --ignore-locus    ignore locus name (merge rows from all loci)
  -s, --stats <list>    stats to summarize (default: all stats found in input)
  -h, --help            show this help message and exit
```

### Description

`piawka sum` re-averages previously computed statistics while preserving numerical correctness. It accumulates numerators and denominators across rows that share the same chromosome + locus, then recomputes the value as `sum(numerator) / sum(denominator)`. For statistics with `finalize_*` functions (e.g., Tajima's D, Watterson's theta), it re-runs the finalization step if dependency statistics are present in the input.

**Use cases**:

- **Piping from `piawka calc`**: when `-b` is not given to `calc`, it outputs per-chromosome windows which are automatically piped to `sum` for genome-wide gene-level results.
- **Merging per-window results** across sliding windows or chromosomes.
- **Re-grouping individuals** into higher-level groups with `-g`.
- **Summarizing by custom regions** with `-b`.

### Examples

**Summarize gene-level stats from a previously computed file:**

```console
$ piawka sum results.bed | head -4 | column -t
#chr        start     end       locus      pop1          pop2      stat  value       numerator  denominator
scaffold_1  10035093  10035276  AL5G20950  CESiberia_2n  LE_2n     dxy   0.0071137   460        64664
scaffold_1  10035093  10035276  AL5G20950  PUWS_4n       .         pi    0.00588993  640        108660
scaffold_1  10035511  10036000  AL2G11890  CESiberia_2n  LE_2n     dxy   0.00630281  3520       558536
```

**Average per-gene stats across all genes (whole scaffold):**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g groups.tsv -s pi,dxy \
  | piawka sum --ignore-locus | head -4 | column -t
```

**Re-group individuals into higher-level groups:**

```console
$ piawka sum -g subpop_to_region.tsv results.bed | head -4 | column -t
```

**Summarize calc output by custom windows:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g groups.tsv -s pi,dxy -b genes.bed \
  | piawka sum -b mRNA.bed | head -4 | column -t
```

---

## win

Partition a VCF or genome into BED windows for use with `piawka calc`.

```
Usage: piawka win [OPTIONS] -v file.vcf.gz | -T targets.bed

    Chunk a VCF file into windows for parallel processing or sub-chromosome resolution.
    Default (slow): output windows containing --lines [100000] VCF lines each.
    Alternatively, output windows of --window-size at --step
    (requires --vcf with ##contig header lines or a --fai file).
    If --targets are supplied instead of --vcf, produces bigger windows with multiple targets in.
    EXAMPLE:
        piawka win [OPTIONS] -v file.vcf.gz | -T targets.bed > win.bed

Options:
  -f, --fai <file>         FAI/TSV file with chr names and sizes (cols 1–2)
  -l, --lines <n>          VCF lines per window (default: 100000)
  -n, --number             enumerate windows in column 4
  -s, --step <bp>          step size for sliding windows (requires --window-size)
  -T, --targets <file>     BED file with small regions to aggregate into larger queries
  -v, --vcf <file>         bgzipped and tabix-indexed VCF file
  -w, --window-size <bp>   window size in bp (requires ##contig lines in VCF or --fai)
  -h, --help               show this help message and exit
```

### Description

`piawka win` produces BED windows that are then passed as the `-b` argument to `piawka calc`. There are three modes:

| Mode | Required input | How windows are defined |
|------|---------------|------------------------|
| VCF lines | `--vcf` | Fixed number of VCF lines per window |
| Fixed bp | `--vcf` with `##contig` headers, or `--fai` | Fixed base-pair width; optionally sliding with `--step` |
| Target aggregation | `--targets` | Nearby small targets merged into ~1 kb blocks |

**VCF-lines mode** (default) is the safest: it guarantees every window has approximately the same number of variant sites regardless of the recombination or variant density landscape.

**Fixed-bp mode** with `--window-size` produces equal-length windows, useful for population-genetic scans where physical scale matters. Requires either `##contig` lines in the VCF header or a FAI file.

**Target aggregation** with `--targets` is used internally by `piawka calc -T` to reduce the number of `tabix` queries when many small regions (e.g., exons) are provided.

### Examples

**100 kb fixed windows from VCF header contigs:**

```console
$ piawka win -v alyrata_scaff_1_10000k-10500k.vcf.gz \
             -w 100000 | head -4
scaffold_1	10000000	10100000
scaffold_1	10100000	10200000
scaffold_1	10200000	10300000
scaffold_1	10300000	10400000
```

**Sliding 50 kb windows with 25 kb step:**

```console
$ piawka win -v alyrata_scaff_1_10000k-10500k.vcf.gz \
             -w 50000 -s 25000 | head -4
scaffold_1	10000000	10050000
scaffold_1	10025000	10075000
scaffold_1	10050000	10100000
scaffold_1	10075000	10125000
```

**Aggregate exon targets for faster queries:**

```console
$ piawka win -T examples/cds.bed | head -4
scaffold_1	10035093	10035377
scaffold_1	10035511	10038165
scaffold_1	10038154	10039653
scaffold_1	10040878	10041907
```
