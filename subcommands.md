---
layout: default
title: Subcommands
nav_order: 2
---

# Subcommands reference
{: .no_toc }

> [!WARNING]
> Code examples here are often piped to `head` to avoid terminal cluttering. Because of that, VCF processing ends abruptly and some error messages may pop up in addition to the outputs shown here. This is not to be worried about.

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
      The only mandatory argument is the VCF file (--vcf). Regions (--bed) and sample grouping (--groups) are recommended.
      Default --stats are pi,dxy; other can be provided as a spaceless comma-seprarted list (check `piawka list` for an overview).
      EXAMPLE:
          piawka calc [OPTIONS] -v file.tsv -b windows.bed -g ( groups.tsv | unite | divide ) > file.bed
  Options:
  -1, --persite       output values for each site
  -b, --bed <arg>     BED file with regions to be analyzed
  -d, --dependencies  output dependencies stats as well (best for piping to
                      `piawka sum`)
  -g, --groups <arg>  either 2-columns sample / group table or
                      keywords "unite" (all samples in one group) or "divide"
                      (each sample is a separate group); defaults to "unite"
  -j, --jobs <arg>    number of parallel jobs to run
  -m, --mult          use populations with multiple alleles at a site
  -q, --quiet         do not output progress and warning messages
  -R, --rand <arg>    randomly use this share of sites, 0.0-1.0
  -s, --stats <arg>   stats to calculate, comma-separated with no spaces as in
                      "pi,dxy,fst"; see `piawka list`
  -T, --targets <arg> BED file with targets (faster for numerous small regions)
  -v, --vcf <arg>     gzipped and tabixed VCF file
  -h, --help          show this help message and exit
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
- Neither: the VCF is processed as a single window per chromosome.

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
              -b genes.bed -g groups.tsv -s pi,dxy -q \
  | head -6 | column -t
#chr        start     end       locus      pop1          pop2              stat  value       numerator  denominator
scaffold_1  10035093  10035276  AL5G20950  CESiberia_2n  LE_2n             dxy   0.0071137   460        64664
scaffold_1  10035093  10035276  AL5G20950  PUWS_4n       .                 pi    0.00588993  640        108660
scaffold_1  10035093  10035276  AL5G20950  LE_2n         PUWS_4n           dxy   0.00881262  1102       125048
scaffold_1  10035093  10035276  AL5G20950  LE_2n         .                 pi    0.00772461  1078       139554
scaffold_1  10035093  10035276  AL5G20950  CESiberia_4n  UKScandinavia_2n  dxy   0.00722154  236        32680
```

**Per-site output:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g groups.tsv -s pi --persite -q \
              -b genes.bed | head -4 | column -t
#chr        start     end       locus      pop1              pop2  stat  value  numerator  denominator
scaffold_1  10035093  10035094  AL5G20950  PUWS_4n           .     pi    0      0          756
scaffold_1  10035093  10035094  AL5G20950  LE_2n             .     pi    0      0          992
scaffold_1  10035093  10035094  AL5G20950  UKScandinavia_2n  .     pi    0      0          90
```

**All samples as one group (default or with `-g unite`):**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -s pi -q | head -3 | column -t
#chr        start     end       locus      pop1         pop2  stat  value       numerator  denominator
scaffold_1  10035093  10035276  AL5G20950  all_samples  .     pi    0.00607559  5592       920404
scaffold_1  10035511  10036000  AL2G11890  all_samples  .     pi    0.0167103   36136      2162494
```

**Parallel jobs with Fst and multiple stats:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -g groups.tsv -s pi,dxy,fst -j 4 -q \
  | head -4 | column -t
#chr        start     end       locus      pop1          pop2   stat  value       numerator  denominator
scaffold_1  10035093  10035276  AL5G20950  CESiberia_2n  LE_2n  dxy   0.0071137   460        64664
scaffold_1  10035093  10035276  AL5G20950  CESiberia_2n  LE_2n  fst   0.130414    0.137013   1.0506
scaffold_1  10035093  10035276  AL5G20950  PUWS_4n       .      pi    0.00588993  640        108660
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

`piawka dist` reads the output of `piawka calc` and builds a pairwise distance matrix. If multiple rows exist for the same population pair (e.g., multiple windows), a weighted average is computed using the denominators. The matrix can be output in PHYLIP or Nexus format.

Only the pairwise statistic specified with `-s` is used; within-group rows are ignored. Diagonal entries are set to 0.

### Example

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              groups.tsv -s dxy -q \
  | piawka dist - | column -t
4
PUWS_4n           0.000000  0.010746  0.014952  0.010989
LE_2n             0.010746  0.000000  0.015055  0.010225
UKScandinavia_2n  0.014952  0.015055  0.000000  0.015170
CESiberia_2n      0.010989  0.010225  0.015170  0.000000
```

NEXUS output:

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              groups.tsv -s dxy -q \
  | piawka dist --nexus -
#NEXUS

BEGIN TAXA;
        DIMENSIONS ntax=4;
        TAXLABELS
 PUWS_4n LE_2n UKScandinavia_2n CESiberia_2n ;
END;

BEGIN DISTANCES;
        DIMENSIONS ntax=4;
        FORMAT TRIANGLE = BOTH LABELS = LEFT;
        Matrix
        PUWS_4n 0.000000 0.010746 0.014952 0.010989
        LE_2n 0.010746 0.000000 0.015055 0.010225
        UKScandinavia_2n 0.014952 0.015055 0.000000 0.015170
        CESiberia_2n 0.010989 0.010225 0.015170 0.000000

        ;
END;
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

`piawka filt` pivots the long-format `piawka calc` output to a wide row per (region, group pair) where all available statistics are present as columns and evaluates the expression. Rows where the expression is true are emitted back in long format. Available variables in the expression:

- **Field names**: `chr`, `start`, `end`, `locus`, `pop1`, `pop2`
- **Statistic names**: any stat name present in the input (e.g., `pi`, `dxy`, `fst`, `miss`)
- **Standard AWK operators**: `>`, `<`, `==`, `!=`, `&&`, `||`, `!`

When a statistic is not present for a particular row, its value is empty (evaluates as 0 in numeric context).

### Examples

**Keep only regions where pi > 0.01 in both groups of a pair:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -b genes.bed -g groups.tsv -s pi,dxy -q \
  | piawka filt -e 'pi > 0.01' \
  | head -6 | column -t
#chr        start     end       locus      pop1          pop2     stat  value      numerator  denominator
scaffold_1  10035511  10036000  AL2G11890  CESiberia_2n  LE_2n    dxy   0.0218687  3508       160412
scaffold_1  10035511  10036000  AL2G11890  LE_2n         .        pi    0.0150035  4716       314326
scaffold_1  10035511  10036000  AL2G11890  CESiberia_2n  .        pi    0.0253604  1988       78390
scaffold_1  10037631  10039653  AL1G31870  CESiberia_2n  LE_2n    dxy   0.0173026  9260       535180
scaffold_1  10037631  10039653  AL1G31870  LE_2n         PUWS_4n  dxy   0.014988   13788      919936
```

**Keep sites with no missing data:**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g groups.tsv -s pi,miss --persite -q \
  | piawka filt -e 'miss == 0' \
  | head -6 | column -t
#chr        start     end       locus                pop1          pop2  stat  value  numerator  denominator
scaffold_1  10000000  10000001  scaffold_1:10000001  LE_2n         .     miss  0      0          34
scaffold_1  10000000  10000001  scaffold_1:10000001  LE_2n         .     pi    0      0          1122
scaffold_1  10000000  10000001  scaffold_1:10000001  CESiberia_2n  .     miss  0      0          16
scaffold_1  10000000  10000001  scaffold_1:10000001  CESiberia_2n  .     pi    0      0          240
scaffold_1  10000000  10000001  scaffold_1:10000001  PUWS_4n       .     miss  0      0          28
```

**Keep sites that are fixed differences between two populations (privately differentially fixed):**

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g groups.tsv -s pi,dxy --persite -q \
  | piawka filt -e 'pi == 0 && dxy > 0' | head -6 | column -t
#chr        start     end       locus                pop1          pop2              stat  value  numerator  denominator
scaffold_1  10000444  10000445  scaffold_1:10000445  LE_2n         UKScandinavia_2n  dxy   1      340        340
scaffold_1  10000444  10000445  scaffold_1:10000445  CESiberia_2n  UKScandinavia_2n  dxy   1      160        160
scaffold_1  10000513  10000514  scaffold_1:10000514  LE_2n         UKScandinavia_2n  dxy   1      340        340
scaffold_1  10000513  10000514  scaffold_1:10000514  CESiberia_2n  UKScandinavia_2n  dxy   1      160        160
scaffold_1  10000950  10000951  scaffold_1:10000951  LE_2n         PUWS_4n           dxy   1      8          8
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
  -i, --ignore-chr      summarize across all chromosomes using only the locus field to match
  -I, --ignore-locus    ignore locus name (merge rows from all loci)
  -s, --stats <list>    stats to summarize (default: all stats found in input)
  -h, --help            show this help message and exit
```

### Description

`piawka sum` re-averages previously computed statistics while preserving numerical correctness. It accumulates numerators and denominators across rows that share the same chromosome + locus, then recomputes the value as `sum(numerator) / sum(denominator)`. For statistics with `finalize_*` functions (e.g., Tajima's D, Watterson's theta), it re-runs the finalization step if dependency statistics are present in the input.

**Use cases**:

- **Piping from `piawka calc`**: when `-b` is not given to `calc`, it outputs results in small arbitrary windows which are automatically piped to `sum` for chromosome-wide results.
- **Merging per-window results** across sliding windows or chromosomes.
- **Re-grouping individuals** into higher-level groups with `-g`.
- **Summarizing by custom regions** with `-b`.

### Examples

**Average per-region stats across the entire scaffolds:**

Here we summarize cluster-wise pi values along the entire VCF file.

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g groups.tsv -s pi -b genes.bed -q \
  | piawka sum --ignore-locus | column -t
#chr        start     end       locus  pop1              pop2  stat  value       numerator  denominator
scaffold_1  10035093  10496929  .      PUWS_4n           .     pi    0.00988641  1217900    123189256
scaffold_1  10035093  10496929  .      LE_2n             .     pi    0.00796761  1464568    183815248
scaffold_1  10035093  10496929  .      UKScandinavia_2n  .     pi    0.0077786   115084     14794952
scaffold_1  10035093  10496929  .      CESiberia_2n      .     pi    0.00711766  282906     39747074
```

**Re-group individuals into higher-level groups:**

The following is a way to get _average genic heterozygosity_ across the same groups.

```console
$ piawka calc -v alyrata_scaff_1_10000k-10500k.vcf.gz \
              -g divide -s pi -b genes.bed -q \
  | piawka sum --ignore-locus -g groups.tsv | column -t
#chr        start     end       locus  pop1              pop2  stat  value       numerator  denominator
scaffold_1  10035093  10496929  .      PUWS_4n           .     pi    0.00796982  112704     14141352
scaffold_1  10035093  10496929  .      LE_2n             .     pi    0.00734821  42128      5733096
scaffold_1  10035093  10496929  .      UKScandinavia_2n  .     pi    0.00581118  9748       1677456
scaffold_1  10035093  10496929  .      CESiberia_2n      .     pi    0.00695219  18812      2705910
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

**Target aggregation** with `--targets` is used internally by `piawka calc -T` to reduce the number of `tabix` queries when many small regions (e.g., single sites) are provided.

### Examples

**100 kb fixed windows from VCF header contigs:**

```console
$ piawka win -v alyrata_scaff_1_10000k-10500k.vcf.gz \
             -w 100000 | head -4
scaffold_1	0      	100000
scaffold_1	100000	200000
scaffold_1	200000	300000
scaffold_1	300000	400000
```

**Sliding 50 kb windows with 25 kb step:**

```console
$ piawka win -v alyrata_scaff_1_10000k-10500k.vcf.gz \
             -w 50000 -s 25000 | head -4
scaffold_1      0       50000
scaffold_1      25000   75000
scaffold_1      50000   100000
scaffold_1      75000   125000
```

**Aggregate exon targets for faster queries:**

```console
$ piawka win -T degenotate/fourfolds.bed | head -4
scaffold_1      10035268        10035948
scaffold_1      10037748        10039496
scaffold_1      10040888        10041816
scaffold_1      10048561        10054581
```
