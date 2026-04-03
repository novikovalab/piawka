---
layout: default
title: Technical details
nav_order: 4
---

# Technical details
{: .no_toc }

<details open markdown="block">
  <summary>Contents</summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---

## SNP retrieval

### How piawka identifies usable sites

For each VCF line, `piawka calc` applies the following minimal filter:

1. **REF must be a single nucleotide** — multi-nucleotide reference alleles are skipped.
2. **At least one non-missing, non-star allele per group** — a site is used for a group only if that group has at least one called, non-`*` allele.
3. **Multiallelic handling** — by default, a site is skipped for a pair of groups if more than two distinct alleles are observed in their combined allele pool. This can be relaxed with `--mult`.

Crucially, piawka does **not** require that `$5` (the ALT field) contains a single-character non-star allele. Most tools that rely on REF/ALT checks mark a site as unusable if the whole-cohort ALT is multi-allelic (e.g., `A,T`). piawka instead checks the alleles *actually observed* in the target groups at runtime.

### The per-group ALT-agnostic advantage

Consider a VCF from 200 samples spanning 5 populations. At a given site the full cohort shows three ALT alleles (`A,T,G`). Standard tools skip this site entirely. But if populations A and B together carry only `A` and `T`, this is a clean biallelic SNP for those two populations. piawka recovers it.

```
Full cohort  →  A,T,G  (multiallelic, skipped by most tools)
Pop A + B    →  A,T    (biallelic, used by piawka)
Pop C + D    →  A,G    (biallelic, used by piawka)
```

This effect is strongest when analysing small subsets of a large population panel. In practice, the number of usable sites can increase by 10–40% over strict REF/ALT filtering when working with multi-population VCFs.

> **Illustration idea**: a diagram showing a multi-allelic site in a full VCF and how different per-group allele subsets yield biallelic comparisons. This would suit a hand-drawn schematic with three population ovals and allele counts inside each.

---

## The statistics framework

### Three lifecycle functions

Every statistic in piawka is defined by up to three AWK functions:

| Function | Called | Purpose |
|----------|--------|---------|
| `initiate_<stat>()` | Once, before reading the VCF | Validate arguments, print warnings, precompute constants |
| `increment_<stat>(i)` or `increment_<stat>(i,j)` | Once per usable VCF line, per group (pair) | Compute per-site numerator and denominator |
| `finalize_<stat>(i)` or `finalize_<stat>(i,j)` | Once per print call | Transform accumulated sums into the final value |

These functions are called indirectly by name via `@increment` / `@finalize` function-pointer dispatch in gawk. If a function is not defined, the default behaviour applies (see below).

### Numerators, denominators, and missing-data-aware averaging

Each `increment_*` function assigns:

- `thisnum[i]["stat"]` — the numerator contribution at this site
- `thisden[i]["stat"]` — the denominator contribution at this site

After the function returns, `calc.awk` accumulates:

```awk
num[i][s] += thisnum[i][s]
den[i][s] += thisden[i][s]
```

At print time, the value is `num[i][s] / den[i][s]`. If no `finalize_*` function is defined, this ratio is computed directly.

**Why this is correct for missing data**: consider `pi`. At a site with 8 called alleles, the denominator is 8×7 = 56. At a site with 4 called alleles (4 missing), the denominator is 4×3 = 12. The region-level pi is:

```
pi = (num_site1 + num_site2 + ...) / (56 + 12 + ...)
```

Sites with fewer called alleles naturally receive less weight — without any explicit imputation or site dropping. This is identical in spirit to the "ratio of averages" approach used in vcftools pixy.

### The finalize_* function

Some statistics cannot be computed simply as `sum(num)/sum(den)`. For example, Watterson's theta requires:

```
theta_w = (S / lines) / a1
```

where `a1` (the harmonic number) was accumulated across segregating sites. The `finalize_theta_w` function reads `num[i]["segr"]` and `num[i]["a1"]` (both accumulated by their own `increment_*` functions) and assembles the final value by overwriting `num[i]["theta_w"]` and `den[i]["theta_w"]` before the ratio is printed.

This late-binding design means that `theta_w` is always correct even when `piawka sum` re-runs it on pre-accumulated numerators/denominators — as long as the dependency statistics (`a1`, `segr`, `lines`) are also present.

> **Illustration idea**: a diagram showing VCF sites flowing into increment functions for pi and segr, with arrows pointing to finalize_theta_w that combines their accumulated values.

---

## Statistic dependencies

### How dependencies work

When a statistic is requested (e.g., `-s tajima`), `piawka` automatically resolves and includes all dependencies listed in the fourth argument of `stats::add_stat(...)`:

```awk
BEGIN {
  stats::add_stat("tajima", "Tajima's D", 0, "a1,a2,segr,pi,lines,miss")
}
```

The `stats::parse_stats()` function recursively parses the dependency list, marks all dependency statistics as "used", and sorts the full list so that each statistic is computed after its dependencies. `--dependencies` also marks them as "printed" so they appear in the output.

The `increment_*` and `finalize_*` functions of dependency statistics are guaranteed to run before those of the dependent statistic because the ordering of `@include` statements in the `piawka` file follows the dependency order.

### Why this matters for piawka sum

When `piawka calc` is called without `-b`, it pipes its output directly to `piawka sum`. `sum` re-reads the pre-accumulated numerators and denominators and re-runs the `finalize_*` functions. For this to work correctly, all dependency statistics must be present in the piped stream — which is why `calc` automatically enables `--dependencies` in this mode.

---

## Adding statistics modules

New statistics can be added without modifying any existing file. Follow the template in `include/example.awk`:

1. **Create** a new file, e.g., `include/mystat.awk`.
2. **Register** the statistic in the `BEGIN` block:
   ```awk
   @namespace "calc"
   BEGIN {
     stats::add_stat(
       "mystat",                   # name used with -s
       "My description",           # shown by piawka list
       0,                          # 0 = within-group, nonzero = between-group
       "pi,lines"                  # comma-separated dependencies (or "")
     )
   }
   ```
3. **Optionally define**:
   - `initiate_mystat()` — run once before VCF processing
   - `increment_mystat(i)` (or `increment_mystat(i,j)` for between-group) — run per site
   - `finalize_mystat(i)` (or `finalize_mystat(i,j)`) — run before printing
4. **Add** `@include "mystat.awk"` to the `piawka` script, **after** all its dependencies.

Inside `increment_*`, access the following per-site variables:

| Variable | Type | Description |
|----------|------|-------------|
| `a[i][x]` | array | allele counts: index *x* = allele ID (integer), value = count |
| `n[i]` | integer | total called alleles in group *i* |
| `miss[i]` | integer | missing allele calls in group *i* |
| `a_ij[x]` | array | (between-group only) union allele counts across groups *i* and *j* |

Inside `finalize_*`, the following are available (in addition to all `increment_*` variables):

| Variable | Description |
|----------|-------------|
| `num[i]["stat"]` | accumulated numerator for any used statistic (including dependencies) |
| `den[i]["stat"]` | accumulated denominator |

If a `finalize_*` function is defined, it should assign to `num[i]["mystat"]` and `den[i]["mystat"]` before returning. The printed value will be `num/den`.

> **Illustration idea**: a flowchart of the three function calls (initiate → increment per site → finalize → print), perhaps with the leech mascot holding a formula card.

---

## Parallel processing

### Architecture

`piawka calc -j N` uses `gawk`'s `@load "fork"` extension to create N child processes:

```
Parent process
├── Reads BED file, writes region lines to per-job FIFO buffers
├── Waits for all N children to finish
└── Assembles output in order from per-job temporary files

Child process 0          Child process 1          ...
├── Reads from buffer_0  ├── Reads from buffer_1
├── Runs tabix per region├── Runs tabix per region
└── Writes to 0.tmp      └── Writes to 1.tmp
```

### Buffer and checkpoint files

All inter-process communication goes through the temporary directory (`/tmp/piawkatmpXXXXXXXX`):

| File | Purpose |
|------|---------|
| `buffer_N.tmp` | FIFO-like file: parent writes region lines; child reads them |
| `M.tmp` | Results from job M, written by child, read by parent |
| `M_started.ckp` | Empty checkpoint: child signals it started job M |
| `M_fin.ckp` | Empty checkpoint: child signals job M is complete |
| `N_done.ckp` | Empty checkpoint: child N has finished all jobs |

The parent polls checkpoint files with `awk::sleep(0.1)` intervals. Once all `_fin.ckp` files are present, the parent reads the corresponding `.tmp` files in order and prints them to stdout. A background shell command automatically removes the temporary directory when the parent process exits.

### What is distributed

Each unit of work sent to a child process is **one BED region** from the `-b` file (or one auto-generated window). Within a child, `tabix` retrieves all VCF lines overlapping that region, and the full statistics computation runs serially for that region. With large BED files (thousands of regions), parallelism provides near-linear speedup up to the number of available CPU cores. With few large regions, the gain is limited.

The output order is preserved: results are always printed in the same order as the input BED regardless of which child finished first.

> **Illustration idea**: a simple pipeline diagram showing parent distributing regions to child boxes, each child writing a tmp file, and the parent reading them back in order.
