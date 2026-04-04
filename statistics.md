---
layout: default
title: Statistics
nav_order: 3
---

# Statistics reference
{: .no_toc }

<details open markdown="block">
  <summary>Contents</summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

This page describes every statistic that `piawka calc` can compute. Run `piawka list` to see the current list; run `piawka list -d` to include dependencies = helper statistics.

Statistics marked **within** are computed per group; those marked **between** are computed per pair of groups.

---

## Within-group statistics

### daf — derived allele frequency

**Scope**: within &emsp; **Dependencies**: none

The frequency of non-REF (derived) alleles at each site, averaged across all SNP sites in the region:

$$\text{daf} = \frac{\sum_\text{sites} (\text{total alleles} - \text{REF allele count})}{\sum_\text{sites} \text{total alleles}}$$

The interpretation as "derived" frequency assumes the VCF is **polarized**: that is, the REF allele corresponds to the ancestral state. Check polarization carefully before interpreting `daf`; see also `theta_low`, `theta_mid`, `theta_high`.

---

### lines — site count

**Scope**: within &emsp; **Dependencies**: none

The number of VCF lines (sites) that were used in the calculation for a given region and group. A site is used if at least one allele was observed in the group. `lines` is primarily a helper statistic used by `tajima` and `theta_w`, but it can be useful on its own as a measure of site availability in the dataset.

---

### maf — minor allele frequency

**Scope**: within &emsp; **Dependencies**: none

The frequency of the least common allele at each site, averaged across all SNP sites:

$$\text{maf} = \frac{\sum_\text{sites} \min_\text{alleles}(\text{allele count})}{\sum_\text{sites} \text{total alleles}}$$

Unlike `daf`, `maf` does not require a polarized VCF.

---

### miss — missingness

**Scope**: within &emsp; **Dependencies**: none

The fraction of missing genotype calls across all sites in the region:

$$\text{miss} = \frac{\sum_\text{sites} \text{missing allele calls}}{\sum_\text{sites} (\text{called alleles} + \text{missing allele calls})}$$

`miss` is useful for quality filtering: `piawka filt -e 'miss < 0.1'` retains only windows with fewer than 10% missing calls. It is also an internal dependency of `tajima`.

---

### pi — nucleotide diversity

**Scope**: within &emsp; **Dependencies**: none

Expected heterozygosity per site (or absolute nucleotide diversity, or $\hat{\theta}_{\pi}$, or Tajima's π): the probability that two randomly drawn alleles differ at a given site, averaged across all sites.

$$\pi = \frac{\sum_\text{sites} \left(n^2 - \sum_\text{alleles} c_a^2\right)}{\sum_\text{sites} n(n-1)}$$

where *n* is the number of called alleles and *c_a* is the count of allele *a* at a site. The formula follows [Korunes & Samuk (2021) *Mol. Ecol. Resources*](https://doi.org/10.1111/1755-0998.13326), which clarifies the correct handling of missing data.

This is the most commonly used summary of within-population genetic variation. It is robust to missing data because sites with fewer called alleles contribute a smaller denominator.

**Literature**: [Tajima (1983) *Genetics* 105:437](https://www.genetics.org/content/105/2/437).

---

### tajima — Tajima's D

**Scope**: within &emsp; **Dependencies**: `a1`, `a2`, `segr`, `pi`, `lines`, `miss`

Tajima's D is the standardized difference between two estimators of the population mutation rate θ:

$$D = \frac{\hat\theta_\pi - \hat\theta_W}{\sqrt{\mathrm{Var}(D)}}$$

- **θ̂_π** = `pi` × `lines` (sum of pairwise differences)
- **θ̂_W** = S / a₁ where S = number of segregating sites and a₁ is the (n-1)-th harmonic number
- **Var(D)** uses the classic Tajima (1989) formula but a₁ and a₂ are computed per-site weighted by sample size (see [Technical details](technical.html#missing-data-aware-averaging) for the missing-data treatment)

| D value | Interpretation |
|---------|----------------|
| D < 0 | Excess of rare variants → possible directional selection or expansion |
| D ≈ 0 | Neutral evolution under constant population size |
| D > 0 | Excess of intermediate-frequency variants → possible balancing selection or bottleneck |

Cannot be computed with `--persite`.

**Literature**: [Tajima (1989) *Genetics* 123:585](https://www.genetics.org/content/123/3/585).

---

### tajimalike — Tajima's D with missing-data correction (experimental)

**Scope**: within &emsp; **Dependencies**: `a1`, `a2`, `segrcorr`, `pi`, `lines`

An experimental variant of Tajima's D that uses `segrcorr` (the expected number of segregating sites given the observed missingness) instead of the raw count. This corrects the numerator for the downward bias introduced by missing data, at the cost of additional approximation. Treat results with caution and compare with `tajima`.

---

### theta_w — Watterson's theta

**Scope**: within &emsp; **Dependencies**: `a1`, `segr`, `lines`

Watterson's estimator of the population mutation rate θ per site:

$$\hat\theta_W = \frac{S / a_1}{\text{lines}}$$

where S is the number of segregating sites and a₁ = Σ (1/k) for k = 1…(n−1) is the (n−1)-th harmonic number. When sample sizes vary across sites (due to missing data), a₁ is averaged across segregating sites weighted by the per-site sample size.

Cannot be computed with `--persite`.

**Literature**: [Watterson (1975) *Theor. Popul. Biol.* 7:256](https://doi.org/10.1016/0040-5809(75)90020-9).

---

### theta_low — low-frequency theta estimator

**Scope**: within &emsp; **Dependencies**: `lines`

A theta estimator based only on sites where the derived allele frequency is in the low range (0 < DAF < 1/3). Computed as:

$$\hat\theta_\text{low} = \frac{\sum_\text{low-freq sites} \text{alt} \cdot w}{\text{lines}}$$

where `alt` = non-REF allele count and *w* = 1/⌊n/3⌋ is a weight that normalizes by the number of frequency classes. Requires a **polarized** VCF. Not defined for multiallelic comparisons.

The relative magnitude of `theta_low`, `theta_mid`, and `theta_high` summarizes the shape of the site frequency spectrum and is sensitive to selection and demographic history.

**Literature**: [Achaz (2009) *Genetics* 183:249](https://doi.org/10.1534/genetics.109.104075).

---

### theta_mid — intermediate-frequency theta estimator

**Scope**: within &emsp; **Dependencies**: `lines`

A theta estimator based only on sites with intermediate derived allele frequency (1/3 ≤ DAF < 2/3). Same weight formula as `theta_low`. Requires a polarized VCF.

**Literature**: [Achaz (2009) *Genetics* 183:249](https://doi.org/10.1534/genetics.109.104075).

---

### theta_high — high-frequency theta estimator

**Scope**: within &emsp; **Dependencies**: `lines`

A theta estimator based only on sites where the derived allele is at high frequency (DAF ≥ 2/3, i.e. REF count ≤ n/3). Same weight formula as `theta_low`. Requires a polarized VCF.

**Literature**: [Achaz (2009) *Genetics* 183:249](https://doi.org/10.1534/genetics.109.104075).

---

## Between-group statistics

### afd — average allele frequency difference

**Scope**: between &emsp; **Dependencies**: none

The average absolute difference in allele frequency between two groups, per site:

$$\text{afd} = \frac{\sum_\text{sites} \sum_\text{alleles} |f_{i,a} - f_{j,a}|}{2 \cdot \sum_\text{sites} n_i n_j / (n_i + n_j)}$$

where *f* is allele frequency and *n* is the number of called alleles. `afd` ranges from 0 (identical allele frequencies) to 1 (complete fixation of different alleles). Unlike `fst`, it is not normalized by within-group diversity and is therefore more directly interpretable as a measure of raw allelic differentiation.

Not reliable in multiallelic comparisons.

**Literature**: [Berner (2019) *Genes* 10(4):308](https://doi.org/10.3390/genes10040308).

---

### dxy — absolute nucleotide divergence

**Scope**: between &emsp; **Dependencies**: none

The probability that two alleles drawn one from each group differ, averaged over all sites:

$$d_{xy} = \frac{\sum_\text{sites} \left( n_i n_j - \sum_\text{alleles} c_{i,a} \cdot c_{j,a} \right)}{\sum_\text{sites} n_i n_j}$$

where *c_{i,a}* is the count of allele *a* in group *i*. `dxy` is an absolute divergence measure: it does not depend on within-group diversity, so it can increase even when `fst` is constant (e.g., when both groups simultaneously become more diverse).

**Literature**: [Nei & Li (1979) *PNAS* 76:5269](https://doi.org/10.1073/pnas.76.10.5269); the estimator formula follows [Korunes & Samuk (2021)](https://doi.org/10.1111/1755-0998.13326).

---

### fst — fixation index (Hudson's estimator)

**Scope**: between &emsp; **Dependencies**: none

The proportion of genetic variation that is explained by group membership, using Hudson's moment estimator:

$$F_{ST} = \frac{\sum_\text{sites} (H_B - H_W)}{\sum_\text{sites} H_B}$$

where H_W = average within-group heterozygosity (average of the two groups' per-site pi) and H_B = between-group heterozygosity (equivalent to dxy at a single site). Averaging numerators and denominators separately across sites gives the correct "ratio of averages" estimator that is robust to differences in sample size.

Ranges from 0 (panmixia) to 1 (complete reproductive isolation).

**Literature**: [Hudson et al. (1992) *Genetics* 132:153](https://www.genetics.org/content/132/1/153); see also [Bhatia et al. (2013)](https://doi.org/10.1101/gr.154831.113) for an evaluation of estimators.

---

### fstwc — fixation index (Weir & Cockerham's estimator)

**Scope**: between &emsp; **Dependencies**: none

Weir & Cockerham's estimator of Fst, following the simplified form of Bhatia et al. (2013) eq. (6):

$$F_{ST}^{WC} = \frac{\sum_\text{sites} \left[ \alpha - 2\beta \cdot \text{mism} \right]}{\sum_\text{sites} \alpha}$$

where α and β are functions of the two group sample sizes and mism is the sum of within-group heterozygosity terms. This estimator accounts for unequal and potentially small sample sizes more accurately than the Hudson estimator in such cases.

**Literature**: [Weir & Cockerham (1984) *Evolution* 38:1358](https://doi.org/10.2307/2408641); simplified form from [Bhatia et al. (2013)](https://doi.org/10.1101/gr.154831.113).

---

### nei — Nei's D standard genetic distance

**Scope**: between &emsp; **Dependencies**: `pi`, `dxy`

Nei's standard genetic distance, computed from the identity measures J_X (= 1 − π_i), J_Y (= 1 − π_j), and J_XY (= 1 − d_xy):

$$D = -\ln \frac{J_{XY}}{\sqrt{J_X \cdot J_Y}}$$

This distance ranges from 0 (identical populations) and increases without bound. It is proportional to divergence time under an infinite-allele neutral model.

**Literature**: [Nei (1972) *Am. Nat.* 106:283](https://doi.org/10.1086/282771).

---

### rho — Ronfort's ρ

**Scope**: between &emsp; **Dependencies**: `pi`

An analog of Fst for polyploid populations that corrects for the ploidy-induced inflation of within-group heterozygosity (H_sp):

$$\rho = \frac{H_{pt} - H_s}{H_{pt} - H_{sp}}$$

where H_s is the average within-group heterozygosity, H_t is the between-group heterozygosity, and H_sp is H_s corrected for ploidy. For diploids, ρ ≈ Fst. For autotetraploids, `rho` is lower than naive `fst` because it removes the upward bias caused by within-individual fixed heterozygosity.

Requires that all samples in a group have the same ploidy. Not reliable in multiallelic comparisons.

**Literature**: [Ronfort et al. (1998) *Genetics* 150:921](https://www.genetics.org/content/150/2/921); [Meirmans & Van Tienderen (2004) *Mol. Ecol.* 4:394](https://doi.org/10.1111/j.1471-8286.2004.00590.x).

---

## Helper statistics

Helper statistics are internal building blocks used by other statistics. They are hidden by default in `piawka list`; use `piawka list -d` to show them. They can be requested with `-s` and included in output with `--dependencies`.

### a1 — 1st harmonic number

**Scope**: within &emsp; **Dependencies**: none

At each segregating site, increments the accumulated a₁ by the (n−1)-th harmonic number H_{n-1} = Σ (1/k) for k = 1…n-1, where *n* is the per-site allele count. Used by `theta_w` and `tajima` to normalize the number of segregating sites.

---

### a2 — 2nd harmonic number

**Scope**: within &emsp; **Dependencies**: none

At each biallelic segregating site, increments the accumulated a₂ by the (n−1)-th second-order harmonic number H_{n-1}^{(2)} = Σ (1/k²) for k = 1…n-1. Used by `tajima` to compute the variance of D.

---

### segr — segregating site count

**Scope**: within &emsp; **Dependencies**: none

Counts the number of sites in a region that are polymorphic (have more than one allele) in the group. Used by `theta_w` and `tajima`.

---

### segrcorr — corrected segregating site count

**Scope**: within &emsp; **Dependencies**: none

An expected count of segregating sites, corrected upward for the missing data at each site. If a site is segregating in *n₁* out of *n₂* possible allele calls (the rest missing), the corrected contribution S_expected is computed from the ratio of harmonic coefficients for n₁ and n₂. Used by `tajimalike`.
