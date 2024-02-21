``piawka`` experimental branch
==========

Here I test some unstable features:

### Options

 - `FST=WC_ORIG` -- exact formula for Weir and Cockerham's Fst from the original paper. It gives the same result as `pixy` in default Fst mode, but I am not sure why the result is so much different from Hudson's Fst and from Bhatia's formula for WC estimator (which gives values way closer to Hudson-based).
 - `TAJLIKE=1` (experimental) : calculate Tajima's D-like statistic for a locus. Does not work with PERSITE=1 or MULT=1 or with mixed-ploidy groups. The only difference from classic formula is adjusting $\theta_W$ for missing data as suggested [here](https://doi.org/10.1534/genetics.112.139949) and [here](https://doi.org/10.4310/SII.2015.v8.n4.a4) (with `PIXY=1`, $\theta_\pi$ is adjusted too). Therefore, Tajima's denominator we use is not fully applicable here as the variance of the statistic is different, so the scaling of the function does not match Tajima's D. It means that absolute values cannot be compared with Tajima's D, and assessment of significance of deviation from zero is more complicated. Nevertheless, the sign of the metric conveys the same meaning as in Tajima's D, and the values are comparable across the genome. This estimator was already implemented in `rehh` R package [here](), but it is advised to use it with at least 80% genotyping rate (i.e. `MIS=0.8`).

