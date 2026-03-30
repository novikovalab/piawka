# This file gives an example of how piawka can be expanded with further statistics.
# It is directly usable after filling the BEGIN block 
# as well as initiate_* and/or increment_* and/or finalize_* functions.

# Key points:
#  - BEGIN has to have a stats::add_stat("mystat") call describing 
#    the scope and dependencies of the statistic if any.
#  - function initiate_mystat() is run in `piawka calc` ONLY ONCE before processing the VCF.
#  - function increment_mystat() is run ONCE PER VCF LINE that has a usable SNP.
#  - function finalize_mystat() is run ONCE PER PRINT CALL for that statistic.
#  - further functions may be supplied for use in the above functions.
#  - values from "dependencies" statistics can be reused in functions
#      if the dependencies are included before mystat in the piawka file.
#
# Once the file is filled, add the @include directive with its name in the piawka file.
# Example is given on the basis of theta_low because it features all three functions. 

@namespace "calc"

# stats::add_stat() should be called in the BEGIN block.
# Args to stats::add_stat should be ordered as:
#   - stat name as to be printed (lowercase by convention),
#   - extended description to be shown in `piawka list` output,
#   - "0" if the stat is within group and any nonzero value if between-groups,
#   - comma-separated spaceless list of dependencies stats (= those whose numerators/denominators are used) 
#
# Args can be omitted or assigned empty string "" to skip them.
#
BEGIN{
  stats::add_stat("theta_low", "Theta estimator based on sites with 0<allele_freq<0.33", 0, "lines")
}

# initiate_* is run once after parsing command line arguments and before reading in the VCF file.
# it is not mandatory to define.
# variables available within initiate_* are:
#   - arg::args[] -- array whose indices are the long versions of options passed with the command,
#                    values are either the option value or "1" if the option is a flag
#   - groups[] -- array whose indices are groups names; values are number of samples in group
#   - pairs[][] -- nested array, pairs of indices stand for pairs of groups used for 
#                  calculation of pairwise stats
#   - ploidy[] -- array whose indices are groups names; values are integers describing ploidy
#   - stats::within[]  -- indexed 1 to stats::n_within in order of @include appearance in `piawka` file, 
#                       values are stats names;
#   - stats::between[]  -- indexed 1 to stats::n_between in order of @include appearance in `piawka` file,
#                        values are stats names;
#   - stats::used[]  -- indexed 1 to stats::n_used == stats::n_within + stats::n_between, 
#                     union of stats::within and stats::between, values are stats names;
#   - stats::printed[] -- indices are stats names to be printed in the output contains printed stats, 
#                         includes helper stats if arg::args["dependencies"]==1

function initiate_theta_low(){
  say("Warning: -s theta_low concerns derived allele frequencies, make sure the VCF is polarized correctly!")
  piawka::assert( arg::args["persite"]!=1, "-s theta_low at a single site is meaningless" )
  piawka::assert( arg::args["mult"]!=1, "-s theta_low is not defined in multiallelic comparisons" )
}

# increment_* is run once per VCF line that is used for the stats calculation.
# it should take one argument for group id if it is within-group (i by convention)
# and two group id args otherwise (i,j by convention).
# The function should assign the values to the following variables:
#   - thisnum[i][*] or thisnum[i,j][*] where * is the statistic name:
#       numerator for the statistic for this group(-s) at this site;
#   - thisden[i][*] or thisden[i,j][*] where * is the statistic name:
#       denominator for the statistic for this group(-s) at this site.
# Unless finalize_* functions are specified, per-site values of thisnum as well as thisden 
# are summed per group (pair) and per statistic and returned as the stat value for the region.
#
# Following variables are available when calculating statistics:
#   - all of what is available for initiate_*;
#   - n[] -- indices = group ids, values = # genotypes for this group at this site;
#   - a[][] -- first indices = group ids, second indices = allele ids, 
#              values = genotype counts with this allele;
#   - a_ij[] -- for pairwise stats, union of a[i] and a[j], indices = allele ids,
#               values = genotype counts for this allele across both groups;
#   - miss[] -- indices = group ids, values = # missing genotype calls for this group & site;
#   - all variables from increment_* functions of dependencies stats (include dependencies first!);
#   - for between-group stats, all variables defined by increment_* of WITHIN-GROUP stats are available.
#
# Feel free to declare other variables in the function, 
# then add them to args list to restrict their scope to the function itself

function increment_theta_low(i,    alt, weight){
  if ( length(a[i]) > 1 && n[i] > 3 && 0 in a[i] ) {
    if (a[i][0] > (2/3)*n[i]) {
      weight=1/int(n[i]/3)
      alt=n[i]-a[i][0]
      thisnum[i]["theta_low"]= alt * weight
    }
  }
  # assigning values to thisden[i]["theta_low"] is omitted because
  # denominator is defined in the finalize_theta_low() function
}

# The finalize_* function, if defined, is called once before printing the results.
# It can be used to alter the values of num[][] and den[][]
# that default to sums of thisnum[][] and thisden[][] respectively.
#
# Available variables:
#   - all of what is available for increment_*;
#   - all variables from finalize_* functions of dependencies stats (include dependencies first!);
#   - for between-group stats, all variables defined by finalize_* of WITHIN-GROUP stats are available.

function finalize_theta_low(i) {
  den[i]["theta_low"]=num[i]["lines"]

