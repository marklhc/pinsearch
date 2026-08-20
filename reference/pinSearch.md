# Search for noninvariant parameters across groups.

The function implements the sequential selection method similar to that
discussed in [Yoon and Millsap
(2007)](https://doi.org/10.1080/10705510701301677). The function
proceeds in the order of metric, scalar (threshold), and strict
invariance. In each stage, invariance constraints in all items, and the
constraint associated with the biggest test statistic above a predefined
threshold is freed, before recomputing the test statistic for the next
constraint to free.

## Usage

``` r
pinSearch(
  config_mod,
  ...,
  type = c("loadings", "intercepts", "thresholds", "residuals", "residual.covariances"),
  inv_test = c("mod", "score", "lrt"),
  sig_level = 0.05,
  control_fdr = FALSE,
  min2 = FALSE,
  effect_size = FALSE,
  progress = FALSE
)
```

## Arguments

- config_mod:

  Syntax of a configural invariance model to be passed to
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html).

- ...:

  Additonal arguments passed to
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html).

- type:

  Character variable indicating the stage of invariance to be searched.
  Currently supported options are (a) for continuous indicators,
  "loadings", "intercepts", "residuals", and "residual.covariances",
  and (b) "loadings", "thresholds", "residual.covariances", in an
  increasingly strict order. A stricter model (e.g.,
  "residual.covariances") will have constraints of all previous stages.

- inv_test:

  Character variable indicating the statistical test to be used for
  specification search. Currently supported options are `"mod"` for
  modification index using
  [`lavaan::modindices()`](https://rdrr.io/pkg/lavaan/man/modificationIndices.html),
  `"score"` for score test statistic using
  [`lavaan::lavTestScore()`](https://rdrr.io/pkg/lavaan/man/lavTestScore.html),
  and (experimental) `"lrt"` for likelihood ratio test statistic using
  [`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html).

- sig_level:

  Significance level used to determine whether the parameter associated
  with the highest modification index should be removed. Default is .05.

- control_fdr:

  Logical; whether to use adjust for false discovery rate for multiple
  testing. If `TRUE`, the method by Benjamini & Gavrilov (2009) will be
  used.

- min2:

  Logical; whether to require at least two invariant items when
  searching for noninvariance.

- effect_size:

  Logical; whether to compute dmacs (two groups) or fmacs (\> two
  groups) effect size or not (default). This is an experimental feature.

- progress:

  Logical; an experimental feature of showing a progress bar if `TRUE`.
  Because the number of steps is unknown until the stopping criteria are
  reached, the progress bar may be inaccurate.

## Value

A list of three elements:

- `Partial Invariance Fit`: A
  [`lavaan::lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan-class.html)
  object containing the final partial invariance model.

- `Non-Invariant Items`: A data frame of non-invariant parameters.

- `effect_size`: Effect size statistics obtained from
  [`pin_effsize()`](https://marklhc.github.io/pinsearch/reference/es_lavaan.md).

## Details

Note that when an item has a non-invariant loading, the corresponding
intercept constraint will automatically be freed, as intercept
difference across groups is sensitive to the location of the zero point
for the latent variable and the item.

For a particular stage of invariance constraints, the Benjamini &
Gavrilov method uses an adjusted \\alpha\\ level of \$\$iq / \[m + 1 -
i(1 - q)\]\$\$ where \\i\\ is the step index in the search, \\m\\ is the
maximum number of constraints that can be freed, and \\q\\ is the
desirable significance level.

Under lavaan \>= 0.7, ordinal (categorical) models in which thresholds
are tied across groups can yield a degenerate stage test (a likelihood
ratio test with `Df diff = 0`). `pinSearch()` identifies the threshold
stage (and any later stage, but not the configural or loadings stages)
with the item intercepts fixed at 0 and the latent mean fixed in the
first group only (free in the remaining groups), so threshold invariance
remains testable, including for 2-group binary data. Any other
degenerate stage test falls back to the per-parameter test, so results
for ordinary multi-category and 3-group ordinal data are unchanged;
continuous models are not affected.

## References

Yoon, M., & Millsap, R. E. (2007). Detecting violations of factorial
invariance using data-based specification searches: A Monte Carlo study.
Structural Equation Modeling: A Multidisciplinary Journal, 14(3),
435-463.

## Examples

``` r
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(MASS)
# Simulate random data
set.seed(14) # for reproducible results
mod_ninv <- "f =~ c(1, 0.8, .6)*y1 + c(0.8, 1.2, 1.0)*y2"
mod_inv <- paste0("1*y", 3:6, collapse = " + ")
sim_mod <- paste(
    paste(mod_ninv, "+", mod_inv),
    "f ~~ c(1, 1.3, 1.5)*f
     f ~ c(0, 0.5, 1.0)*1
     y1 ~ c(0, 0.3, 0)*1
     y3 ~ c(0.3, 0, -0.3)*1
     y1 ~~ c(1, .5, 1)*y1",
    sep = "\n"
)
# The uniqueness of each item is assumed to be 1.0
dat_sim <- simulateData(sim_mod, sample.nobs = c(100, 100, 100))
# Fit configural model:
sam_config <- paste(
    paste0("f =~ ", paste0("y", 1:6, collapse = " + "))
)
pinSearch(sam_config,
    data = dat_sim, group = "group",
    type = "intercepts"
)
#> $`Partial Invariance Fit`
#> lavaan 0.7-2 ended normally after 60 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        58
#>   Number of equality constraints                    19
#> 
#>   Number of observations per group:                   
#>     1                                              100
#>     2                                              100
#>     3                                              100
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                47.688
#>   Degrees of freedom                                42
#>   P-value (Chi-square)                           0.252
#>   Test statistic for each group:
#>     1                                           12.480
#>     2                                           20.956
#>     3                                           14.252
#> 
#> $`Non-Invariant Items`
#>   group lhs rhs       type
#> 1     1   f  y1   loadings
#> 2     1  y3     intercepts
#> 3     3  y1     intercepts
#> 4     2  y3     intercepts
#> 
```
