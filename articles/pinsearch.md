# (P)artial (in)variance specification search

This package implements the sequential specification search method for
identifying parameters violating measurement invariance across groups.
The method is discussed in, for example, [Yoon and Millsap
(2007)](https://doi.org/10.1080/10705510701301677), and is a popular
method for identifying a partial invariance model.

## Workflow

A suggested workflow for using the main function
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md)
is

1.  Define the [`lavaan`](https://www.lavaan.ugent.be/) syntax for the
    configural invariance model.
2.  Conduct a
    [`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md)
    with specific test statistic and stopping criteria. The default is
    to use modification indices as in Yoon and Millsap (2007), but the
    score test is also supported. The function returns a final partial
    invariance model and a list of noninvariant parameters in the search
    order.
3.  Obtain effect size statistics for the noninvariant items.

## Example

Here we use the classic data set by Holzinger & Swineford, with the goal
of searching for a partial strict invariance model.

``` r

library(pinsearch)
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
mod <-
    " visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9 "
# Output the final partial invariance model, and the noninvariant items
ps1 <-
    pinSearch(mod,
        data = HolzingerSwineford1939,
        group = "school", type = "residuals",
        inv_test = "mod"  # use modification indices for search (default)
    )
```

Here are the noninvariant parameters:

``` r

ps1$`Non-Invariant Items`
#>   group lhs rhs       type
#> 1     1  x3     intercepts
#> 2     1  x7     intercepts
```

And the final partial invariance model

``` r

ps1$`Partial Invariance Fit` |>
    summary()
#> lavaan 0.7-2 ended normally after 64 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        66
#>   Number of equality constraints                    25
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                               147.260
#>   Degrees of freedom                                67
#>   P-value (Chi-square)                           0.000
#>   Test statistic for each group:
#>     Pasteur                                     79.438
#>     Grant-White                                 67.823
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Expected
#>   Information saturated (h1) model          Structured
#> 
#> 
#> Group 1 [Pasteur]:
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   visual =~                                           
#>     x1      (.p1.)    0.881    0.093    9.474    0.000
#>     x2      (.p2.)    0.549    0.085    6.437    0.000
#>     x3      (.p3.)    0.726    0.084    8.599    0.000
#>   textual =~                                          
#>     x4      (.p4.)    0.945    0.069   13.652    0.000
#>     x5      (.p5.)    1.064    0.077   13.832    0.000
#>     x6      (.p6.)    0.882    0.065   13.554    0.000
#>   speed =~                                            
#>     x7      (.p7.)    0.564    0.071    7.981    0.000
#>     x8      (.p8.)    0.674    0.073    9.189    0.000
#>     x9      (.p9.)    0.603    0.070    8.666    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   visual ~~                                           
#>     textual           0.436    0.089    4.908    0.000
#>     speed             0.332    0.110    3.012    0.003
#>   textual ~~                                          
#>     speed             0.315    0.097    3.240    0.001
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1      (.25.)    4.910    0.093   52.673    0.000
#>    .x2      (.26.)    6.072    0.079   76.983    0.000
#>    .x3                2.487    0.090   27.772    0.000
#>    .x4      (.28.)    2.784    0.086   32.195    0.000
#>    .x5      (.29.)    4.029    0.096   41.811    0.000
#>    .x6      (.30.)    1.927    0.081   23.746    0.000
#>    .x7                4.432    0.082   53.865    0.000
#>    .x8      (.32.)    5.568    0.073   75.821    0.000
#>    .x9      (.33.)    5.411    0.071   76.698    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1      (.10.)    0.636    0.101    6.303    0.000
#>    .x2      (.11.)    1.101    0.100   10.973    0.000
#>    .x3      (.12.)    0.724    0.085    8.529    0.000
#>    .x4      (.13.)    0.383    0.047    8.096    0.000
#>    .x5      (.14.)    0.435    0.057    7.613    0.000
#>    .x6      (.15.)    0.353    0.042    8.336    0.000
#>    .x7      (.16.)    0.738    0.076    9.696    0.000
#>    .x8      (.17.)    0.479    0.073    6.592    0.000
#>    .x9      (.18.)    0.580    0.069    8.386    0.000
#>     visual            1.000                           
#>     textual           1.000                           
#>     speed             1.000                           
#> 
#> 
#> Group 2 [Grant-White]:
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   visual =~                                           
#>     x1      (.p1.)    0.881    0.093    9.474    0.000
#>     x2      (.p2.)    0.549    0.085    6.437    0.000
#>     x3      (.p3.)    0.726    0.084    8.599    0.000
#>   textual =~                                          
#>     x4      (.p4.)    0.945    0.069   13.652    0.000
#>     x5      (.p5.)    1.064    0.077   13.832    0.000
#>     x6      (.p6.)    0.882    0.065   13.554    0.000
#>   speed =~                                            
#>     x7      (.p7.)    0.564    0.071    7.981    0.000
#>     x8      (.p8.)    0.674    0.073    9.189    0.000
#>     x9      (.p9.)    0.603    0.070    8.666    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   visual ~~                                           
#>     textual           0.506    0.123    4.125    0.000
#>     speed             0.634    0.161    3.929    0.000
#>   textual ~~                                          
#>     speed             0.418    0.134    3.115    0.002
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1      (.25.)    4.910    0.093   52.673    0.000
#>    .x2      (.26.)    6.072    0.079   76.983    0.000
#>    .x3      (.27.)    1.951    0.114   17.044    0.000
#>    .x4      (.28.)    2.784    0.086   32.195    0.000
#>    .x5      (.29.)    4.029    0.096   41.811    0.000
#>    .x6      (.30.)    1.927    0.081   23.746    0.000
#>    .x7      (.31.)    3.992    0.099   40.135    0.000
#>    .x8      (.32.)    5.568    0.073   75.821    0.000
#>    .x9      (.33.)    5.411    0.071   76.698    0.000
#>     visual            0.062    0.146    0.423    0.672
#>     textual           0.608    0.129    4.722    0.000
#>     speed            -0.127    0.157   -0.807    0.420
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1      (.10.)    0.636    0.101    6.303    0.000
#>    .x2      (.11.)    1.101    0.100   10.973    0.000
#>    .x3      (.12.)    0.724    0.085    8.529    0.000
#>    .x4      (.13.)    0.383    0.047    8.096    0.000
#>    .x5      (.14.)    0.435    0.057    7.613    0.000
#>    .x6      (.15.)    0.353    0.042    8.336    0.000
#>    .x7      (.16.)    0.738    0.076    9.696    0.000
#>    .x8      (.17.)    0.479    0.073    6.592    0.000
#>    .x9      (.18.)    0.580    0.069    8.386    0.000
#>     visual            0.856    0.208    4.114    0.000
#>     textual           0.980    0.183    5.369    0.000
#>     speed             1.401    0.326    4.294    0.000
```

## Effect Size

The package supports the $`d_\text{MACS}`$ effect size statistics by
[Nye and Drasgow (2011)](https://doi.org/10.1037/a0022955) for
noninvariant items, and the generalization to an $`f_\text{bias}`$
statistic for more than two groups.

These can be obtained by setting the `effect_size = TRUE` argument:

``` r

pinSearch(mod,
    data = HolzingerSwineford1939,
    group = "school", type = "residuals",
    inv_test = "mod",
    effect_size = TRUE
)
#> $`Partial Invariance Fit`
#> lavaan 0.7-2 ended normally after 64 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        66
#>   Number of equality constraints                    25
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                               147.260
#>   Degrees of freedom                                67
#>   P-value (Chi-square)                           0.000
#>   Test statistic for each group:
#>     Pasteur                                     79.438
#>     Grant-White                                 67.823
#> 
#> $`Non-Invariant Items`
#>   group lhs rhs       type
#> 1     1  x3     intercepts
#> 2     1  x7     intercepts
#> 
#> $effect_size
#>       x3-visual  x7-speed
#> dmacs 0.4866004 0.4161171
```

Note that in the above output, the effect sizes are printed three times
for `x3` and `x7`, which have non-invariant intercepts. This is because
there are three latent variables in the model, and the function allows
for the possibility of cross-loadings. When each item is only loaded on
one latent variable, user can ignore the duplicated values (we will try
to fix this in future releases).

Another possibility is to use the
[`pin_effsize()`](https://marklhc.github.io/pinsearch/reference/es_lavaan.md)
function, which also works for partial invariance model not obtained
from
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md).

``` r

# Fit the partial strict model by hand
pstrict_fit <- cfa(mod, data = HolzingerSwineford1939,
                   group = "school",
                   group.equal = c("loadings", "intercepts",
                                   "residuals"),
                   group.partial = c("x3~1", "x7~1"))
# Effect sizes
pin_effsize(pstrict_fit)
#>       x3-visual  x7-speed
#> dmacs 0.4866011 0.4161184
```

## Cross-Loadings

The function also supports models with cross-loadings. For example,

``` r

mod2 <- "
   visual =~ x1 + x2 + x3 + x9
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"
ps2 <-
    pinSearch(mod2,
        data = HolzingerSwineford1939,
        group = "school", type = "residuals",
        inv_test = "mod",
        effect_size = TRUE
    )
ps2
#> $`Partial Invariance Fit`
#> lavaan 0.7-2 ended normally after 55 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        68
#>   Number of equality constraints                    26
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                               113.828
#>   Degrees of freedom                                66
#>   P-value (Chi-square)                           0.000
#>   Test statistic for each group:
#>     Pasteur                                     67.913
#>     Grant-White                                 45.915
#> 
#> $`Non-Invariant Items`
#>   group lhs rhs       type
#> 1     1  x3     intercepts
#> 2     1  x7     intercepts
#> 
#> $effect_size
#>       x3-visual  x7-speed
#> dmacs  0.472144 0.4069861
```

It just happened in this examplle that the final results are not
affected.

## Score Test

One can also use the score test instead of the modification indices, as
[starting in version 0.6-17, `lavaan` recommends using that for equality
constraints](https://lavaan.ugent.be/history/dot6.html#version-0.6-17).
This can be achieved with `inv_test = "score"`.

``` r

pinSearch(mod2,
    data = HolzingerSwineford1939,
    group = "school", type = "residuals",
    inv_test = "score",
    effect_size = TRUE
)
#> $`Partial Invariance Fit`
#> lavaan 0.7-2 ended normally after 59 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        68
#>   Number of equality constraints                    25
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                               109.770
#>   Degrees of freedom                                65
#>   P-value (Chi-square)                           0.000
#>   Test statistic for each group:
#>     Pasteur                                     66.273
#>     Grant-White                                 43.497
#> 
#> $`Non-Invariant Items`
#>   group lhs rhs       type
#> 1     2  x3     intercepts
#> 2     2  x7     intercepts
#> 3     2  x7  x7  residuals
#> 
#> $effect_size
#>       x3-visual  x7-speed
#> dmacs 0.4717683 0.4071137
```

The test is slightly more sensitive and identifies an additional
non-invariant uniqueness for item `x7`.

## Likelihood Ratio Test (LRT)

LRT is generally more accurate, but more computationally expensive.

``` r

pinSearch(mod2,
    data = HolzingerSwineford1939,
    group = "school", type = "residuals",
    inv_test = "lrt",
    effect_size = TRUE
)
#> $`Partial Invariance Fit`
#> lavaan 0.7-2 ended normally after 59 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        68
#>   Number of equality constraints                    25
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                               109.770
#>   Degrees of freedom                                65
#>   P-value (Chi-square)                           0.000
#>   Test statistic for each group:
#>     Pasteur                                     66.273
#>     Grant-White                                 43.497
#> 
#> $`Non-Invariant Items`
#>   group lhs rhs       type
#> 1     1  x3     intercepts
#> 2     1  x7     intercepts
#> 3     1  x7  x7  residuals
#> 
#> $effect_size
#>       x3-visual  x7-speed
#> dmacs 0.4717683 0.4071137
```

In this example, the results using LRT are the same as the score test.

## Controlling for False Discovery Rate

The function also supports adjusting for multiplicity by using a more
stringent threshold for flagging a parameter as noninvariant, using the
multistage false discovery rate procedure by [Benjamini & Gavrilov
(2009)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-3/issue-1/A-simple-forward-selection-procedure-based-on-false-discovery-rate/10.1214/08-AOAS194.pdf)
for forward selection. It is up for debate whether controlling for false
discovery rate is as important in measurement/factorial invariance
analyses as in other context, and it probably depends on whether one is
interested in identifying non-invariant items or in inferences for the
latent parameters. If this is of interest, one can use the
`control_fdr = TRUE` option.

``` r

pinSearch(mod2,
    data = HolzingerSwineford1939,
    group = "school", type = "residuals",
    inv_test = "score",
    effect_size = TRUE,
    control_fdr = TRUE
)
#> $`Partial Invariance Fit`
#> lavaan 0.7-2 ended normally after 55 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        68
#>   Number of equality constraints                    26
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                               113.828
#>   Degrees of freedom                                66
#>   P-value (Chi-square)                           0.000
#>   Test statistic for each group:
#>     Pasteur                                     67.913
#>     Grant-White                                 45.915
#> 
#> $`Non-Invariant Items`
#>   group lhs rhs       type
#> 1     2  x3     intercepts
#> 2     2  x7     intercepts
#> 
#> $effect_size
#>       x3-visual  x7-speed
#> dmacs  0.472144 0.4069861
```

The uniqueness for `x7` is no longer flagged.
