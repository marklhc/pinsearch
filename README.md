
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pinsearch

<!-- badges: start -->
<!-- badges: end -->

The goal of pinsearch is to automate the process of performing
specification search in identifying noninvariant items to arrive at a
partial factorial invariance model.

## Installation

<!-- You can install the released version of pinsearch from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("pinsearch") -->
<!-- ``` -->
<!-- And the development version from [GitHub](https://github.com/) with: -->

You can install the development version of pinsearch from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("marklhc/pinsearch")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(pinsearch)
library(lavaan)
#> This is lavaan 0.6-9
#> lavaan is FREE software! Please report any bugs.
HS.model <- '  visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
# Output the final partial invariance model, and the noninvariant items
pinSearch(HS.model, data = HolzingerSwineford1939, 
          group = "school", type = "intercepts")
#> $`Partial Invariance Fit`
#> lavaan 0.6-9 ended normally after 69 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        66
#>   Number of equality constraints                    16
#>                                                       
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#>                                                       
#> Model Test User Model:
#>                                                       
#>   Test statistic                               129.422
#>   Degrees of freedom                                58
#>   P-value (Chi-square)                           0.000
#>   Test statistic for each group:
#>     Pasteur                                     71.170
#>     Grant-White                                 58.253
#> 
#> $`Non-Invariant Items`
#>   lhs rhs group       type
#> 1  x3         1 intercepts
#> 2  x7         1 intercepts
# Compute dmacs effect size (added in version 0.1.2)
pinSearch(HS.model, data = HolzingerSwineford1939, 
          group = "school", type = "intercepts",
          effect_size = TRUE)
#> $`Partial Invariance Fit`
#> lavaan 0.6-9 ended normally after 69 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        66
#>   Number of equality constraints                    16
#>                                                       
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#>                                                       
#> Model Test User Model:
#>                                                       
#>   Test statistic                               129.422
#>   Degrees of freedom                                58
#>   P-value (Chi-square)                           0.000
#>   Test statistic for each group:
#>     Pasteur                                     71.170
#>     Grant-White                                 58.253
#> 
#> $`Non-Invariant Items`
#>   lhs rhs group       type
#> 1  x3         1 intercepts
#> 2  x7         1 intercepts
#> 
#> $dmacs
#>                        x3-visual x7-visual x3-textual x7-textual  x3-speed
#> Pasteur vs Grant-White 0.4824515 0.4161323  0.4824515  0.4161323 0.4824515
#>                         x7-speed
#> Pasteur vs Grant-White 0.4161323
```
