
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pinvsearch

<!-- badges: start -->

<!-- badges: end -->

The goal of pinvsearch is to automate the process of performing
specification search in identifying noninvariant items to arrive at a
partial factorial invariance model.

## Installation

<!-- You can install the released version of pinvsearch from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("pinvsearch") -->

<!-- ``` -->

<!-- And the development version from [GitHub](https://github.com/) with: -->

You can install the development version of pinvsearch from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("marklhc/pinvsearch")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(pinvsearch)
library(lavaan)
#> This is lavaan 0.6-5
#> lavaan is BETA software! Please report any bugs.
HS.model <- '  visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
# Output the final partial invariance model, and the noninvariant items
pinvSearch(HS.model, data = HolzingerSwineford1939, 
           group = "school", type = "intercepts")
#> $`Partial Invariance Fit`
#> lavaan 0.6-5 ended normally after 69 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of free parameters                         66
#>   Number of equality constraints                    16
#>   Row rank of the constraints matrix                16
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
#>   item group       type
#> 1   x3     1 intercepts
#> 2   x7     1 intercepts
```
