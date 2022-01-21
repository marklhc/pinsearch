
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

## Example 2

A simulated example with ordinal data

``` r
# Simulate data
library(MASS)
num_obs <- 500
lambda1 <- seq(.9, .6, length.out = 7)
lambda2 <- c(lambda1[1], 1, lambda1[3:7])
cov1 <- tcrossprod(lambda1) + diag(.5, 7)
dimnames(cov1) <- list(paste0("yy", 1:7), paste0("yy", 1:7))
thres1 <- rbind(seq(-1.5, 1.5, length.out = 7))
thres2 <- rbind(c(thres1[1], 0.25, thres1[3:6], 0.3))
mean1 <- rep(0, 7)
ystar1 <- mvrnorm(num_obs, mu = mean1, Sigma = cov1)
y1 <- ystar1
cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 7)
dimnames(cov2) <- dimnames(cov1)
mean2 <- lambda1 * .4
ystar2 <- mvrnorm(num_obs, mu = mean2, Sigma = cov2)
y2 <- ystar2
# Ordinal indicators
thres1 <- rbind(seq(-1.5, 0, length.out = 7),
                seq(-0.5, 0.25, length.out = 7),
                rep(1, 7))
thres2 <- rbind(c(thres1[1, 1], -0.5, thres1[1, 3:6], -0.5),
                thres1[2,],
                c(rep(1, 3), rep(0.5, 2), rep(1, 2)))
for (j in seq_len(ncol(ystar1))) {
    y1[, j] <- findInterval(ystar1[, j], thres1[, j])
}
for (j in seq_len(ncol(ystar2))) {
    y2[, j] <- findInterval(ystar2[, j], thres2[, j])
}
df <- rbind(cbind(y1, group = 1), cbind(y2, group = 2))
```

``` r
pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
          data = df, group = "group", type = "thresholds",
          ordered = paste0("yy", 1:7),
          effect_size = TRUE)
#> Unique variances are constrained to 1 for identification
#> $`Partial Invariance Fit`
#> lavaan 0.6-9 ended normally after 43 iterations
#> 
#>   Estimator                                       DWLS
#>   Optimization method                           NLMINB
#>   Number of model parameters                        58
#>   Number of equality constraints                    24
#>                                                       
#>   Number of observations per group:                   
#>     1                                              500
#>     2                                              500
#>                                                       
#> Model Test User Model:
#>                                               Standard      Robust
#>   Test Statistic                                21.980      32.588
#>   Degrees of freedom                                50          50
#>   P-value (Chi-square)                           1.000       0.973
#>   Scaling correction factor                                  0.605
#>   Shift parameter for each group:                                 
#>       1                                                     -1.875
#>       2                                                     -1.875
#>        simple second-order correction                             
#>   Test statistic for each group:
#>     1                                           12.160      18.228
#>     2                                            9.820      14.360
#> 
#> $`Non-Invariant Items`
#>   lhs rhs group       type
#> 1 yy2  t1     1 thresholds
#> 2 yy7  t1     2 thresholds
#> 3 yy4  t3     1 thresholds
#> 4 yy5  t3     1 thresholds
#> 
#> $dmacs
#>            yy2-f    yy4-f     yy5-f     yy7-f
#> 1 vs 2 0.2446916 0.191393 0.1395843 0.1999788
```
