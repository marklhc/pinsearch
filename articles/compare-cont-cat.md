# Comparing Effect Size for Continuous and Categorical Data

In this brief vignette, we compare the noninvariance effect size for
treating Likert-scale items as either continuous or ordinal.

``` r

library(pinsearch)
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(MASS, include.only = "mvrnorm")
library(difR)
```

Simulate a five-point example

``` r

set.seed(2049)
num_obs <- 500
lambda1 <- seq(.9, .6, length.out = 7)
lambda2 <- c(lambda1[1], 1, lambda1[3:7])
cov1 <- tcrossprod(lambda1) + diag(.5, 7)
dimnames(cov1) <- list(paste0("yy", 1:7), paste0("yy", 1:7))
thres1 <- rbind(seq(-1.5, 0, length.out = 7),
                seq(-0.5, 0.25, length.out = 7),
                rep(1, 7),
                seq(2, 2.5, length.out = 7))
thres2 <- rbind(c(thres1[1, 1], -0.5, thres1[1, 3:6], -0.5),
                thres1[2,],
                c(rep(1, 3), rep(0.5, 2), rep(1, 2)),
                thres1[4,])
mean1 <- rep(0, 7)
cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 7)
dimnames(cov2) <- dimnames(cov1)
mean2 <- lambda1 * .4
# Parameters:
# Lambda
(lam <- rbind(lambda1, lambda2))
#>         [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> lambda1  0.9 0.85  0.8 0.75  0.7 0.65  0.6
#> lambda2  0.9 1.00  0.8 0.75  0.7 0.65  0.6
# Thresholds
(nu <- rbind(setNames(c(thres1), nm = rep(1:7, each = 4)),
             setNames(c(thres2), nm = rep(1:7, each = 4))))
#>         1    1 1 1     2      2 2        2  3     3 3        3     4      4   4
#> [1,] -1.5 -0.5 1 2 -1.25 -0.375 1 2.083333 -1 -0.25 1 2.166667 -0.75 -0.125 1.0
#> [2,] -1.5 -0.5 1 2 -0.50 -0.375 1 2.083333 -1 -0.25 1 2.166667 -0.75 -0.125 0.5
#>         4    5 5   5        5     6     6 6        6    7    7 7   7
#> [1,] 2.25 -0.5 0 1.0 2.333333 -0.25 0.125 1 2.416667  0.0 0.25 1 2.5
#> [2,] 2.25 -0.5 0 0.5 2.333333 -0.25 0.125 1 2.416667 -0.5 0.25 1 2.5
```

Population ES

``` r

(dpop <- dmacs_ordered(
  nu,
  loadings = lam,
  thetas = 1,
  link = "probit",
  pooled_item_sd = 1,
  latent_mean = 0,
  latent_sd = 1
))
#>       [,1]      [,2] [,3]     [,4]      [,5] [,6]      [,7]
#> dmacs    0 0.2555382    0 0.143817 0.1449711    0 0.1696835
```

``` r

ystar1 <- mvrnorm(num_obs, mu = mean1, Sigma = cov1)
y1 <- ystar1
for (j in seq_len(ncol(ystar1))) {
    y1[, j] <- findInterval(ystar1[, j], thres1[, j])
}
ystar2 <- mvrnorm(num_obs, mu = mean2, Sigma = cov2)
y2 <- ystar2
for (j in seq_len(ncol(ystar2))) {
    y2[, j] <- findInterval(ystar2[, j], thres2[, j])
}
df <- rbind(cbind(y1, group = 1), cbind(y2, group = 2))
```

``` r

contfit <- cfa(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
               data = df, group = "group",
               group.equal = c("loadings", "intercepts"),
               group.partial = c("yy2~1", "yy4~1",
                                 "yy5~1", "yy7~1"))
(dcont <- pin_effsize(contfit))
#>           yy2-f     yy4-f     yy5-f     yy7-f
#> dmacs 0.1515229 0.1010677 0.0850543 0.2321431

ps_cat <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                    data = df, group = "group", type = "thresholds",
                    ordered = paste0("yy", 1:7))
#> Warning: lavaan->lav_model_vcov():  
#>    The variance-covariance matrix of the estimated parameters (vcov) does not 
#>    appear to be positive definite! The smallest eigenvalue (= 4.614582e-17) 
#>    is close to zero. This may be a symptom that the model is not identified.
(dcat <- pin_effsize(ps_cat$`Partial Invariance Fit`))
#>           yy2-f     yy4-f     yy5-f     yy7-f
#> dmacs 0.1912212 0.1400285 0.1219134 0.2069972
```

``` r

data.frame(
  rbind(continuous = dcont[1, ],
        categorical = dcat[1, ],
        population = dpop[c(2, 4, 5, 7)])
)
#>                 yy2.f     yy4.f     yy5.f     yy7.f
#> continuous  0.1515229 0.1010677 0.0850543 0.2321431
#> categorical 0.1912212 0.1400285 0.1219134 0.2069972
#> population  0.2555382 0.1438170 0.1449711 0.1696835
```

Compared to Mantel-Haenszel (MH)

``` r

# Dichotomous data (0, 1, 2 recoded to 0; 3, 4 recoded to 1)
df2 <- df
df2[, 1:7] <- as.integer(df2[, 1:7] >= 3)
```

``` r

ps2_cat <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                     data = df2, group = "group", type = "thresholds",
                     ordered = paste0("yy", 1:7))
#> Warning: lavaan->lav_model_vcov():  
#>    The variance-covariance matrix of the estimated parameters (vcov) does not 
#>    appear to be positive definite! The smallest eigenvalue (= 6.895669e-17) 
#>    is close to zero. This may be a symptom that the model is not identified.
(dcat <- pin_effsize(ps2_cat$`Partial Invariance Fit`))
#>           yy4-f    yy5-f
#> dmacs 0.3463871 0.324468
```

``` r

difMH(df2, group = "group", focal.name = 1, purify = TRUE)
#> 
#> Detection of Differential Item Functioning using Mantel-Haenszel method 
#> with continuity correction and with item purification
#> 
#> Results based on asymptotic inference 
#>  
#> Convergence reached after 2 iterations
#> 
#> Matching variable: test score 
#>  
#> No set of anchor items was provided 
#>  
#> No p-value adjustment for multiple comparisons 
#>  
#> Mantel-Haenszel Chi-square statistic: 
#>  
#>     Stat.   P-value    
#> yy1  0.0005  0.9815    
#> yy2  0.0530  0.8179    
#> yy3  0.2962  0.5863    
#> yy4 25.9472  0.0000 ***
#> yy5 18.0483  0.0000 ***
#> yy6  1.2175  0.2699    
#> yy7  0.0211  0.8844    
#> 
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1  
#> 
#> Detection threshold: 3.8415 (significance level: 0.05)
#> 
#> Items detected as DIF items: 
#>     
#>  yy4
#>  yy5
#> 
#>  
#> Effect size (ETS Delta scale): 
#>  
#> Effect size code: 
#>  'A': negligible effect 
#>  'B': moderate effect 
#>  'C': large effect 
#>  
#>     alphaMH deltaMH  
#> yy1  0.9806  0.0461 A
#> yy2  1.0707 -0.1606 A
#> yy3  1.1541 -0.3368 A
#> yy4  2.8198 -2.4362 C
#> yy5  2.4575 -2.1130 C
#> yy6  0.7535  0.6652 A
#> yy7  1.0566 -0.1294 A
#> 
#> Effect size codes: 0 'A' 1.0 'B' 1.5 'C' 
#>  (for absolute values of 'deltaMH') 
#>  
#> Output was not captured!
```

The results are consistent.
