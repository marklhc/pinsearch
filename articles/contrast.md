# Using Contrasts

``` r

library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(pinsearch)
```

This vignette demonstrates the use of contrasts in the
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md)
function based on a real data set.

## Load Data

``` r

data(lui_sim)
```

## Descriptive Statistics

``` r

summary(lui_sim)
#>      class1          class2          class3          class4     
#>  Min.   :1.000   Min.   :1.000   Min.   :1.000   Min.   :1.000  
#>  1st Qu.:2.000   1st Qu.:1.000   1st Qu.:1.000   1st Qu.:1.000  
#>  Median :2.000   Median :2.000   Median :2.000   Median :1.000  
#>  Mean   :1.779   Mean   :1.577   Mean   :1.703   Mean   :1.432  
#>  3rd Qu.:2.000   3rd Qu.:2.000   3rd Qu.:2.000   3rd Qu.:2.000  
#>  Max.   :2.000   Max.   :2.000   Max.   :2.000   Max.   :2.000  
#>      class5          class6          class7          class8         class9     
#>  Min.   :1.000   Min.   :1.000   Min.   :1.000   Min.   :1.00   Min.   :1.000  
#>  1st Qu.:1.000   1st Qu.:1.000   1st Qu.:1.000   1st Qu.:1.00   1st Qu.:1.000  
#>  Median :1.000   Median :1.000   Median :2.000   Median :2.00   Median :2.000  
#>  Mean   :1.379   Mean   :1.223   Mean   :1.655   Mean   :1.58   Mean   :1.638  
#>  3rd Qu.:2.000   3rd Qu.:1.000   3rd Qu.:2.000   3rd Qu.:2.00   3rd Qu.:2.000  
#>  Max.   :2.000   Max.   :2.000   Max.   :2.000   Max.   :2.00   Max.   :2.000  
#>     class10         class11         class12         class13     
#>  Min.   :1.000   Min.   :1.000   Min.   :1.000   Min.   :1.000  
#>  1st Qu.:1.000   1st Qu.:1.000   1st Qu.:1.000   1st Qu.:1.000  
#>  Median :2.000   Median :2.000   Median :2.000   Median :1.000  
#>  Mean   :1.513   Mean   :1.722   Mean   :1.653   Mean   :1.356  
#>  3rd Qu.:2.000   3rd Qu.:2.000   3rd Qu.:2.000   3rd Qu.:2.000  
#>  Max.   :2.000   Max.   :2.000   Max.   :2.000   Max.   :2.000  
#>     class14         class15          group      
#>  Min.   :1.000   Min.   :1.000   Min.   :1.000  
#>  1st Qu.:1.000   1st Qu.:1.000   1st Qu.:2.000  
#>  Median :2.000   Median :1.000   Median :4.000  
#>  Mean   :1.742   Mean   :1.436   Mean   :3.615  
#>  3rd Qu.:2.000   3rd Qu.:2.000   3rd Qu.:5.000  
#>  Max.   :2.000   Max.   :2.000   Max.   :6.000
```

## Specification search

### Configural invariance

``` r

config_mod <- "
f1 =~ class1 + class2 + class3 + class4 + class5 + class6 + 
      class7 + class8 + class9 + class10 + class11 + class12 + 
      class13 + class14 + class15
"
config_fit <- cfa(config_mod, data = lui_sim, group = "group", ordered = TRUE)
# Release covariance between class1 and class2
config_fit2 <- update(config_fit, c(config_mod, "class1 ~~ class2"))
anova(config_fit, config_fit2)
#> 
#> Scaled and Shifted Chi-Squared Difference Test (method = "satorra.2000")
#> 
#> lavaan->lavTestLRT():  
#>    lavaan NOTE: The "Chisq" column contains standard test statistics, not the 
#>    robust test that should be reported per model. A robust difference test is 
#>    a function of two standard (not robust) statistics.
#> 
#>              Df AIC BIC  Chisq Chisq diff  RMSEA Df diff Pr(>Chisq)    
#> config_fit2 534         354.80                                         
#> config_fit  540         411.35     36.784 0.2357       6   1.94e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Forward search

``` r

knitr::kable(ps[[2]])
```

| group | lhs     | rhs     | type       |
|------:|:--------|:--------|:-----------|
|     6 | f1      | class3  | loadings   |
|     6 | f1      | class6  | loadings   |
|     2 | f1      | class15 | loadings   |
|     6 | f1      | class8  | loadings   |
|     1 | f1      | class13 | loadings   |
|     2 | f1      | class11 | loadings   |
|     4 | f1      | class7  | loadings   |
|     4 | f1      | class4  | loadings   |
|     5 | f1      | class5  | loadings   |
|     2 | f1      | class14 | loadings   |
|     1 | f1      | class3  | loadings   |
|     3 | f1      | class9  | loadings   |
|     2 | f1      | class5  | loadings   |
|     5 | f1      | class6  | loadings   |
|     5 | f1      | class15 | loadings   |
|     4 | f1      | class8  | loadings   |
|     1 | f1      | class12 | loadings   |
|     2 | f1      | class8  | loadings   |
|     4 | class14 | t1      | thresholds |
|     3 | class3  | t1      | thresholds |
|     2 | class9  | t1      | thresholds |
|     5 | class8  | t1      | thresholds |
|     2 | class10 | t1      | thresholds |
|     2 | class12 | t1      | thresholds |
|     1 | class10 | t1      | thresholds |
|     4 | class1  | t1      | thresholds |
|     1 | class15 | t1      | thresholds |
|     4 | class2  | t1      | thresholds |
|     3 | class8  | t1      | thresholds |

## Effect Size

``` r

# Obtain fmacs effect size
(f_omni <- pin_effsize(ps[[1]]))
#>       class1-f1 class2-f1 class3-f1  class4-f1  class5-f1 class6-f1  class7-f1
#> fmacs 0.1199234 0.1105238 0.2381024 0.09092262 0.04957459 0.1858849 0.09464019
#>       class8-f1 class9-f1 class10-f1 class11-f1 class12-f1 class13-f1
#> fmacs 0.1744948 0.1149935   0.126812  0.0353763  0.1087859  0.0708557
#>       class14-f1 class15-f1
#> fmacs  0.1860244   0.142588
# fmacs by gender
(f_gender <- pin_effsize(ps[[1]], group_factor = c(1, 1, 1, 2, 2, 2)))
#>        class1-f1  class2-f1 class3-f1  class4-f1  class5-f1  class6-f1
#> fmacs 0.03579041 0.03298516 0.1443592 0.02713531 0.03031336 0.04858438
#>       class7-f1  class8-f1 class9-f1 class10-f1 class11-f1 class12-f1
#> fmacs 0.0282448 0.03521963  0.083309 0.04836829 0.01813981 0.01888187
#>       class13-f1 class14-f1 class15-f1
#> fmacs 0.02702558 0.03136801 0.04376062
# fmacs by ethnicity
(f_eth <- pin_effsize(ps[[1]], group_factor = c(1, 2, 3, 1, 2, 3)))
#>       class1-f1  class2-f1 class3-f1  class4-f1  class5-f1 class6-f1  class7-f1
#> fmacs 0.0580995 0.05354566 0.1249914 0.04404945 0.01246479 0.1131955 0.04585051
#>        class8-f1  class9-f1 class10-f1 class11-f1 class12-f1 class13-f1
#> fmacs 0.08975193 0.07570749 0.07851749 0.02558737 0.08675476 0.04387133
#>       class14-f1 class15-f1
#> fmacs 0.07596325  0.1148166
# interaction (using contrast matrix)
contr <- local({
    gen <- factor(c("F", "M"))
    contrasts(gen) <- contr.sum(length(gen))
    eth <- factor(1:3)
    contrasts(eth) <- contr.sum(length(eth))
    model.matrix(~ gen * eth, data = expand.grid(eth = eth, gen = gen))
})
(f_int <- pin_effsize(ps[[1]], contrast = contr[, 5:6, drop = FALSE]))
#>       class1-f1  class2-f1 class3-f1  class4-f1  class5-f1 class6-f1  class7-f1
#> fmacs 0.0580995 0.05354566 0.1641827 0.04404945 0.04276027 0.1131955 0.04585051
#>       class8-f1  class9-f1 class10-f1 class11-f1 class12-f1 class13-f1
#> fmacs 0.1310042 0.07570749 0.07851749 0.02558737 0.08675476 0.04387133
#>       class14-f1 class15-f1
#> fmacs  0.1208021 0.07983613
```

``` r

# Render as table
item_names <- gsub("class", replacement = "Item ", colnames(f_omni)) |>
    gsub(pattern = "-f1", replacement = "")
t(rbind(f_omni, f_gender, f_eth, f_int)) |>
    as.data.frame() |>
    `dimnames<-`(list(
        item_names,
        c("Overall", "Gender", "Ethnicity", "Gender x Ethnicity")
    )) |>
    knitr::kable(digits = 2, caption = "$f_\\text{MACS}$ effect sizes")
```

|         | Overall | Gender | Ethnicity | Gender x Ethnicity |
|:--------|--------:|-------:|----------:|-------------------:|
| Item 1  |    0.12 |   0.04 |      0.06 |               0.06 |
| Item 2  |    0.11 |   0.03 |      0.05 |               0.05 |
| Item 3  |    0.24 |   0.14 |      0.12 |               0.16 |
| Item 4  |    0.09 |   0.03 |      0.04 |               0.04 |
| Item 5  |    0.05 |   0.03 |      0.01 |               0.04 |
| Item 6  |    0.19 |   0.05 |      0.11 |               0.11 |
| Item 7  |    0.09 |   0.03 |      0.05 |               0.05 |
| Item 8  |    0.17 |   0.04 |      0.09 |               0.13 |
| Item 9  |    0.11 |   0.08 |      0.08 |               0.08 |
| Item 10 |    0.13 |   0.05 |      0.08 |               0.08 |
| Item 11 |    0.04 |   0.02 |      0.03 |               0.03 |
| Item 12 |    0.11 |   0.02 |      0.09 |               0.09 |
| Item 13 |    0.07 |   0.03 |      0.04 |               0.04 |
| Item 14 |    0.19 |   0.03 |      0.08 |               0.12 |
| Item 15 |    0.14 |   0.04 |      0.11 |               0.08 |

$`f_\text{MACS}`$ effect sizes {.table}
