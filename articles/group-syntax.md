# Group-Specific Model Syntax

In a factorial invariance analysis, the groups sometimes do not share
the same set of observed variables. A common instance is a shorter or
longer scale form: one group answers four items and another answers
five. A single shared configural model cannot be written for such data,
because it would reference an item that does not exist in every group.

`lavaan` handles this with **group-specific model syntax**: a `group:`
block defines a separate model for each group.
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md)
supports this syntax, so a partial invariance specification search can
still be run when the item sets differ across groups.

## The syntax

Write the configural model as a string (or a character vector) with one
`group: N` block per group; each `group: N` starts on its own line:

``` r

mod <- c(
  "group: 1", "F =~ y1 + y2 + y3 + y4",         # e.g. the short form
  "group: 2", "F =~ y1 + y2 + y3 + y4 + y5")    # e.g. the long form
```

The number of `group:` blocks, and their order, must match the groups
used to split the data. As usual,
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md)
passes `config_mod` straight through to
\[[`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html)\], so you
still pass `group = "..."` (via `...`) to split the data; `group: N`
then refers to the **Nth level** of that grouping factor.

## Example

We simulate a single five-item trait. Item `y5` is answered only by the
“long” form, and item `y3` loads more weakly in the “long” form, so the
search has something to find.

``` r

library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(pinsearch)
library(MASS, include.only = "mvrnorm")
set.seed(8)

lamS <- c(.90, .85, .80, .75)              # short form: y1-y4
lamL <- c(.90, .85, .50, .75, .70)         # long form:  y1-y5 (y3 weakened)
sigS <- tcrossprod(lamS) + diag(1 - lamS^2)
sigL <- tcrossprod(lamL) + diag(1 - lamL^2)

n  <- 500
gS <- mvrnorm(n, rep(0, 4), sigS)          # short form
gL <- mvrnorm(n, rep(0, 5), sigL)          # long form
dS <- as.data.frame(gS);  names(dS) <- paste0("y", 1:4)
dL <- as.data.frame(gL);  names(dL) <- paste0("y", 1:5)
df <- rbind(cbind(dS, y5 = NA, group = "short"),
            cbind(dL,          group = "long"))
df$group <- factor(df$group, levels = c("short", "long"))
```

A single model `F =~ y1 + y2 + y3 + y4 + y5` cannot describe both
groups, because `y5` is absent for the “short” form. With the
group-specific syntax, both forms are handled in one call:

``` r

mod <- c(
  "group: 1", "F =~ y1 + y2 + y3 + y4",          # short: no y5
  "group: 2", "F =~ y1 + y2 + y3 + y4 + y5")     # long:  full form
ps <- pinSearch(mod,
    data = df, group = "group",
    type = "intercepts"        # search loadings, then intercepts
)
ps$`Non-Invariant Items`
#>   group lhs rhs     type
#> 1     2   F  y3 loadings
```

The specification search flags the loading of `y3` in the long form. The
final partial invariance model is

``` r

summary(ps$`Partial Invariance Fit`)
#> lavaan 0.7-2 ended normally after 27 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        29
#>   Number of equality constraints                     6
#> 
#>   Number of observations per group:                   
#>     short                                          500
#>     long                                           500
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 5.945
#>   Degrees of freedom                                11
#>   P-value (Chi-square)                           0.877
#>   Test statistic for each group:
#>     short                                        3.474
#>     long                                         2.472
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Expected
#>   Information saturated (h1) model          Structured
#> 
#> 
#> Group 1 [short]:
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   F =~                                                
#>     y1      (.p1.)    0.924    0.035   26.247    0.000
#>     y2      (.p2.)    0.864    0.035   24.933    0.000
#>     y3      (.p3.)    0.805    0.038   21.131    0.000
#>     y4      (.p4.)    0.788    0.035   22.590    0.000
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     F                 0.000                           
#>    .y1      (.11.)   -0.055    0.045   -1.224    0.221
#>    .y2      (.12.)   -0.035    0.043   -0.800    0.424
#>    .y3      (.13.)   -0.088    0.045   -1.969    0.049
#>    .y4      (.14.)   -0.065    0.042   -1.559    0.119
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     F                 1.000                           
#>    .y1                0.207    0.023    9.140    0.000
#>    .y2                0.284    0.025   11.490    0.000
#>    .y3                0.350    0.028   12.698    0.000
#>    .y4                0.450    0.033   13.645    0.000
#> 
#> 
#> Group 2 [long]:
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   F =~                                                
#>     y1      (.p1.)    0.924    0.035   26.247    0.000
#>     y2      (.p2.)    0.864    0.035   24.933    0.000
#>     y3                0.482    0.044   10.919    0.000
#>     y4      (.p4.)    0.788    0.035   22.590    0.000
#>     y5                0.776    0.048   16.024    0.000
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     F                 0.054    0.067    0.808    0.419
#>    .y1      (.11.)   -0.055    0.045   -1.224    0.221
#>    .y2      (.12.)   -0.035    0.043   -0.800    0.424
#>    .y3               -0.070    0.045   -1.570    0.116
#>    .y4      (.14.)   -0.065    0.042   -1.559    0.119
#>    .y5               -0.092    0.052   -1.758    0.079
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     F                 0.995    0.101    9.889    0.000
#>    .y1                0.195    0.023    8.330    0.000
#>    .y2                0.341    0.029   11.934    0.000
#>    .y3                0.703    0.046   15.244    0.000
#>    .y4                0.432    0.032   13.357    0.000
#>    .y5                0.603    0.043   14.043    0.000
```

Effect sizes for the flagged item follow the usual workflow:

``` r

pinSearch(mod, data = df, group = "group",
          type = "intercepts", effect_size = TRUE)
#> $`Partial Invariance Fit`
#> lavaan 0.7-2 ended normally after 27 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        29
#>   Number of equality constraints                     6
#> 
#>   Number of observations per group:                   
#>     short                                          500
#>     long                                           500
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 5.945
#>   Degrees of freedom                                11
#>   P-value (Chi-square)                           0.877
#>   Test statistic for each group:
#>     short                                        3.474
#>     long                                         2.472
#> 
#> $`Non-Invariant Items`
#>   group lhs rhs     type
#> 1     2   F  y3 loadings
#> 
#> $effect_size
#>            y3-F
#> dmacs 0.3293811
```

## Related

### Keeping at least two invariant items

When a group has only a few items, the search could, in principle, free
so many loadings that fewer than two invariant indicators remained.
`min2 = TRUE` caps how many items may be freed during the search:

``` r

pinSearch(mod, data = df, group = "group",
          type = "loadings", min2 = TRUE)
```

### Ordered items

The same `group:` syntax works for ordered categorical indicators; just
add `ordered = ...` (and a `parameterization`, if needed) as you would
for a model without group blocks.

### Not the per-group-value idiom

Do not confuse this with lavaan’s `c(.)` per-group values (for example
`F =~ c(1.0, 0.9)*y1`, as in the
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md)
examples). That idiom still assumes the *same* variables in every group;
the `group:` block syntax is for *different* observed variables or
relationships.

## References

Yoon, M., & Millsap, R. E. (2007). Detecting violations of factorial
invariance using data-based specification searches: A Monte Carlo study.
*Structural Equation Modeling: A Multidisciplinary Journal, 14*(3),
435-463.
