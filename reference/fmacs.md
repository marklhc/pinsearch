# Compute fMACS effect size for two or more groups.

`fmacs` returns the fMACS effect size statistics given a set of loadings
and intercepts.

## Usage

``` r
fmacs(
  intercepts,
  loadings = NULL,
  pooled_item_sd,
  num_obs = NULL,
  weights = 0 * intercepts + 1,
  group_factor = NULL,
  contrast = contr.sum(nrow(intercepts)),
  latent_mean = 0,
  latent_sd = 1,
  item_weights = NULL
)

fmacs_ordered(
  thresholds,
  loadings,
  thetas = 1,
  num_obs = NULL,
  weights = 0 * loadings + 1,
  group_factor = NULL,
  contrast = contr.sum(nrow(thresholds)),
  link = c("probit", "logit"),
  pooled_item_sd = NULL,
  latent_mean = 0,
  latent_sd = 1,
  item_weights = NULL
)
```

## Arguments

- intercepts:

  A \\G \times p\\ matrix of measurement intercepts.

- loadings:

  A \\G \times p\\ matrix of factor loadings, where p is the number of
  items and G is the number of groups.

- pooled_item_sd:

  A numeric vector of length p of the pooled standard deviation (SD) of
  the items across groups.

- num_obs:

  A vector of length \\G\\ of sample sizes. If not `NULL`, the weights
  will be proportional to sample sizes, assuming the same weights across
  items.

- weights:

  A \\G \times p\\ matrix of weights. Default assumes equal weights
  across groups.

- group_factor:

  A vector of length \\G\\ indicating grouping for contrast. For
  example, `c(1, 1, 2)` means contrasting Group 1 & 2 vs. Group 3. The
  default is to not combine any groups, meaning the omnibus effect is
  computed.

- contrast:

  A \\p \times k\\ contrast matrix where `colSums(contrast)` = 0.
  Default is `contr.sum(p)` if `group_factor` is not specified.

- latent_mean:

  latent factor mean for the reference group. Default to 0.

- latent_sd:

  latent factor SD for the reference group. Default to 1.

- item_weights:

  Default is `NULL`. Otherwise, one can specify a vector of length \\p\\
  of weights; if so, test-level fMACS will be computed.

- thresholds:

  A matrix with G rows for measurement thresholds. The matrix must have
  column names indicating to which item index each column corresponds.

- thetas:

  Not currently used.

- link:

  Link function for the model (probit or logit).

## Value

A 1 x p matrix of fMACS effect size.

## Details

The \\f\_\text{MACS}\\ effect size is defined as \$\$f\_{\text{MACS}, i}
= \frac{1}{\mathit{SD}\_{iP}} \sqrt{\int \[(\nu\_{ij} - \bar{\nu}\_j) +
(\lambda\_{ij} - \bar{\lambda}\_j) \eta\]^2 f(\eta) d \eta}\$\$ where
\\\lambda\\ is the loading and \\\nu\\ is the intercept, and j indexes
group. The effect size reflects the square root of the ratio between the
variance in observed item score due to measurement noninvariance and the
variance of the observed item scores. \\f\_\text{MACS}\\ is analogous to
the Cohen's f effect size. When there are two groups with equal sample
sizes, \\f\_\text{MACS}\\ = \\f\_\text{MACS}\\ / 2

## Examples

``` r
lambda <- rbind(c(.7, .8, .7, .9),
                c(.7, .8, .7, .8),
                c(.8, .7, .7, .5))
nu <- rbind(c(0, .5, 0, 1),
            c(0, .2, 0, 1.1),
            c(0, .3, 0, 1.2))
fmacs(lambda,
      loadings = nu,
      pooled_item_sd = c(1, 1, 1, 1),
      latent_mean = 0,
      latent_sd = 1)
#>             [,1]      [,2] [,3]      [,4]
#> fmacs 0.04714045 0.1333333    0 0.1885618
# With contrast (Group 1 & 2 vs. Group 3)
fmacs(lambda,
      loadings = nu,
      pooled_item_sd = c(1, 1, 1, 1),
      group_factor = c(1, 1, 2),
      latent_mean = 0,
      latent_sd = 1)
#>             [,1]       [,2] [,3]      [,4]
#> fmacs 0.04714045 0.05270463    0 0.1795055
# Thresholds
lambda <- rbind(c(.8, .5, .7, .5),
                c(.8, .5, .4, .6),
                c(.8, .7, .7, .5))
tau <- rbind(c(-0.5, 0, 1, -0.3, 0.1, 0.5, -0.5, 1.5),
             c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1, 1.5),
             c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1, 0.5))
# three thresholds for items 1 and 2; one threshold for items 3 and 4
colnames(tau) <- c(1, 1, 1, 2, 2, 2, 3, 4)
fmacs_ordered(tau,
              loadings = lambda,
              pooled_item_sd = c(1, 1, 1, 1),
              latent_mean = 0,
              latent_sd = 1)
#>       [,1]       [,2]       [,3]      [,4]
#> fmacs    0 0.06906871 0.08713051 0.1169068
# With contrast (Group 1 & 2 vs. Group 3)
fmacs_ordered(tau,
              loadings = lambda,
              pooled_item_sd = c(1, 1, 1, 1),
              group_factor = c(1, 2, 1),
              latent_mean = 0,
              latent_sd = 1)
#>       [,1]       [,2]       [,3]       [,4]
#> fmacs    0 0.04029493 0.06394493 0.05383212
```
