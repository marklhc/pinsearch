# Compute dMACS effect size described in Nye & Drasgow (2011) for two groups.

`dmacs` returns the dMACS effect size statistics given a set of loadings
and intercepts.

## Usage

``` r
dmacs(
  intercepts,
  loadings = NULL,
  pooled_item_sd = NULL,
  latent_mean = 0,
  latent_sd = 1,
  uniqueness = NULL,
  ns = NULL,
  item_weights = NULL
)

dmacs_ordered(
  thresholds,
  loadings,
  thetas = 1,
  link = c("probit", "logit"),
  pooled_item_sd = NULL,
  latent_mean = 0,
  latent_sd = 1,
  item_weights = NULL
)
```

## Arguments

- intercepts:

  A \\2 \times p\\ matrix of measurement intercepts.

- loadings:

  A \\2 \times p\\ matrix of factor loadings, where p is the number of
  items.

- pooled_item_sd:

  A numeric vector of length p of the pooled standard deviation (SD) of
  the items across groups.

- latent_mean:

  latent factor mean for the reference group. Default to 0.

- latent_sd:

  latent factor SD for the reference group. Default to 1.

- uniqueness:

  A vector of length \\p\\ of uniqueness.

- ns:

  A vector of length \\p\\ of sample sizes.

- item_weights:

  Default is `NULL`. Otherwise, one can specify a vector of length \\p\\
  of weights; if so, test-level dMACS will be computed.

- thresholds:

  A matrix with two rows for measurement thresholds. The matrix must
  have column names indicating to which item index each column
  corresponds.

- thetas:

  Not currently used.

- link:

  Link function for the model (probit or logit).

## Value

A 1 x p matrix of dMACS effect size. If `item_weights` is not `NULL`,
\\p\\ = 1.

## Details

The \\d\_\text{MACS}\\ effect size is defined as (Nye & Drasgow, 2011,
p. 968) \$\$d\_{\text{MACS}, i} = \frac{1}{\mathit{SD}\_{iP}} \sqrt{\int
\[(\nu\_{iR} - \nu{iF}) + (\lambda\_{iR} - \lambda\_{iF}) \eta\]^2
f(\eta) d \eta}\$\$ where \\\lambda\\ is the loading and \\\nu\\ is the
intercept, F and R denote the focal and the reference group. The effect
size reflects the standardized mean difference on an item due to
measurement noninvariance, and is analogous to the Cohen's d effect
size.

## References

Nye, C. & Drasgow, F. (2011). Effect size indices for analyses of
measurement equivalence: Understanding the practical importance of
differences between groups. Journal of Applied Psychology, 96(5),
966-980.

## Examples

``` r
lambdaf <- c(.8, .5, .7, .5)
lambdar <- c(.8, .5, .4, .6)
nuf <- c(0.1, 0, 0.2, 0)
nur <- c(0.2, 0, 0, 0)
dmacs(rbind(nuf, nur),
      loadings = rbind(lambdaf, lambdar),
      pooled_item_sd = c(1, 1, 1, 1),
      latent_mean = 0,
      latent_sd = 1)
#>       [,1] [,2]      [,3] [,4]
#> dmacs  0.1    0 0.3605551  0.1
dmacs(rbind(nuf, nur),
      loadings = rbind(lambdaf, lambdar),
      pooled_item_sd = c(1, 1, 1, 1),
      latent_mean = 0,
      latent_sd = 1,
      item_weights = c(1, 1, 1, 1))
#>        item_sum
#> dmacs 0.1118034
# Thresholds
lambda <- rbind(c(.8, .5, .7, .5),
                c(.8, .5, .4, .6))
tau <- rbind(c(-0.5, 0, 1, -0.3, 0.1, 0.5, -0.5, 1.5),
             c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1, 1.5))
# three thresholds for items 1 and 2; one threshold for items 3 and 4
colnames(tau) <- c(1, 1, 1, 2, 2, 2, 3, 4)
dmacs_ordered(tau,
              loadings = lambda,
              pooled_item_sd = c(1, 1, 1, 1),
              latent_mean = 0,
              latent_sd = 1)
#>       [,1]       [,2]      [,3]       [,4]
#> dmacs    0 0.01711254 0.2019787 0.02276356
```
