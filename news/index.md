# Changelog

## pinsearch (development version)

## pinsearch 0.1.7

CRAN release: 2026-11-19

- The test-level effect size in
  [`dmacs()`](https://marklhc.github.io/pinsearch/reference/dmacs.md)
  and
  [`dmacs_ordered()`](https://marklhc.github.io/pinsearch/reference/dmacs.md)
  now handles missing group sample sizes.

- Fixed the denominator of the test-level effect size in
  [`dmacs()`](https://marklhc.github.io/pinsearch/reference/dmacs.md),
  [`dmacs_ordered()`](https://marklhc.github.io/pinsearch/reference/dmacs.md),
  [`fmacs()`](https://marklhc.github.io/pinsearch/reference/fmacs.md),
  and
  [`fmacs_ordered()`](https://marklhc.github.io/pinsearch/reference/fmacs.md)
  (via `item_weights`): it now uses the pooled standard deviation of the
  weighted item sum, including the cross-item covariances induced by the
  common factor and the correct (quadratic) weight scaling.
  [`es_lavaan()`](https://marklhc.github.io/pinsearch/reference/es_lavaan.md)/[`pin_effsize()`](https://marklhc.github.io/pinsearch/reference/es_lavaan.md)
  compute it exactly from the sample covariances for continuous
  indicators; for ordered indicators the observed-scale variance-only
  fallback is used, as `sampstat$cov` is on the latent-response scale
  there. A new `pooled_test_sd` argument allows passing this value
  directly (it takes precedence over the sample-based value).

- The test-level effect size in
  [`dmacs()`](https://marklhc.github.io/pinsearch/reference/dmacs.md)
  and
  [`fmacs()`](https://marklhc.github.io/pinsearch/reference/fmacs.md)
  now works when `loadings` is not supplied (as it already did at the
  item level).

- Added a `NEWS.md` file to track changes to the package.
