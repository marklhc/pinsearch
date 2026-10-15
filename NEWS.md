# pinsearch 0.1.6.2

* Fixed the denominator of the test-level effect size in `dmacs()`,
  `dmacs_ordered()`, `fmacs()`, and `fmacs_ordered()` (via `item_weights`):
  it now uses the pooled standard deviation of the weighted item sum,
  including the cross-item covariances induced by the common factor and the
   correct (quadratic) weight scaling. `es_lavaan()`/`pin_effsize()` compute
   it exactly from the sample covariances for continuous indicators; for
   ordered indicators the observed-scale variance-only fallback is used,
   as `sampstat$cov` is on the latent-response scale there. A new
   `pooled_test_sd` argument allows passing this value directly (it takes
   precedence over the sample-based value).

* The test-level effect size in `dmacs()` and `fmacs()` now works when
  `loadings` is not supplied (as it already did at the item level).

* Added a `NEWS.md` file to track changes to the package.
