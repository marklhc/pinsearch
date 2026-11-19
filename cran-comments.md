## R CMD check results

0 errors | 0 warnings | 0 notes

This is a new release (0.1.7) of the existing package pinsearch.

Changes since 0.1.6:
- Fixed the denominator of the test-level effect size (item_weights) and
  added a `pooled_test_sd` argument to dmacs()/fmacs() and ordered variants.
- Test-level effect size now works when `loadings` is not supplied.
- Test-level dMACS now handles missing group sample sizes.
