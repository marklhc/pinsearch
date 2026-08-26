# Adjusted alpha level for multiple testing.

Computes the adjusted significance level at step \\i\\ of \\m\\ tests in
the sequential selection procedure of Benjamini & Gavrilov (2009), as
used by
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md)
when `control_fdr = TRUE`:

## Usage

``` r
fdr_alpha(i, m, q = 0.05)
```

## Arguments

- i:

  Step index (rank) of the current test in the sequential search.

- m:

  Total number of tests (maximum number of constraints that can be
  freed).

- q:

  Target false discovery rate. Default is .05.

## Value

A numeric scalar, the adjusted alpha level at step \\i\\.

## Details

\$\$\frac{iq}{m + 1 - i(1 - q)}\$\$ where \\i\\ is the step index in the
search, \\m\\ is the maximum number of constraints that can be freed,
and \\q\\ is the target significance level.

Note that \\i\\ is defined for \\1 \<= i \<= m\\, with \\m \>= 1\\ and
\\0 \< q \< 1\\. These bounds are maintained by the internal logic of
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md),
but external callers must ensure them on their own.

## References

Benjamini, Y. & Gavrilov, N. M. (2009). Sequential selection procedures
for testing dependent hypotheses.

## Examples

``` r
fdr_alpha(1, 10)
#> [1] 0.004975124
```
