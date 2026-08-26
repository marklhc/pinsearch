# Map invariance types to lavaan parTable operators.

Maps the invariance type strings used by
[`pinSearch()`](https://marklhc.github.io/pinsearch/reference/pinSearch.md)
to the corresponding values of the `op` column of a `lavaan` parameter
table
([`lavaan::parTable()`](https://rdrr.io/pkg/lavaan/man/parTable.html)):
`"loadings"` maps to `"=~"`, `"intercepts"` to `"~1"`, `"thresholds"` to
`"|"`, and both `"residuals"` and `"residual.covariances"` to `"~~"`.

## Usage

``` r
type2op(type)
```

## Arguments

- type:

  A character string indicating the invariance type. Must be one of
  `"loadings"`, `"intercepts"`, `"thresholds"`, `"residuals"`, or
  `"residual.covariances"`.

## Value

A character string, the corresponding lavaan operator.

## Examples

``` r
type2op("loadings")
#> [1] "=~"
```
