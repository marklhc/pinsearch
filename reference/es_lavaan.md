# Item-level effect size for non-invariance

For two groups, the function uses
[`dmacs()`](https://marklhc.github.io/pinsearch/reference/dmacs.md) to
compute \\d\_\text{MACS}\\. For more than two groups, the function uses
[`fmacs()`](https://marklhc.github.io/pinsearch/reference/fmacs.md) to
compute \\f\_\text{MACS}\\, a generalisation of \\d\_\text{MACS}\\
similar to the Cohen's \\f\\ effect size.

## Usage

``` r
es_lavaan(object, ...)

pin_effsize(object, ...)
```

## Arguments

- object:

  A CFA model of class
  [`lavaan::lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan-class.html)
  fitted by [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html)

- ...:

  Additional arguments passed to
  [`dmacs()`](https://marklhc.github.io/pinsearch/reference/dmacs.md) or
  [`fmacs()`](https://marklhc.github.io/pinsearch/reference/fmacs.md)

## Value

A matrix of 1 row showing the effect size values for each non-invariant
item on each latent variable.
