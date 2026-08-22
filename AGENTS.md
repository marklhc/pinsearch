# AGENTS.md

R package `pinsearch`: specification search for partial factorial invariance, built on top of `lavaan` model syntax.

## Commands

- Run all tests: `Rscript -e 'testthat::test_local()'` (or `devtools::test()`)
- Run one test file: `Rscript -e 'testthat::test_local(filter = "dmacs")` (filter = test file stem)
- Regenerate docs/NAMESPACE: `devtools::document()`
- Full check (what CI runs): `Rscript -e 'rcmdcheck::rcmdcheck()'`
- Rebuild README: edit `README.Rmd` then `devtools::build_readme()` — `README.md` is generated, never edit it directly
- No linter or formatter is configured; CI is only R-CMD-check + pkgdown

## Architecture

- `R/pinSearch.R` — main entrypoint `pinSearch()`. Model strings are passed straight to `lavaan::cfa()` (stage fits go through the `fit_stage()` helper); it sequentially frees invariance constraints (loadings → intercepts/thresholds → residuals/residual covariances) per Yoon & Millsap (2007). `std.lv = TRUE` is injected by default for continuous models. Model strings may be given as a character vector and are collapsed internally; `group:`-prefixed block syntax selects the group-as-block identification path.
- `R/dmacs.R` / `R/fmacs.R` — effect sizes: `dmacs()`/`dmacs_ordered()` (two groups), `fmacs()`/`fmacs_ordered()` (>2 groups). `pin_effsize()` (R/dmacs.R:465) aliases `es_lavaan()`. Both `effect_size` and `progress` args are documented as experimental.
- `R/remove_cons.R` — frees constraints inside a lavaan `parTable`; `R/helper.R` reads lavaan S4 slots.
- `data-raw/DATASET.R` regenerates the `lui_sim` dataset in `data/` (needs `haven`, `dplyr`, downloads from OSF); `data-raw` is never built.

## Testing quirks

- Test fixtures are pre-saved in `tests/testthat/testdata/*.rds`, loaded by `helper-testdata.R`. The regeneration code is gated behind `if (FALSE)` in that file — flip to `TRUE` and run it (needs `MASS`, `withr`) to regenerate, then flip back.
- Tests fit real lavaan models (simulated data up to n = 10000/group), so the full suite is slow; prefer `test_local(filter = ...)` while iterating.
- `difR` is in Suggests and used only by the vignette `vignettes/compare-cont-cat.Rmd`.

## Repo conventions

- Default branch is `master` (CI also triggers on `main`); PRs come from `marklhc/issueNN` branches.
- After a fix, expect two follow-up commits: bump `Version` in `DESCRIPTION` (dev versions use a 4th field, e.g. `0.1.4.3`) and re-render the README.
- Vignettes live in `vignettes/*.Rmd` with pre-saved outputs (`.rds`); `vignettes/*.R` and `*.html` are gitignored.
- Dependencies: `lavaan` and `stats` are hard Imports — do not add to Imports casually; new optional deps go in Suggests.

## Test environment note

- `test-categorical.R` had 4 failures/errors on lavaan 0.7-2 (a mis-specified ordinal threshold stage produced a degenerate, `Pr(>Chisq) = NA` stage LRT). Fixed on `master` (PR #21) by tying `group.equal = "intercepts"` from the threshold stage onward. It now passes; if it fails, treat it as a real regression rather than an environment issue.

## Related packages

- Companion package `longcfa` (github.com/marklhc/longcfa, longitudinal CFA helpers) declares `Imports: pinsearch (>= 0.1.6)` and reuses the exported `fdr_alpha()` and `type2op()`. Keep those two functions' signatures/semantics stable (they are part of the public API now) or coordinate the change with longcfa. CRAN order: pinsearch must be archived before longcfa.

## Open work (as of 2026-08)

- None. #15 "Longitudinal CFA" was closed in favor of the `longcfa` package (the `marklhc/issue15` draft branch is obsolete; its semTools equivalence test was ported to longcfa).

Recently merged: #22 (refactored model-fitting calls in `pinSearch()` into `fit_stage()`), #16 "group syntax" (`min2`, LRT fixes, group-specific model syntax), #7/#4 multi-latent effect sizes (PR #20), #21 (lavaan 0.7 compat).
