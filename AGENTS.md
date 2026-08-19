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

- `R/pinSearch.R` — main entrypoint `pinSearch()`. Model strings are passed straight to `lavaan::cfa()`; it sequentially frees invariance constraints (loadings → intercepts/thresholds → residuals/residual covariances) per Yoon & Millsap (2007). `std.lv = TRUE` is injected by default for continuous models.
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

## Open work (as of 2026-08)

- #15 "Longitudinal CFA" (`marklhc/issue15` branch): `longcfa()`/`longcfa_syntax()` draft with tests + vignette; still needs categorical items. Deferred to a later wave.
- #7/#4 multi-latent effect sizes: fixed in local branch `marklhc/issue7` (per-(item, latent) columns, corrected latent scaling, regression tests) — not yet pushed; PR pending.
- #16 (draft PR) "group syntax": local branch `group-specific-syntax` is rebased onto `master` (worktree `/tmp/opencode/gss`), `min2` + LRT fixes + group syntax verified end-to-end (also covers #14's differing observed variables scenario) — not yet pushed; un-draft the PR after pushing. Note: its `remove_cons()` needs `op2col()`, which lives in `R/remove_cons.R`.
- `cat` branch is fully merged (stale, safe to delete once push access is available).
