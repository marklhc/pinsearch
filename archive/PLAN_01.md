# Plan 01 — lavaan 0.7-2 compatibility fix (master-level, merge first)

> This document is self-contained. A new session should be able to execute it
> end-to-end without any prior context. Read the **Context**, **Current state**,
> and **Decisions** sections before starting.

## Objective

Make `pinsearch` work under the current CRAN **lavaan 0.7-2** (published
2026-07-16), which breaks the ordinal (threshold) invariance-search stage gate
in `pinSearch()`. Ship this as a **master-level PR, merged before the two open
PRs are rebased**, so that both open PRs' `R-CMD-check` go fully green.

This is the single remaining thing keeping R-CMD-check red on both open PRs
(#16 `group-specific-syntax`, #20 `marklhc/issue7`).

## Context — why this is needed (root cause, already confirmed)

All 4 pre-existing `tests/testthat/test-categorical.R` failures occur under
lavaan 0.7-2 / R 4.6.1 and are **not** caused by any branch-specific change —
they reproduce on unmodified `master`.

Root cause (confirmed by direct probing of the installed lavaan 0.7-2):
the invariance-stage gate in `pinSearch()` uses `lavaan::lavTestLRT(base, new)`.
For ordinal models where a stage adds `group.equal = "thresholds"`, lavaan
0.7-2 reports the tied-threshold constraints as **not reducing the model's df**
(tied parameters keep a `free` index and are counted), so the LRT degenerates:
`Df diff = 0` and `Chisq diff ≈ 0` (or `Pr(>Chisq) = NA` when `Chisq diff` is
nonzero). That produces exactly the 4 observed failures:

1. `test-categorical.R:22` (binary, thresholds) — `Chisq diff≈0 && Df diff==0`
   triggers the "All free X are noninvariant!" branch → threshold stage is
   **skipped** → detection returns **empty** (`ps1[[2]]$lhs` is `character(0)`).
2. `test-categorical.R:40` (noninvariant uniqueness) — same stage-skip → `yy7`
   not detected.
3. `test-categorical.R:44` (noninvariant unique covariances) — `Chisq diff`
   nonzero but `Df diff == 0` → `Pr(>Chisq)` is `NA` → `if (NA >= sig_level)`
   hard-errors: *"missing value where TRUE/FALSE needed"* (this is the one that
   aborts the whole test run).
4. `test-categorical.R:97` — ~7% shift in the hard-coded fMACS 4-vector
   (estimator drift; detection itself is correct in this 3-group case).

Key facts that shape the fix:
- The per-parameter `mod` test (default `inv_test`) **still detects the planted
  non-invariances correctly** under 0.7-2 (proven: the 3-group case finds the
  right items). So falling back to the per-parameter test is sound.
- The package only reads `fit@test$standard$df` (still a valid scalar in 0.7-2)
  and never uses `parTable` `std.*` columns — so **no other code exposure** from
  the 0.7-2 `@test` restructure / `parTable` column removal.
- Continuous (non-ordinal) models are **unaffected** (their LRT stays
  well-defined, `p` is not `NA`) — so the rest of the suite is a regression
  safety net that must stay bit-identical.

## Decisions (already made — do not re-litigate)

- **Version strategy:** *adapt the code* to lavaan ≥ 0.7. Do **not** pin lavaan
  and do **not** add a version cap in `DESCRIPTION` (the fix is backward
  compatible via `is.na()` guards, so lavaan < 0.7 must stay byte-identical).
- **Branching:** a **separate master-level PR**, merged **first**, then rebase
  #16 and #20.
- **Orchestration:** flat (orchestrator launches each subagent itself and
  commits between phases), NOT architect-led delegation — this keeps clean
  one-logical-change-per-commit boundaries (the architect and doc subagents both
  touch `R/pinSearch.R`, so architect-led would merge their diffs).
- **Commits on the branch:** 4 separate commits (code / docs / tests / AGENTS).
- Subagents **never** commit or push. Only the orchestrator touches git.

## Current state (starting point)

- Repo: `/home/marklai/pinsearch`. `master` = `033de11`, `DESCRIPTION`
  `Version: 0.1.4.3`.
- Branches/PRs already in flight (do **not** touch these until Phase 4):
  - `marklhc/issue7` (PR #20) = local `11ac20d` (0.1.4.7), 2 commits ahead of
    `origin`.
  - `group-specific-syntax` (PR #16) = worktree `/tmp/opencode/gss` at local
    `754f06c` (0.1.5.1), 1 commit ahead of `origin`.
- The 4 categorical failures are present on `master` under current CRAN.
- Subagents (defined in `~/.config/opencode/agents/`, **do not edit them**):
  - `r-architect` — edits `R/` + `DESCRIPTION` (NOT `man/`/`NAMESPACE`); owns
    strategy, cross-file consistency, dependency mgmt; runs `devtools::check()`
    as final gate.
  - `r-tester` — edits **only** `tests/testthat/`; runs `devtools::load_all()`
    then `devtools::test()`; must not "fix" a test to match broken behavior
    (report bugs to architect).
  - `r-doc` — edits **only** `R/*` roxygen `#'` blocks, regenerates `man/`+
    `NAMESPACE` via `devtools::document()`; runs `document()`+`test()`, not
    `check()`.

## Execution

### Phase 0 — Branch (orchestrator)
- From `master` (`033de11`), create branch `marklhc/lavaan-0.7-compat`.
- Verify clean worktree first (`git status`). No commits yet.

### Phase 1 — Code (subagent: `r-architect`) — then orchestrator commit 1
Instruct `r-architect` (code slice only; tell it explicitly: *"I am launching
@R-Tester and @R-Doc separately — do NOT delegate or edit tests/docs; implement
only the R/ and DESCRIPTION changes; verify with `devtools::load_all()` before
returning"*):

- **Stage gate in `R/pinSearch.R` (~lines 265-274).** Divert **only** when the
  LRT is degenerate (`Pr(>Chisq)` is `NA`) so lavaan < 0.7 is byte-identical:
  - Add `&& !is.na(lrt_p)` to the
    `if (isTRUE(all.equal(lrt_base_new[2, "Chisq diff"], 0)) && df_diff == 0)`
    "All free X are noninvariant!" branch (that branch is the *degenerate* case
    under 0.7-2, not genuine invariance).
  - Gate the `if (lrt_base_new[2, "Pr(>Chisq)"] >= sig_level)` accept-stage
    branch on `!is.na(lrt_p)` as well.
  - The `else` falling through to the existing per-parameter `fn_get_inv`
    freeing loop now handles both `p < sig` (old behavior) **and** `p = NA`
    (degenerate); that loop already accepts-when-nothing-significant and sets
    `base_fit <- new_fit` — do not restructure the loop.
  - Guard the FDR path: when `control_fdr` is `TRUE` and `df_diff < 1`, set
    `p_enter <- sig_level` instead of `fdr_alpha(num_free, m = df_diff, q =
    sig_level)` (FDR undefined at 0 df; `fdr_alpha(i, m=0)` degenerates to 1,
    over-permissive).
- **NA guard in `get_invlrt()` in `R/remove_cons.R` (~lines 110-142).** Single-
  constraint ordinal LRTs also return `Df diff = 0` under 0.7-2, so `Pr(>Chisq)`
  can be `NA`. Before `which.min`, filter `NA` `Pr(>Chisq)` rows out of
  `lrt_mat`/`mis`; if none remain, `return(NULL)` (stage accepted, no error).
  Apply to the active path (and the `if (FALSE)` block for consistency).
- **Cross-file consistency:** confirm the gate's fall-through calls the guarded
  `get_inv*` functions correctly.
- **Dependency mgmt:** bump `DESCRIPTION` `Version` `0.1.4.3 → 0.1.4.4`. Do not
  add any lavaan version constraint.
- Verify `devtools::load_all()` builds before returning.

Orchestrator: inspect diff (must be `R/pinSearch.R`, `R/remove_cons.R`,
`DESCRIPTION` only; no `man/`/`NAMESPACE` churn), then
**commit 1**: `Fix invariance stage-gate for lavaan >= 0.7`.

### Phase 2 — Docs (subagent: `r-doc`) — then orchestrator commit 2
Instruct `r-doc` (scope: roxygen on `pinSearch` + `devtools::document()`; note
the 4 categorical tests are still red at this point — expected, values are
updated by @R-Tester next):
- Add a user-facing note to `pinSearch`'s roxygen (e.g. an `@details` paragraph
  or an existing `@note`-style spot): under **lavaan ≥ 0.7**, group-constrained
  ordinal (threshold) models can yield a degenerate LRT (`Df diff = 0`); the
  search then falls back to the per-parameter test. Note the `inv_test = "lrt"`
  + ordinal caveat (degrades gracefully when the single-constraint LRT is
  degenerate).
- `devtools::document()`; confirm **only** the intended `man/pinSearch.Rd`
  (and, if the note adds/removes an `@importFrom`, `NAMESPACE`) changed.
- Run `devtools::test()` to confirm the roxygen change broke nothing (the 4
  categorical failures are expected-pending here, not caused by docs).

Orchestrator: inspect diff (roxygen block in `R/pinSearch.R` + regenerated
`man/`/`NAMESPACE` only), then
**commit 2**: `Document lavaan >= 0.7 compatibility`.

### Phase 3 — Tests (subagent: `r-tester`) — then orchestrator commit 3
Instruct `r-tester` (scope: `tests/testthat/test-categorical.R` only; run
`devtools::load_all()` then `devtools::test()` until green — this run now covers
all `R/` changes):
- **Regenerate the hard-coded fMACS vector at `test-categorical.R:97-99`:**
  run the fixed `pinSearch(...)` on the `dfo` 3-group fixture, capture
  `c(pin_effsize(ps5[[1]]))` (the 4 values), and replace the literal. Keep the
  version-independent cross-check vs `ps5_re` (line 100).
- **Confirm item-set expectations** at lines 22, 40, 44 now pass under restored
  detection. If detection legitimately differs from the old literal, adjust the
  expectation **and** report why (do not paper over a real regression).
- **Add a regression test** for the degenerate path: assert no hard `if(NA)`
  error on an `NA Pr(>Chisq)` stage and that the per-parameter fallback fires
  (guards against re-introducing the error). Match existing fixture/convention
  style (reuse `df1`/`dfo` fixtures; no new `Suggests`).
- **Regression safety net:** confirm the **continuous** tests (the rest of the
  suite) are green with **no new numeric drift** (the fix must not move
  continuous values).

Orchestrator: inspect diff (tests only), then
**commit 3**: `Update categorical test expectations for lavaan 0.7`.

### Phase 4 — Final gate, AGENTS, ship (orchestrator + `r-architect`)
1. **Re-invoke `r-architect`** to run `devtools::check()` over the whole
   committed change (its assigned final gate). All 4 categorical tests must be
   green; treat any **new** WARNING/NOTE relative to the pre-existing baseline
   as blocking.
2. Orchestrator updates `AGENTS.md` "Test environment note": record the
   lavaan 0.7-2 incompatibility and that it is now **resolved in-code** (so a
   future agent doesn't re-chase the 4 categorical failures).
   **commit 4**: `Note lavaan >= 0.7 resolution in AGENTS.md`.
3. `git push -u origin marklhc/lavaan-0.7-compat`; open PR to `master`
   (describe the lavaan-0.7 degeneracy in the body; optionally file a GitHub
   issue and reference it).
4. **Merge to master first**, then rebase `marklhc/issue7` (0.1.4.7) and
   `group-specific-syntax` (0.1.5.1) onto the fixed master (both versions stay
   ahead; the 3 pending commits ride along in the rebases).
   - The `group-specific-syntax` PR has the most `pinSearch` overlap (its
     `group_as_block` logic sits just above this gate) — resolve conflicts
     carefully.
5. Verify both PRs' R-CMD-check now goes fully green (CI tests the
   PR-vs-master merge ref).

## Out of scope (do not do here)
- Pushing/PR-ing `marklhc/issue7` or `group-specific-syntax` (happens in
  Phase 4 step 4, as rebases of the existing PRs).
- `#15 Longitudinal CFA`, the `add-agents-md` branch, and any lavaan *pinning*.

## Verification checklist (done when all true)
- [ ] `marklhc/lavaan-0.7-compat` has exactly 4 commits (code/docs/tests/AGENTS).
- [ ] No `man/` or `NAMESPACE` churn beyond what `r-doc`'s intended note
      regenerated.
- [ ] All 4 `test-categorical.R` tests pass; rest of the suite green with no
      new numeric drift on continuous tests.
- [ ] `devtools::check()` clean (no new WARNINGs/NOTEs vs baseline).
- [ ] Compat PR merged to `master` **before** #16/#20 rebase.
- [ ] After rebase: R-CMD-check green on both #16 and #20.
