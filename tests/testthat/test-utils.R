test_that("type2op() maps invariance types to lavaan operators", {
    expect_identical(type2op("loadings"), "=~")
    expect_identical(type2op("intercepts"), "~1")
    expect_identical(type2op("thresholds"), "|")
    expect_identical(type2op("residuals"), "~~")
    expect_identical(type2op("residual.covariances"), "~~")
})

test_that("type2op() errors on invalid input", {
    expect_error(type2op("loadings!"))
    expect_error(type2op("Thresholds"))
})

test_that("fdr_alpha() computes the right number", {
    # Hand-computed: 3 * .05 / (10 + 1 - 3 * .95) = .15 / 8.15
    expect_equal(fdr_alpha(3, m = 10, q = .05), 3 * .05 / (10 + 1 - 3 * .95))
    # Degenerate case (single test): .05 / (1 + 1 - 1 * (1 - .05)) = .05/1.05
    # = 1/21 (~0.0476); it comes out to 1/21, not q itself
    expect_equal(fdr_alpha(1, m = 1, q = .05), .05 / 1.05)
    # Boundary (last step): 10 * .05 / (10 + 1 - 10 * (1 - .05))
    # = .5 / (11 - 9.5) = .5 / 1.5 = 1/3
    expect_equal(fdr_alpha(10, m = 10, q = .05), 1/3)
})

test_that("fdr_alpha() is strictly increasing in i", {
    alphas <- fdr_alpha(1:20, m = 20, q = .05)
    expect_true(!any(diff(alphas) <= 0))
})

test_that("fdr_alpha() is vectorized in i", {
    alphas <- fdr_alpha(1:5, m = 20)
    expect_length(alphas, 5)
    expect_equal(alphas, vapply(1:5, fdr_alpha, numeric(1), m = 20))
})

test_that("type2op and fdr_alpha are exported", {
    expect_true("type2op" %in% getNamespaceExports("pinsearch"))
    expect_true("fdr_alpha" %in% getNamespaceExports("pinsearch"))
})
