library(lavaan)

test_that("dmacs() computes the right number", {
  d1 <- dmacs(rbind(c(0.5, 1), c(0.5, 1)), pooled_item_sd = 1)
  d2 <- dmacs(rbind(c(0.3, 0.8), c(0.5, 1)), pooled_item_sd = 1)
  d3 <- dmacs(rbind(0.8, 0.5), loadings = rbind(1, 0.5),
              pooled_item_sd = 2,
              latent_mean = 0.5, latent_sd = 2)
  expect_equal(as.vector(d1), c(0, 0))
  expect_equal(as.vector(d2), c(0.2, 0.2))
  expect_equal(as.vector(d3),
               sqrt(0.3^2 + 2 * 0.3 * 0.5 * 0.5 + 0.5^2 * (0.5^2 + 2^2)) / 2)
})

# Initialize an example
lambda1 <- seq(.9, .7, length.out = 5)
lambda2 <- c(.4, 1, lambda1[-(1:2)])
cov1 <- tcrossprod(lambda1) + diag(.5, 5)
dimnames(cov1) <- list(paste0("y", 1:5), paste0("y", 1:5))
nu1 <- seq(-0.5, 0.5, length.out = 5)
nu2 <- c(nu1[1:3], 0, 0)
mean1 <- nu1

test_that("es_lavaan() works properly for invariant data", {
  cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 5)
  dimnames(cov2) <- list(paste0("y", 1:5), paste0("y", 1:5))
  mean2 <- lambda1 * .4 + nu1
  scalar1 <- cfa(' f =~ y1 + y2 + y3 + y4 + y5',
                 sample.cov = list(cov1, cov2),
                 sample.mean = list(mean1, mean2),
                 sample.nobs = c(100, 100),
                 group.equal = c("loadings", "intercepts"))
  expect_length(es_lavaan(scalar1), 0)
})

test_that("es_lavaan() works properly for noninvariant data", {
  cov2 <- tcrossprod(lambda2) * 1.3 + diag(.5, 5)
  dimnames(cov2) <- list(paste0("y", 1:5), paste0("y", 1:5))
  mean2 <- lambda2 * .4 + nu2
  ps3 <- pinSearch(' f =~ y1 + y2 + y3 + y4 + y5',
                   sample.cov = list(cov1, cov2),
                   sample.mean = list(mean1, mean2),
                   sample.nobs = c(10000, 10000),
                   type = "intercepts")
  dmacs_ps3 <- es_lavaan(ps3$`Partial Invariance Fit`)
  expect_length(dmacs_ps3, 4)
  expect_gt(dmacs_ps3[1, "y1-f"], dmacs_ps3[1, "y2-f"])
  expect_lt(dmacs_ps3[1, "y4-f"], dmacs_ps3[1, "y5-f"])
})

# Ordered items
lambda <- rbind(c(1.323, 0.875), c(1.323, 0.875))
thres <- rbind(c(-2.211, -0.728, 1.468, -0.014, 0.404, 1.438),
               c(-2.211, -0.728, 1.468, -0.635, 0.404, 1.438))
colnames(thres) <- rep(1:2, each = 3)
test_that("dmacs_ordered() computes a sensible number", {
    d4 <- dmacs_ordered(thres, loadings = lambda, pooled_item_sd = 1)
    expect_equal(c(d4), c(0, 0.19029), tolerance = 0.0001)
})

test_that("dmacs_ordered() works for binary items", {
    thres_bin <- rbind(c(-2.211, 1.438),
                       c(-2.211, 1.7))
    colnames(thres_bin) <- c(1, 2)
    d5 <- dmacs_ordered(thres_bin, loadings = lambda, pooled_item_sd = 1)
    expect_length(d5, n = 2)
})
