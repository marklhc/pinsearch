library(lavaan)

# Initialize an example
lambda1 <- seq(.9, .7, length.out = 5)
lambda2 <- c(.4, 1, lambda1[-(1:2)])
cov1 <- tcrossprod(lambda1) + diag(.5, 5)
dimnames(cov1) <- list(paste0("y", 1:5), paste0("y", 1:5))
nu1 <- seq(-0.5, 0.5, length.out = 5)
nu2 <- c(nu1[1:3], 0, 0)
mean1 <- nu1

test_that("pinSearch() works properly for invariant data", {
  cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 5)
  dimnames(cov2) <- list(paste0("y", 1:5), paste0("y", 1:5))
  mean2 <- lambda1 * .4 + nu1
  ps1 <- pinSearch(' f =~ y1 + y2 + y3 + y4 + y5',
                   sample.cov = list(cov1, cov2),
                   sample.mean = list(mean1, mean2),
                   sample.nobs = c(100, 100))
  expect_equal(nrow(ps1[[2]]), 0)
})

test_that("pinSearch() works properly for noninvariant loadings", {
  cov2 <- tcrossprod(lambda2) * 1.3 + diag(.5, 5)
  dimnames(cov2) <- list(paste0("y", 1:5), paste0("y", 1:5))
  mean2 <- lambda2 * .4 + nu1
  ps2 <- pinSearch(' f =~ y1 + y2 + y3 + y4 + y5',
                   sample.cov = list(cov1, cov2),
                   sample.mean = list(mean1, mean2),
                   sample.nobs = c(10000, 10000))
  expect_identical(ps2[[2]]$item, c("y1", "y2"))
})

test_that("pinSearch() works properly for noninvariant intercepts", {
  cov2 <- tcrossprod(lambda2) * 1.3 + diag(.5, 5)
  dimnames(cov2) <- list(paste0("y", 1:5), paste0("y", 1:5))
  mean2 <- lambda2 * .4 + nu2
  ps3 <- pinSearch(' f =~ y1 + y2 + y3 + y4 + y5',
                   sample.cov = list(cov1, cov2),
                   sample.mean = list(mean1, mean2),
                   sample.nobs = c(10000, 10000),
                   type = "intercepts")
  expect_identical(ps3[[2]]$item, c("y1", "y2", "y5", "y4"))
  expect_identical(ps3[[2]]$type, rep(c("loadings", "intercepts"), each = 2))
})

