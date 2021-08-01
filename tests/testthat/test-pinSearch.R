library(lavaan)

# Initialize an example
lambda1 <- seq(.9, .7, length.out = 5)
lambda2 <- c(lambda1[1], 1, lambda1[3:5])
cov1 <- tcrossprod(lambda1) + diag(.5, 5)
dimnames(cov1) <- list(paste0("y", 1:5), paste0("y", 1:5))
nu1 <- seq(-0.5, 0.5, length.out = 5)
nu2 <- c(nu1[1], 0.25, nu1[3:4], 0.3)
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
  expect_identical(ps2[[2]]$rhs, c("y2"))
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
  expect_identical(ps3[[2]]$lhs, c("f", "y5"))
  expect_identical(ps3[[2]]$rhs, c("y2", ""))
  expect_identical(ps3[[2]]$type, c("loadings", "intercepts"))
})

test_that("pinSearch() works properly for three groups", {
  cov2 <- tcrossprod(lambda2) * 1.3 + diag(.5, 5)
  mean2 <- lambda2 * .4 + nu2
  cov3 <- tcrossprod(lambda2) * 0.8 + diag(.5, 5)
  dimnames(cov2) <- dimnames(cov3) <- list(paste0("y", 1:5), paste0("y", 1:5))
  nu3 <- c(-0.8, nu1[2:5])
  mean3 <- lambda2 * -.3 + nu3
  ps4 <- pinSearch(' f =~ y1 + y2 + y3 + y4 + y5',
                   sample.cov = list(cov1, cov2, cov3),
                   sample.mean = list(mean1, mean2, mean3),
                   sample.nobs = c(10000, 10000, 10000),
                   type = "intercepts")
  expect_identical(ps4[[2]]$lhs, c("f", "y5", "y2", "y1"))
  expect_identical(ps4[[2]]$rhs, c("y2", "", "", ""))
  expect_identical(ps4[[2]]$type, rep(c("loadings", "intercepts"), c(1, 3)))
})

test_that("pinSearch() works properly for noninvariant residual covariances", {
  cov4 <- tcrossprod(lambda1) + diag(.5, 5)
  lambda2 <- c(.4, 1, lambda1[-(1:2)])
  cov4[1, 2] <- cov4[2, 1] <- cov4[1, 2] + 0.1
  cov4[1, 5] <- cov4[5, 1] <- cov4[1, 5] - 0.1
  lambda5 <- c(.4, lambda1[2:4], 1)
  cov5 <- tcrossprod(lambda5) * 1.3 + diag(c(.5, .75, .5, .6, .5))
  cov5[1, 2] <- cov5[2, 1] <- cov5[1, 2] + 0.2
  cov5[1, 5] <- cov5[5, 1] <- cov5[1, 5] - 0.2
  dimnames(cov4) <- dimnames(cov5) <- list(paste0("y", 1:5), paste0("y", 1:5))
  mean4 <- mean1
  mean5 <- lambda2 * .4 + nu2
  ps5 <- pinSearch(' f =~ y1 + y2 + y3 + y4 + y5
                     y1 ~~ y2
                     y1 ~~ y5 ',
                   sample.cov = list(cov4, cov5),
                   sample.mean = list(mean4, mean5),
                   sample.nobs = c(10000, 10000),
                   type = "residual.covariances")
  expect_identical(ps5[[2]]$lhs, c("f", "f", "y2", "y2", "y4", "y1"))
  expect_identical(ps5[[2]]$rhs, c("y1", "y5", "", "y2", "y4", "y5"))
  expect_identical(ps5[[2]]$type,
                   rep(c("loadings", "intercepts", "residuals",
                         "residual.covariances"), c(2, 1, 2, 1)))
})

