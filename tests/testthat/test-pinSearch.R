library(lavaan)

test_that("pinSearch() works properly for invariant data", {
  lambda1 <- seq(.9, .7, length.out = 5)
  cov1 <- tcrossprod(lambda1) + diag(.5, 5)
  cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 5)
  dimnames(cov1) <- dimnames(cov2) <-
    list(paste0("y", 1:5), paste0("y", 1:5))
  nu1 <- seq(-0.5, 0.5, length.out = 5)
  mean1 <- nu1
  mean2 <- lambda1 * .4 + nu1
  ps1 <- pinSearch(' f =~ y1 + y2 + y3 + y4 + y5',
                   sample.cov = list(cov1, cov2),
                   sample.mean = list(mean1, mean2),
                   sample.nobs = c(100, 100))
  expect_equal(nrow(ps1[[2]]), 0)
})
