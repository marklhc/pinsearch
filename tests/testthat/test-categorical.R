library(lavaan)
library(MASS)

# Initialize a binary example
lambda1 <- seq(.9, .6, length.out = 7)
lambda2 <- c(lambda1[1], 1, lambda1[3:7])
cov1 <- tcrossprod(lambda1) + diag(.5, 7)
dimnames(cov1) <- list(paste0("yy", 1:7), paste0("yy", 1:7))
thres1 <- rbind(seq(-1.5, 1.5, length.out = 7))
thres2 <- rbind(c(thres1[1], 0.25, thres1[3:6], 0.3))
mean1 <- rep(0, 7)
ystar1 <- mvrnorm(200, mu = mean1, Sigma = cov1)
y1 <- ystar1
for (j in seq_len(ncol(ystar1))) {
    y1[, j] <- findInterval(ystar1[, j], thres1[, j])
}

test_that("pinSearch() works properly for invariant data", {
  cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 7)
  dimnames(cov2) <- dimnames(cov2)
  mean2 <- lambda1 * .4
  ystar2 <- mvrnorm(200, mu = mean2, Sigma = cov2)
  y2 <- ystar2
  for (j in seq_len(ncol(ystar2))) {
      y2[, j] <- findInterval(ystar2[, j], thres2[, j])
  }
  df <- rbind(cbind(y1, group = 1), cbind(y2, group = 2))
  ps1 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                   data = df, group = "group", type = "thresholds",
                   ordered = paste0("yy", 1:7))
  expect_equal(nrow(ps1[[2]]), 0)
})
