library(lavaan)
library(MASS)

# Initialize a binary example
set.seed(2049)
lambda1 <- seq(.9, .6, length.out = 7)
lambda2 <- c(lambda1[1], 1, lambda1[3:7])
cov1 <- tcrossprod(lambda1) + diag(.5, 7)
dimnames(cov1) <- list(paste0("yy", 1:7), paste0("yy", 1:7))
thres1 <- rbind(seq(-1.5, 1.5, length.out = 7))
thres2 <- rbind(c(thres1[1], 0.25, thres1[3:6], 0.3))
mean1 <- rep(0, 7)
ystar1 <- mvrnorm(500, mu = mean1, Sigma = cov1)
y1 <- ystar1
for (j in seq_len(ncol(ystar1))) {
    y1[, j] <- findInterval(ystar1[, j], thres1[, j])
}
cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 7)
dimnames(cov2) <- dimnames(cov1)
mean2 <- lambda1 * .4
ystar2 <- mvrnorm(500, mu = mean2, Sigma = cov2)
y2 <- ystar2
for (j in seq_len(ncol(ystar2))) {
    y2[, j] <- findInterval(ystar2[, j], thres2[, j])
}
df <- rbind(cbind(y1, group = 1), cbind(y2, group = 2))

test_that("pinSearch() works properly for binary data", {
  ps1 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                   data = df, group = "group", type = "thresholds",
                   ordered = paste0("yy", 1:7))
  expect_identical(ps1[[2]]$lhs, c("yy7", "yy2"))
})

test_that("pinSearch() works properly for noninvariant uniqueness", {
    # Unique variances should shift to loadings and intercepts
    cov2 <- tcrossprod(lambda1) * 1.3 + diag(c(rep(.5, 6), 2))
    dimnames(cov2) <- dimnames(cov1)
    ystar2 <- mvrnorm(500, mu = mean2, Sigma = cov2)
    y2 <- ystar2
    for (j in seq_len(ncol(ystar2))) {
        y2[, j] <- findInterval(ystar2[, j], thres1[, j])
    }
    df <- rbind(cbind(y1, group = 1), cbind(y2, group = 2))
    ps2 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                     data = df, group = "group", type = "thresholds",
                     ordered = paste0("yy", 1:7))
    expect_identical(ps2[[2]]$lhs, c("yy7"))
})

test_that("pinSearch() works properly for noninvariant unique covariances", {
    cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 7)
    cov2[2, 3] <- cov2[3, 2] <- 1.2
    dimnames(cov2) <- dimnames(cov1)
    ystar2 <- mvrnorm(500, mu = mean2, Sigma = cov2)
    y2 <- ystar2
    for (j in seq_len(ncol(ystar2))) {
        y2[, j] <- findInterval(ystar2[, j], thres2[, j])
    }
    df <- rbind(cbind(y1, group = 1), cbind(y2, group = 2))
    ps3 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                     data = df, group = "group", type = "residual.covariances",
                     ordered = paste0("yy", 1:7))
    expect_identical(ps3[[2]]$lhs, c("yy7", "yy2"))
    ps4 <- pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7
                       yy2 ~~ yy3
                       yy4 ~~ yy5 ',
                     data = df, group = "group", type = "residual.covariances",
                     ordered = paste0("yy", 1:7))
    expect_identical(ps4[[2]][3, "lhs"], "yy2")
    expect_identical(ps4[[2]][3, "rhs"], "yy3")
})

test_that("pinSearch() gives error when type = 'residuals' for ordered items", {
    expect_error(
        pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
              data = df, group = "group", type = "residuals",
              ordered = paste0("yy", 1:7))
    )
})

test_that("type = 'thresholds' gives error for continuous items", {
    expect_error(
        pinSearch(' f =~ yy1 + yy2 + yy3 + yy4 + yy5 + yy6 + yy7 ',
                  data = df, group = "group", type = "thresholds")
    )
})

# Ordinal indicators . . .
