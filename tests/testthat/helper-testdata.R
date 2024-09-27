library(lavaan)
library(MASS)
library(withr, include.only = "with_seed")

simy <- function(nobs, meany, covy, thresy, seed = 2049) {
    ystar <- withr::with_seed(seed,
        MASS::mvrnorm(nobs, mu = meany, Sigma = covy)
    )
    y <- ystar
    for (j in seq_len(ncol(ystar))) {
        y[, j] <- findInterval(ystar[, j], thresy[, j])
    }
    y
}

# Initialize a binary example
num_obs <- 500
lambda1 <- seq(.9, .6, length.out = 7)
lambda2 <- c(lambda1[1], 1, lambda1[3:7])
cov1 <- tcrossprod(lambda1) + diag(.5, 7)
dimnames(cov1) <- list(paste0("yy", 1:7), paste0("yy", 1:7))
thres1 <- rbind(seq(-1.5, 1.5, length.out = 7))
thres2 <- rbind(c(thres1[1], 0.25, thres1[3:6], 0.3))
mean1 <- rep(0, 7)
y1 <- simy(num_obs, meany = mean1, covy = cov1, thresy = thres1)
cov2 <- tcrossprod(lambda1) * 1.3 + diag(.5, 7)
mean2 <- lambda1 * .4

cov3 <- tcrossprod(lambda1) * 1.3 + diag(c(rep(.5, 6), 2))

cov4 <- cov2
cov4[2, 3] <- cov4[3, 2] <- 1.2

dimnames(cov4) <- dimnames(cov3) <- dimnames(cov2) <- dimnames(cov1)

# Ordinal indicators
thres1o <- rbind(seq(-1.5, 0, length.out = 7),
                 seq(-0.5, 0.25, length.out = 7),
                 rep(1, 7))
thres2o <- rbind(c(thres1o[1, 1], -0.5, thres1o[1, 3:6], -0.5),
                 thres1o[2, ],
                 c(rep(1, 3), rep(0.5, 2), rep(1, 2)))
y1o <- simy(num_obs, meany = mean1, covy = cov1, thresy = thres1o)
y2o <- simy(num_obs, meany = mean2, covy = cov2, thresy = thres1o)
y3o <- simy(num_obs, meany = mean2, covy = cov1, thresy = thres2o)
dfo <- rbind(cbind(y1o, group = 1), cbind(y2o, group = 2),
             cbind(y3o, group = 3))