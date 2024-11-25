lambda <- rbind(c(.7, .8, .7),
                c(.7, .8, .7),
                c(.8, .7, .7))
nu <- rbind(c(0, .5, 0),
            c(0, .2, 0),
            c(0, .3, 0))
tau <- rbind(c(-0.5, 0, 1, -0.3, 0.1, 0.5, -0.5),
             c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1),
             c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1))
colnames(tau) <- c(1, 1, 1, 2, 2, 2, 3)

test_that("fmacs() equals dmacs() / 2 in two groups", {
    f1 <- fmacs(nu[c(1, 3),], loadings = lambda[c(1, 3),], pooled_item_sd = 1)
    d1 <- dmacs(nu[c(1, 3),], loadings = lambda[c(1, 3),], pooled_item_sd = 1)
    f2 <- fmacs_ordered(tau[1:2,], loadings = lambda[1:2,],
                        pooled_item_sd = 1.5)
    d2 <- dmacs_ordered(tau[1:2,], loadings = lambda[1:2,],
                        pooled_item_sd = 1.5)
    expect_equal(as.vector(d1), as.vector(f1) * 2, tolerance = 0.0001)
    expect_equal(as.vector(d2), as.vector(f2) * 2, tolerance = 0.0001)
})

test_that("Noninvariant items cancelled out at test level", {
    lam <- rbind(
        c(.7, .8, .7),
        c(.7, .7, .8),
        c(.8, .7, .7)
    )
    nu <- rbind(
        c(-.5, 0, 0),
        c(0, 0, -.5),
        c(0, -.5, 0)
    )
    f1 <- fmacs(nu,
                loadings = lambda,
                pooled_item_sd = 2)
    f2 <- fmacs(
        nu,
        loadings = lambda,
        pooled_item_sd = 2,
        item_weights = c(1, 1, 1)
    )
    expect_true(f1[[1]] == f1[[2]])
    expect_equal(f2, 0, ignore_attr = TRUE)
})

test_that("fmacs() is larger with more different parameters", {
    f3 <- fmacs(nu, loadings = lambda, pooled_item_sd = 1)
    expect_equal(f3[3], 0)
    lambda2 <- rbind(lambda[1:2,], c(.5, .4, .6))
    f4 <- fmacs(nu, loadings = lambda2, pooled_item_sd = 1)
    expect_true(all(f4 > f3))
    nu2 <- rbind(c(-0.3, 0.6, 0.1), nu[2:3,])
    f5 <- fmacs(nu2, loadings = lambda2, pooled_item_sd = 1)
    expect_true(all(f5 > f4))
})

test_that("fmacs_ordered(..., group_factor) works", {
    f6 <- fmacs_ordered(tau, loadings = lambda,
                        pooled_item_sd = 1.5)
    f6g <- fmacs_ordered(tau, loadings = lambda,
                         pooled_item_sd = 1.5,
                         group_factor = c(1, 1, 2))
    expect_equal(f6[1], f6g[1])
    expect_true(all(f6g[-1] < f6[-1]))
    # f should be larger with more weights to group 3 for item 1
    f7g <- fmacs_ordered(tau, loadings = lambda,
                         pooled_item_sd = 1.5,
                         num_obs = c(100, 50, 100),
                         group_factor = c(1, 1, 2))
    expect_gt(f7g[1], f6g[1])
    # f should be larger with less weights to group 2 for item 3
    expect_gt(f7g[3], f6g[3])
    # f close to zero for item 3 with no weights for group 1
    f8g <- fmacs_ordered(tau, loadings = lambda,
                         pooled_item_sd = 1,
                         num_obs = c(1, 1e5 - 1, 1e5),
                         group_factor = c(1, 1, 2))
    expect_lt(f8g[3], 1e-3)
})

test_that("Error without 'pooled_sd' argument", {
    expect_error(fmacs(nu, loadings = lambda))
    # Can compute for ordered . . .
})

# fmacs(matrix(c(9, 5, 7, 11)), pooled_item_sd = 1,
#       num_obs = c(10, 5, 16, 9))

test_that("fmacs() works with contrast", code = {
    # Compare to results from ANOVA
    num_obs <- 5
    err <- rnorm(num_obs)
    err <- (err - mean(err))
    err <- err / sqrt(mean(err^2))
    mu2 <- c(5, 7, 9, 11, 6, 14, 8, 10)
    group <- rep(LETTERS[1:4], each = num_obs * 2)
    group2 <- rep(rep(1:2, each = num_obs), length(mu2) / 2)
    y2 <- rep(mu2, each = num_obs) + err
    aov2 <- aov(y2 ~ group * group2)
    f2_aov <- summary(aov2)[[1]]$`Sum Sq`[1:3] / sum(aov2$residuals^2)
    # Overall
    quick_f <- function(...) {
        fmacs(matrix(mu2), num_obs = rep(num_obs, length(mu2)),
              pooled_item_sd = 1, ...)
    }
    f1 <- quick_f()
    expect_equal(as.numeric(f1^2), sum(f2_aov))
    # Contrast matrix
    fac <- factor(group)
    contrasts(fac) <- contr.sum(nlevels(fac))
    fac2 <- factor(group2)
    contrasts(fac2) <- contr.sum(nlevels(fac2))
    contr <- unique(model.matrix(~ fac * fac2))
    f2 <- quick_f(contrast = contr[, -1])
    expect_equal(f1, f2)
    # Main effect for group
    f3 <- quick_f(contrast = contr[, 2:4])
    f3g <- quick_f(group_factor = c(1, 1, 2, 2, 3, 3, 4, 4))
    expect_equal(f3, f3g)
    # Main effect for group
    f4 <- quick_f(contrast = contr[, 5])
    f4g <- quick_f(group_factor = c(1, 2, 1, 2, 1, 2, 1, 2))
    expect_equal(f4, f4g)
    # Interaction
    f5 <- quick_f(contrast = contr[, 6:8])
    expect_equal(as.numeric(c(f3, f4, f5)^2), f2_aov)
    expect_equal(sum(c(f3, f4, f5)^2), as.numeric(f1^2))
})

test_that("fmacs_ordered() works with contrast", code = {
    lambda <- rbind(
        c(.7, .8, .8),
        c(.7, .8, .7),
        c(.8, .7, .7),
        c(.5, .3, .8)
    )
    tau <- rbind(
        c(-0.5, 0, 1, -0.3, 0.1, 0.5, 0),
        c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1),
        c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1),
        c(-0.5, 0, 1, -0.5, 0.3, 0.5, 0)
    )
    colnames(tau) <- c(1, 1, 1, 2, 2, 2, 3)
    f9 <- fmacs_ordered(tau, loadings = lambda,
                        pooled_item_sd = 1.5)
    f10 <- fmacs_ordered(tau, loadings = lambda,
                         pooled_item_sd = 1.5,
                         contrast = c(1, 1, -1, -1))
    expect_true(all(f9 > f10))
    expect_equal(f10[3], 0)
})

test_that(
    "Error when `group_factor` has incorrect length",
    code = {
        expect_error(
            fmacs_ordered(tau,
                loadings = lambda,
                pooled_item_sd = 1.5,
                group_factor = c(1, 1, 2, 3)
            )
        )
    }
)