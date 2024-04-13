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
    expect_equal(f8g[3], 0, tolerance = 1e-6)
})

test_that("Error without 'pooled_sd' argument", {
    expect_error(fmacs(nu, loadings = lambda))
    # Can compute for ordered . . .
})
