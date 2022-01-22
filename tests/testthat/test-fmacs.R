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
ns <- rbind(c(100, 50, 100),
            c(100, 50, 100),
            c(100, 50, 100))

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


