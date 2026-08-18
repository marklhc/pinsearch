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

test_that("`pin_effsize()` invariant with scaling", {
    cov2 <- tcrossprod(lambda2) * 1.3 + diag(.5, 5)
    dimnames(cov2) <- list(paste0("y", 1:5), paste0("y", 1:5))
    mean2 <- lambda2 * .4 + nu2
    pf1 <- cfa(' f =~ y1 + y2 + y3 + y4 + y5',
               sample.cov = list(cov1, cov2),
               sample.mean = list(mean1, mean2),
               sample.nobs = c(10000, 10000),
               std.lv = TRUE,
               group.equal = c("loadings", "intercepts"),
               group.partial = c("f=~y1", "f=~y2",
                                 "y1~1", "y2~1", "y4~1", "y5~1"))
    pf2 <- cfa(' f =~ NA * y1 + y2 + y3 + y4 + y5
                 f ~~ c(1.5, NA) * f
                 f ~ c(-0.5, NA) * 1 ',
               sample.cov = list(cov1, cov2),
               sample.mean = list(mean1, mean2),
               sample.nobs = c(10000, 10000),
               group.equal = c("loadings", "intercepts"),
               group.partial = c("f=~y1", "f=~y2",
                                 "y1~1", "y2~1", "y4~1", "y5~1"))
    expect_equal(pin_effsize(pf1),
                 pin_effsize(pf2),
                 tolerance = 0.00001)
})

test_that("`pin_effsize()` works for scaling indicator", {
    pf1 <- cfa(' f =~ c(1, .9) * x1 + x2 + x3 ',
               data = HolzingerSwineford1939,
               std.lv = TRUE,
               group = "school",
               group.equal = c("loadings", "intercepts"),
               group.partial = c("f=~x1"))
    pin_es1 <- pin_effsize(pf1)
    expect_equal(dim(pin_es1), c(1, 1))
})

test_that("Noninvariant items cancelled out at test level", {
    lam <- rep(1, 4)
    cov1 <- tcrossprod(lam) + diag(1, 4)
    dimnames(cov1) <- list(paste0("x", 1:4), paste0("x", 1:4))
    mean1 <- c(1, 0, 0, 0)
    mean2 <- lam * .4 + c(0, 1, 0, 0)
    pf1 <- cfa(' f =~ x1 + x2 + x3 + x4 ',
               sample.cov = list(cov1, cov1),
               sample.mean = list(mean1, mean2),
               sample.nobs = c(10000, 10000),
               std.lv = TRUE,
               group = "school",
               group.equal = c("loadings", "intercepts"),
               group.partial = c("x1~1", "x2~1"))
    expect_true(all(pin_effsize(pf1) > 0.7))
    expect_equal(
        pin_effsize(pf1, item_weights = c(1, 1, 1, 1)),
        0,
        ignore_attr = TRUE
    )
    expect_gt(
        pin_effsize(pf1, item_weights = 4:1), 0
    )
    expect_error(
        pin_effsize(pf1, item_weights = rep(1, 5))
    )
})

# Multidimensional models (see #7, #4)
# Build a population (sample.cov / sample.mean) pair for the two-factor
# example used in the multidimensional es_lavaan tests. Group 2 (and 3)
# have lower loadings of y2 on both factors and a higher intercept.
# y4 and b3 are invariant extra indicators for identification headroom.
mk_md_covs <- function(ngroups = 2) {
    lv <- c("a", "b")
    ov <- c("y1", "y2", "y3", "y4", "b1", "b2", "b3")
    psi <- matrix(c(1, -.4, -.4, 1.25), 2, dimnames = list(lv, lv))
    # Latent means are 0: the fitted model cannot estimate them (no mean
    # structure in md_model), so keep them 0 in the data generation as well.
    alpha <- c(a = 0, b = 0)
    nu1 <- c(y1 = 0, y2 = .1, y3 = .2, y4 = .15, b1 = 0, b2 = .1, b3 = .25)
    L1 <- matrix(0, 7, 2, dimnames = list(ov, lv))
    L1[c("y1", "y2", "y3", "y4"), "a"] <- c(1, .9, .7, .65)
    # b~b1 (and a~y1) are the reference loadings, fixed at 1 in the fit.
    L1[c("b1", "b2", "b3"), "b"] <- c(1, .8, .75)
    L1["y2", "b"] <- .8
    mkgroup <- function(L, nu) {
        cov <- as.matrix(L %*% psi %*% t(L) + diag(1, 7))
        dimnames(cov) <- list(ov, ov)
        list(cov = cov, mean = as.numeric(L %*% alpha + nu))
    }
    g1 <- mkgroup(L1, nu1)
    out <- list(cov = list(g1$cov), mean = list(g1$mean), nobs = 300)
    if (ngroups > 1) {
        L2 <- L1
        L2["y2", "a"] <- .6
        L2["y2", "b"] <- .3
        nu2 <- nu1
        nu2["y2"] <- .3
        g2 <- mkgroup(L2, nu2)
        out$cov <- c(out$cov, list(g2$cov))
        out$mean <- c(out$mean, list(g2$mean))
        out$nobs <- c(out$nobs, 300)
    }
    if (ngroups > 2) {
        L3 <- L1
        L3["y2", "a"] <- .85
        L3["y2", "b"] <- .5
        nu3 <- nu1
        nu3["y2"] <- .2
        g3 <- mkgroup(L3, nu3)
        out$cov <- c(out$cov, list(g3$cov))
        out$mean <- c(out$mean, list(g3$mean))
        out$nobs <- c(out$nobs, 300)
    }
    out
}

# No latent variable mean structure: with group.equal, free latent means
# trigger a spurious "vcov not positive definite" warning from lavaan.
# Latent means default to 0 in this case.
md_model <- "a =~ y1 + y2 + y3 + y4
             b =~ b1 + b2 + b3 + y2"

test_that("es_lavaan() computes correct values in multidimensional models (two groups)", {
    dat <- mk_md_covs(2)
    fit0 <- cfa(md_model,
        sample.cov = dat$cov,
        sample.mean = dat$mean,
        sample.nobs = dat$nobs,
        group.equal = c("loadings", "intercepts")
    )
    pt0 <- parTable(fit0)
    pt0$id <- seq_len(nrow(pt0))
    pt0 <- remove_cons(pt0, lhs = "a", rhs = "y2", group = 2, op = "=~")
    pt0 <- remove_cons(pt0, lhs = "b", rhs = "y2", group = 2, op = "=~")
    pt0 <- remove_cons(pt0, lhs = "y2", rhs = "", group = 2, op = "~1")
    fit <- cfa(pt0,
        sample.cov = dat$cov,
        sample.mean = dat$mean,
        sample.nobs = dat$nobs
    )
    out <- expect_no_warning(es_lavaan(fit))
    expect_identical(colnames(out), c("y2-a", "y2-b"))
    # Expected values computed from the fitted estimates
    pars <- lavInspect(fit, "est")
    ns <- lavInspect(fit, "nobs")
    ss <- lavInspect(fit, "sampstat")
    sdy2 <- sqrt(sum(vapply(ss,
        function(x) x$cov["y2", "y2"], numeric(1)
    ) * prop.table(ns)))
    dint <- pars[[2]]$nu["y2", ] - pars[[1]]$nu["y2", ]
    expected <- vapply(c("a", "b"), function(lv) {
        dload <- pars[[2]]$lambda["y2", lv] - pars[[1]]$lambda["y2", lv]
        mu <- pars[[1]]$alpha[lv, ]
        sd2 <- pars[[1]]$psi[lv, lv]
        sqrt(dint^2 + 2 * dint * dload * mu + dload^2 * (sd2 + mu^2)) / sdy2
    }, numeric(1))
    expect_equal(as.numeric(out[1, ]), unname(expected), tolerance = 1e-6)
})

test_that("es_lavaan() computes correct values in multidimensional models (three groups)", {
    dat <- mk_md_covs(3)
    fit0 <- cfa(md_model,
        sample.cov = dat$cov,
        sample.mean = dat$mean,
        sample.nobs = dat$nobs,
        group.equal = c("loadings", "intercepts")
    )
    pt0 <- parTable(fit0)
    pt0$id <- seq_len(nrow(pt0))
    pt0 <- remove_cons(pt0, lhs = "a", rhs = "y2", group = 3, op = "=~")
    pt0 <- remove_cons(pt0, lhs = "b", rhs = "y2", group = 3, op = "=~")
    pt0 <- remove_cons(pt0, lhs = "y2", rhs = "", group = 3, op = "~1")
    pt0 <- remove_cons(pt0, lhs = "a", rhs = "y2", group = 2, op = "=~")
    pt0 <- remove_cons(pt0, lhs = "b", rhs = "y2", group = 2, op = "=~")
    pt0 <- remove_cons(pt0, lhs = "y2", rhs = "", group = 2, op = "~1")
    fit <- cfa(pt0,
        sample.cov = dat$cov,
        sample.mean = dat$mean,
        sample.nobs = dat$nobs
    )
    out <- expect_no_warning(es_lavaan(fit))
    expect_identical(colnames(out), c("y2-a", "y2-b"))
    # Reference call of fmacs() with the correct latent parameters
    pars <- lavInspect(fit, "est")
    lam <- do.call(rbind, lapply(pars, function(x) x$lambda["y2", ]))
    nu <- do.call(rbind, lapply(pars, function(x) x$nu["y2", ]))
    ns <- lavInspect(fit, "nobs")
    ss <- lavInspect(fit, "sampstat")
    sdy2 <- sqrt(sum(vapply(ss,
        function(x) x$cov["y2", "y2"], numeric(1)
    ) * prop.table(ns)))
    ref <- fmacs(intercepts = t(rep(1, 2)) %x% nu,
        loadings = lam,
        pooled_item_sd = rep(sdy2, 2),
        latent_mean = as.numeric(pars[[1]]$alpha),
        latent_sd = sqrt(diag(pars[[1]]$psi)),
        num_obs = ns
    )
    expect_equal(as.numeric(out), as.numeric(ref), tolerance = 1e-6)
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

test_that("Error without 'pooled_sd' argument", {
    expect_error(dmacs(rbind(c(0.3, 0.8), c(0.5, 1))))
    # Can compute for ordered . . .
})

test_that("Pooled SD computed correctly using uniqueness", {
    cov2 <- tcrossprod(lambda2) * 1.3 + diag(.5, 5)
    dimnames(cov2) <- list(paste0("y", 1:5), paste0("y", 1:5))
    mean2 <- lambda2 * .4 + nu2
    pscalar1 <- cfa(' f =~ y1 + y2 + y3 + y4 + y5',
                    sample.cov = list(cov1, cov2),
                    sample.mean = list(mean1, mean2),
                    sample.nobs = c(1000000, 1000000),
                    group.equal = c("loadings", "intercepts"),
                    group.partial = c("f=~y1", "f=~y2",
                                      "y1~1", "y2~1", "y4~1", "y5~1"),
                    std.lv = TRUE)
    expect_equal(
        as.numeric(es_lavaan(pscalar1)),
        dmacs(rbind(nu1, nu2), loadings = rbind(lambda1, lambda2),
              uniqueness = matrix(.5, nrow = 2, ncol = 5),
              latent_sd = c(1, sqrt(1.3)),
              ns = c(100, 100))[c(1:2, 4:5)],
        tolerance = 0.0001
    )
})
