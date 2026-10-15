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

test_that("Test-level dMACS uses the variance of the item sum", {
    lambda <- rbind(c(.8, .5, .7, .5), c(.8, .5, .4, .6))
    nu <- rbind(c(0.1, 0, 0.2, 0), c(0.2, 0, 0, 0))
    # pooled variance of the (unweighted) sum score: the sum of the item
    # variances plus the common-factor covariance, pooled over groups
    cov_part <- mean((rowSums(lambda)^2 - rowSums(lambda^2)))
    # cross-check the covariance part against the definitional formula:
    # the pooled variance of the sum score is w' Sigma w with
    # Sigma_g = lambda_g lambda_g' + diag(1 - lambda_g^2) (latent SD 1
    # and item SD 1, so the uniqueness is 1 - lambda_g^2)
    sigma <- lapply(1:2, function(g) tcrossprod(lambda[g, ]) +
                      diag(1 - lambda[g, ]^2))
    w_sum <- rep(1, 4)
    expect_equal(sqrt(4 + cov_part), sqrt(mean(vapply(sigma,
        function(S) as.numeric(crossprod(w_sum, S %*% w_sum)),
        numeric(1)))))
    out <- dmacs(nu, loadings = lambda, pooled_item_sd = rep(1, 4),
                 item_weights = rep(1, 4))
    A <- diff(rowSums(nu))
    B <- diff(rowSums(lambda))
    expect_equal(out, matrix(sqrt(A^2 + B^2) / sqrt(4 + cov_part), 1,
                             dimnames = list("dmacs", "item_sum")))
    # cross-item covariances must not be omitted: the denominator is
    # larger than the square root of the sum of item variances
    expect_gt(sqrt(4 + cov_part), 2)

    # non-unit weights: quadratic weight scaling and covariances
    w <- c(2, 1, 1, 1) * 4 / 5
    out2 <- dmacs(nu, loadings = lambda, pooled_item_sd = rep(1, 4),
                  item_weights = c(2, 1, 1, 1))
    A2 <- w %*% (nu[2, ] - nu[1, ])
    B2 <- w %*% (lambda[2, ] - lambda[1, ])
    # per group: 2 * sum_{i<j} w_i w_j lambda_i lambda_j
    # = (lambda' w)^2 - sum_i w_i^2 lambda_i^2
    cov_part2 <- mean((lambda %*% w)^2 - as.numeric(lambda^2 %*% w^2))
    expect_equal(sqrt(sum(w^2) + cov_part2), sqrt(mean(vapply(sigma,
        function(S) as.numeric(crossprod(w, S %*% w)),
        numeric(1)))))
    expect_equal(out2, matrix(sqrt(A2^2 + B2^2) / sqrt(sum(w^2) + cov_part2),
                              1, dimnames = list("dmacs", "item_sum")))

    # an explicit pooled_test_sd is used as is
    out3 <- dmacs(nu, loadings = lambda, pooled_item_sd = rep(1, 4),
                  item_weights = rep(1, 4), pooled_test_sd = 2.5)
    expect_equal(out3, matrix(sqrt(A^2 + B^2) / 2.5, 1,
                              dimnames = list("dmacs", "item_sum")))

    # unequal group sizes: the covariance part is pooled with
    # group-proportional weights from ns
    ns <- c(100, 300)
    out4 <- dmacs(nu, loadings = lambda, pooled_item_sd = rep(1, 4),
                  item_weights = rep(1, 4), ns = ns)
    cov_part3 <- as.numeric((ns / sum(ns)) %*%
                              (rowSums(lambda)^2 - rowSums(lambda^2)))
    expect_equal(out4, matrix(sqrt(A^2 + B^2) / sqrt(4 + cov_part3), 1,
                              dimnames = list("dmacs", "item_sum")))

    # a vector latent_sd is reduced to its first element, as in the
    # item-level integral
    out5 <- dmacs(nu, loadings = lambda, pooled_item_sd = rep(1, 4),
                  latent_sd = c(1, 2), item_weights = rep(1, 4))
    expect_equal(out5, out)
})

test_that("Test-level effect size works without loadings", {
    nu <- rbind(c(0.3, 0.8), c(0.5, 1))
    out <- dmacs(nu, pooled_item_sd = 1, item_weights = c(1, 1))
    # No loadings: no common-factor covariance, so the denominator is the
    # square root of the sum of item variances.
    expect_equal(out, matrix(0.4 / sqrt(2), 1,
                             dimnames = list("dmacs", "item_sum")))
    out2 <- fmacs(nu, pooled_item_sd = 1, item_weights = c(1, 1))
    expect_equal(as.vector(out2), as.vector(out) / 2)
})

test_that("es_lavaan() test-level dMACS uses the pooled sum-score SD", {
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
    # The fitted model gives d(nu) = (-1, 1, 0, 0), d(lambda) = 0, with
    # item variances of 2 and covariances of 1 in both groups.
    w <- 4:1 * 4 / 10
    varT <- (sum(w))^2 + sum(w^2)
    d <- sum(w * c(-1, 1, 0, 0)) / sqrt(varT)
    expect_equal(pin_effsize(pf1, item_weights = 4:1),
                 matrix(abs(d), 1, dimnames = list("dmacs", "item_sum")),
                 tolerance = 0.001)
})

test_that("es_lavaan() accepts a user-supplied pooled_test_sd", {
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
    w <- 4:1 * 4 / 10
    # a user-supplied value takes precedence over the sample-based one
    expect_equal(
        pin_effsize(pf1, item_weights = 4:1, pooled_test_sd = 3),
        matrix(abs(sum(w * c(-1, 1, 0, 0))) / 3, 1,
               dimnames = list("dmacs", "item_sum")),
        tolerance = 0.001
    )
    # without item_weights it is accepted (and ignored)
    expect_equal(
        pin_effsize(pf1, pooled_test_sd = 3),
        pin_effsize(pf1),
        ignore_attr = TRUE
    )
})

test_that("es_lavaan() test level stays on the observed scale for ordered items", {
    set.seed(42)
    n <- 500
    # y2|t2 differs across groups (DIF); everything else is invariant.
    mkgrp <- function(th2b) {
        eta <- rnorm(n)
        u1 <- 0.9 * eta + rnorm(n, sd = sqrt(1 - 0.81))
        u2 <- 0.7 * eta + rnorm(n, sd = sqrt(1 - 0.49))
        u3 <- 0.8 * eta + rnorm(n, sd = sqrt(1 - 0.64))
        data.frame(y1 = findInterval(u1, c(-Inf, -0.5, 0.3, Inf)),
                   y2 = findInterval(u2, c(-Inf, -0.5, th2b, Inf)),
                   y3 = findInterval(u3, c(-Inf, 0, 0.8, Inf)))
    }
    d1 <- mkgrp(0.5)
    d2 <- mkgrp(0.1)
    dat <- rbind(cbind(d1, grp = 1), cbind(d2, grp = 2))
    for (v in c("y1", "y2", "y3")) {
        dat[[v]] <- factor(dat[[v]], ordered = TRUE)
    }
    fit <- cfa("f =~ y1 + y2 + y3",
               data = dat,
               group = "grp",
               ordered = c("y1", "y2", "y3"),
               group.equal = c("loadings", "thresholds"),
               group.partial = "y2|t2")
    out <- es_lavaan(fit, item_weights = c(1, 1, 1))
    # Reference: dmacs_ordered() with the fitted parameters and the
    # observed-scale variance-only fallback denominator (sampstat$cov is
    # on the latent-response scale for ordered models, so it must not be
    # used for the sum-score SD).
    pars <- lavInspect(fit, "est")
    ss <- lavInspect(fit, "sampstat")
    ns <- lavInspect(fit, "nobs")
    items <- c("y1", "y2", "y3")
    alpha1 <- if (is.null(pars[[1]]$alpha)) 0 else as.numeric(pars[[1]]$alpha)
    sv1 <- sqrt(pars[[1]]$psi["f", "f"])
    ref <- dmacs_ordered(
        thresholds = to_mat_thresholds(pars, items, 1),
        loadings = to_mat_loadings(pars, items),
        pooled_item_sd = pooledsd_sampstat(ss, items, ns, ordered = TRUE),
        latent_mean = alpha1,
        latent_sd = sv1,
        item_weights = c(1, 1, 1)
    )
    expect_equal(out, ref, tolerance = 1e-6)
    # and the latent-response-scale denominator (sampstat$cov) would give
    # a different value
    sd_latent <- sqrt(vapply(ss, function(s) {
        cov <- s$cov[items, items, drop = FALSE]
        as.numeric(crossprod(rep(1, 3), cov %*% rep(1, 3)))
    }, numeric(1)) %*% prop.table(ns))
    # Independently recover the observed-scale fallback from empirical
    # category scores. Use population variances (divisor n), matching the
    # category probabilities represented by the sample thresholds.
    group_vars <- vapply(list(d1, d2), function(d) {
        sum(vapply(d[items], function(x) mean(x^2) - mean(x)^2,
                   numeric(1)))
    }, numeric(1))
    sd_obs <- sqrt(sum(group_vars * prop.table(ns)))
    expect_equal(sd_obs,
                 sqrt(sum(pooledsd_sampstat(ss, items, ns, ordered = TRUE)^2)),
                 tolerance = 1e-6)
    # The ordered dispatch also accepts an observed-scale user override.
    expect_equal(es_lavaan(fit, item_weights = c(1, 1, 1), pooled_test_sd = 3),
                 ref * sd_obs / 3, tolerance = 1e-6)
    expect_false(isTRUE(all.equal(
        as.numeric(out),
        as.numeric(ref) * sd_obs / sd_latent
    )))
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
        mu <- if (is.null(pars[[1]]$alpha)) 0 else pars[[1]]$alpha[lv, ]
        sd2 <- pars[[1]]$psi[lv, lv]
        sqrt(dint^2 + 2 * dint * dload * mu + dload^2 * (sd2 + mu^2)) / sdy2
    }, numeric(1))
    expect_equal(as.numeric(out[1, ]), unname(expected), tolerance = 1e-6)
})

test_that("es_lavaan() works with item_weights in multidimensional models", {
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
    items <- c("y1", "y2", "y3", "y4", "b1", "b2", "b3")
    w <- rep(1, length(items))
    out <- es_lavaan(fit, item_weights = w)
    expect_identical(colnames(out), c("item_sum-a", "item_sum-b"))
    # Reference: per-latent-variable dmacs() calls with the same weights
    pars <- lavInspect(fit, "est")
    ns <- lavInspect(fit, "nobs")
    ss <- lavInspect(fit, "sampstat")
    var_mat <- do.call(cbind,
        lapply(ss, function(x) diag(x$cov)[items])
    )
    item_sd <- sqrt(var_mat %*% prop.table(ns))
    test_sd <- sqrt(vapply(ss,
        function(x) as.numeric(crossprod(w, x$cov[items, items, drop = FALSE] %*% w)),
        numeric(1)) %*% prop.table(ns))
    expected <- vapply(c("a", "b"), function(lv) {
        lam_k <- do.call(rbind, lapply(pars,
            function(x) x$lambda[items, lv, drop = TRUE]
        ))
        nu_k <- do.call(rbind, lapply(pars,
            function(x) x$nu[items, , drop = TRUE]
        ))
        dmacs(
            intercepts = nu_k,
            loadings = lam_k,
            pooled_item_sd = item_sd,
            latent_mean = if (is.null(pars[[1]]$alpha)) 0 else pars[[1]]$alpha[lv, ],
            latent_sd = sqrt(pars[[1]]$psi[lv, lv]),
            item_weights = w,
            pooled_test_sd = test_sd
        )[1, ]
    }, numeric(1))
    expect_equal(as.numeric(out[1, ]), unname(expected), tolerance = 1e-6)
    # Exercise the separate multidimensional forwarding path.
    expect_equal(es_lavaan(fit, item_weights = w, pooled_test_sd = 3),
                 out * as.numeric(test_sd) / 3, tolerance = 1e-6)
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
        latent_mean = if (is.null(pars[[1]]$alpha)) {
            rep(0, 2)
        } else {
            as.numeric(pars[[1]]$alpha)
        },
        latent_sd = sqrt(diag(pars[[1]]$psi)),
        num_obs = ns
    )
    expect_equal(as.numeric(out), as.numeric(ref), tolerance = 1e-6)
})

test_that("es_lavaan() works with item_weights (three groups)", {
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
    items <- c("y1", "y2", "y3", "y4", "b1", "b2", "b3")
    w <- rep(1, length(items))
    out <- es_lavaan(fit, item_weights = w)
    expect_identical(colnames(out), c("item_sum-a", "item_sum-b"))
    # Reference: per-latent-variable fmacs() calls with the same weights
    pars <- lavInspect(fit, "est")
    ns <- lavInspect(fit, "nobs")
    ss <- lavInspect(fit, "sampstat")
    var_mat <- do.call(cbind,
        lapply(ss, function(x) diag(x$cov)[items])
    )
    item_sd <- sqrt(var_mat %*% prop.table(ns))
    test_sd <- sqrt(vapply(ss,
        function(x) as.numeric(crossprod(w, x$cov[items, items, drop = FALSE] %*% w)),
        numeric(1)) %*% prop.table(ns))
    expected <- vapply(c("a", "b"), function(lv) {
        lam_k <- do.call(rbind, lapply(pars,
            function(x) x$lambda[items, lv, drop = TRUE]
        ))
        nu_k <- do.call(rbind, lapply(pars,
            function(x) x$nu[items, , drop = TRUE]
        ))
        fmacs(
            intercepts = nu_k,
            loadings = lam_k,
            pooled_item_sd = item_sd,
            latent_mean = if (is.null(pars[[1]]$alpha)) 0 else pars[[1]]$alpha[lv, ],
            latent_sd = sqrt(pars[[1]]$psi[lv, lv]),
            item_weights = w,
            num_obs = ns,
            pooled_test_sd = test_sd
        )[1, ]
    }, numeric(1))
    expect_equal(as.numeric(out[1, ]), unname(expected), tolerance = 1e-6)
    # The override must also reach fmacs() through the group-size wrapper.
    expect_equal(es_lavaan(fit, item_weights = w, pooled_test_sd = 3),
                 out * as.numeric(test_sd) / 3, tolerance = 1e-6)
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

test_that("Test-level ordered SD: fallback and pooled_test_sd", {
    s <- c(1.5, 1)
    # The fallback uses only the variance part of the weighted sum
    # (cross-item covariances are not identifiable from the inputs),
    # with quadratically normalized weights
    wn <- c(2, 1) * 2 / sum(c(2, 1))
    den <- sqrt(sum(wn^2 * s^2))
    d_w <- dmacs_ordered(thres, loadings = lambda, pooled_item_sd = s,
                         item_weights = c(2, 1))
    d_w_explicit <- dmacs_ordered(thres, loadings = lambda,
                                  pooled_item_sd = s, item_weights = c(2, 1),
                                  pooled_test_sd = den)
    expect_equal(as.vector(d_w), as.vector(d_w_explicit))
    # an explicit pooled_test_sd is used as is
    d_w_1 <- dmacs_ordered(thres, loadings = lambda, pooled_item_sd = s,
                           item_weights = c(2, 1), pooled_test_sd = 1)
    expect_equal(as.vector(d_w), as.vector(d_w_1) / den)
    f_w <- fmacs_ordered(thres, loadings = lambda, pooled_item_sd = s,
                         item_weights = c(2, 1))
    f_w_explicit <- fmacs_ordered(thres, loadings = lambda,
                                  pooled_item_sd = s, item_weights = c(2, 1),
                                  pooled_test_sd = den)
    expect_equal(as.vector(f_w), as.vector(f_w_explicit))
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
