#' Compute dMACS effect size described in Nye & Drasgow (2011) for two groups.
#'
#' `dmacs` returns the dMACS effect size statistics given a set of loadings
#'   and intercepts.
#'
#' The \eqn{d_\text{MACS}} effect size is defined as
#'   (Nye & Drasgow, 2011, p. 968)
#' \deqn{d_{\text{MACS}, i} = \frac{1}{\mathit{SD}_{iP}}
#'   \sqrt{\int [(\nu_{iR} - \nu{iF}) +
#'               (\lambda_{iR} - \lambda_{iF}) \eta]^2 f(\eta) d \eta}}
#'   where \eqn{\lambda} is the loading and \eqn{\nu} is the intercept, F and R
#'   denote the focal and the reference group. The effect size reflects the
#'   standardized mean difference on an item due to measurement noninvariance,
#'   and is analogous to the Cohen's d effect size.
#'
#' @param loadings A \eqn{2 \times p} matrix of factor loadings, where p
#'   is the number of items.
#' @param intercepts A \eqn{2 \times p} matrix of measurement intercepts.
#' @param pooled_item_sd A numeric vector of length p of the pooled standard
#'   deviation (SD) of the items across groups.
#' @param latent_mean latent factor mean for the reference group. Default to 0.
#' @param latent_sd latent factor SD for the reference group. Default to 1.
#' @param uniqueness A vector of length \eqn{p} of uniqueness.
#' @param ns A vector of length \eqn{p} of sample sizes.
#' @param item_weights Default is `NULL`. Otherwise, one can specify a vector
#'   of length \eqn{p} of weights; if so, test-level dMACS will be computed.
#'
#' @return A 1 x p matrix of dMACS effect size. If `item_weights` is not
#'   `NULL`, \eqn{p} = 1.
#' @references Nye, C. & Drasgow, F. (2011). Effect size indices for
#'   analyses of measurement equivalence: Understanding the practical
#'   importance of differences between groups.
#'   Journal of Applied Psychology, 96(5), 966-980.
#' @examples
#' lambdaf <- c(.8, .5, .7, .5)
#' lambdar <- c(.8, .5, .4, .6)
#' nuf <- c(0.1, 0, 0.2, 0)
#' nur <- c(0.2, 0, 0, 0)
#' dmacs(rbind(nuf, nur),
#'       loadings = rbind(lambdaf, lambdar),
#'       pooled_item_sd = c(1, 1, 1, 1),
#'       latent_mean = 0,
#'       latent_sd = 1)
#' dmacs(rbind(nuf, nur),
#'       loadings = rbind(lambdaf, lambdar),
#'       pooled_item_sd = c(1, 1, 1, 1),
#'       latent_mean = 0,
#'       latent_sd = 1,
#'       item_weights = c(1, 1, 1, 1))
#' @export
dmacs <- function(intercepts, loadings = NULL,
                  pooled_item_sd = NULL,
                  latent_mean = 0, latent_sd = 1,
                  uniqueness = NULL, ns = NULL,
                  item_weights = NULL) {
    if (nrow(intercepts) != 2) {
        stop("Number of rows of loadings must be 2")
    }
    if (!is.null(item_weights)) {
        intercepts <- matrix(
            apply(intercepts, MARGIN = 1, FUN = wsum,
                  w = item_weights),
            ncol = 1)
        loadings <- matrix(
            apply(loadings, MARGIN = 1, FUN = wsum,
                  w = item_weights),
            ncol = 1)
        pooled_sd <- sqrt(wsum(pooled_item_sd^2, w = item_weights))
        cn_out <- "item_sum"
    } else {
        pooled_sd <- pooled_item_sd
        cn_out <- colnames(loadings)
    }
    if (!is.null(loadings)) {
        dloading <- diff(loadings)
    } else {
        dloading <- 0
    }
    psi <- latent_sd[1]^2
    dintercept <- diff(intercepts)
    integral <- dintercept^2 + 2 * dintercept * dloading * latent_mean +
        dloading^2 * (psi + latent_mean^2)
    if (is.null(pooled_sd) && !is.null(uniqueness)) {
        message("Pooled item SD is computed based on the input parameters")
        pooled_sd <- sqrt(
            implied_pooledvar_linear(ns, loadings, latent_sd, uniqueness)
        )
    }
    out <- sqrt(integral) / pooled_sd
    rownames(out) <- "dmacs"
    colnames(out) <- cn_out
    suppress_zero_loadings(out, loadings = loadings)
}

wsum <- function(x, w, ...) {
    weighted.mean(x, w = w, ...) * length(x)
}

#' @rdname dmacs
#' @param thresholds A matrix with two rows for measurement thresholds. The
#'     matrix must have column names indicating to which item index each column
#'     corresponds.
#' @param link Link function for the model (probit or logit).
#' @param thetas Not currently used.
#' @examples
#' # Thresholds
#' lambda <- rbind(c(.8, .5, .7, .5),
#'                 c(.8, .5, .4, .6))
#' tau <- rbind(c(-0.5, 0, 1, -0.3, 0.1, 0.5, -0.5, 1.5),
#'              c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1, 1.5))
#' # three thresholds for items 1 and 2; one threshold for items 3 and 4
#' colnames(tau) <- c(1, 1, 1, 2, 2, 2, 3, 4)
#' dmacs_ordered(tau,
#'               loadings = lambda,
#'               pooled_item_sd = c(1, 1, 1, 1),
#'               latent_mean = 0,
#'               latent_sd = 1)
#' @export
dmacs_ordered <- function(thresholds, loadings,
                          thetas = 1,
                          link = c("probit", "logit"),
                          pooled_item_sd = NULL,
                          latent_mean = 0, latent_sd = 1,
                          item_weights = NULL) {
    if (ncol(thresholds) < ncol(loadings)) {
        stop("number of thresholds should be at least the number of loadings.")
    }
    if (is.null(colnames(thresholds))) {
        stop("thresholds should have column names to indicate which thresholds",
             " are for which items")
    }
    stopifnot(nrow(thresholds) == 2, nrow(loadings) == 2)
    link <- match.arg(link)
    expected_item <- function(a, cs, eta, li = link) {
        # Check the parameterization in lavaan
        nc <- length(cs)
        pfun <- switch(li, logit = stats::plogis, probit = stats::pnorm)
        probs <- pfun(tcrossprod(rep(a, each = nc), eta) - cs)
        c(colSums(probs))
    }

    expected_item_diff <- function(j, thresholds, thres_list, loadings, eta) {
        th <- thresholds[, thres_list[[j]], drop = FALSE]
        (expected_item(loadings[1, j], cs = th[1, ], eta = eta) -
            expected_item(loadings[2, j], cs = th[2, ], eta = eta))
    }
    # item_names <- colnames(loadings)
    thres_list <- split(seq_len(ncol(thresholds)),
                        as.numeric(colnames(thresholds)))
    if (!is.null(item_weights)) {
        item_weights <- item_weights / sum(item_weights) * length(item_weights)
        pooled_sd <- sqrt(wsum(pooled_item_sd^2, w = item_weights))
        cn_out <- "item_sum"
        integrals <- stats::integrate(
            function(x) {
                out <- 0
                for (j in seq_along(thres_list)) {
                    out <- out + expected_item_diff(
                        j, thresholds, thres_list, loadings, eta = x
                    ) * item_weights[j]
                }
                out^2 * stats::dnorm(x, mean = latent_mean, sd = latent_sd)
            },
            lower = -Inf, upper = Inf
        )$value
    } else {
        pooled_sd <- pooled_item_sd
        cn_out <- colnames(loadings)
        integrals <- rep(NA, length(thres_list))
        for (j in seq_along(integrals)) {
            integrals[j] <- stats::integrate(
                function(x) {
                    expected_item_diff(
                        j, thresholds, thres_list, loadings, eta = x
                    )^2 *
                        stats::dnorm(x, mean = latent_mean, sd = latent_sd)
                },
                lower = -Inf, upper = Inf
            )$value
        }
    }
    out <- matrix(sqrt(integrals) / pooled_sd, nrow = 1)
    # if (!is.null(rownames(loadings))) {
    #     rownames(out) <- paste(rownames(loadings), collapse = " vs ")
    # }
    rownames(out) <- "dmacs"
    colnames(out) <- cn_out
    suppress_zero_loadings(out)
}

dmacs_pairwise <-
    function(loading_mat,
             intercept_mat,
             pooled_item_sd,
             latent_mean = 0,
             latent_sd = 1) {
        gp_names <- rownames(loading_mat)
        if (is.null(gp_names)) {
            gp_names <- seq_len(nrow(loading_mat))
        }
        ds <- utils::combn(
            gp_names,
            m = 2,
            FUN = function(x) {
                dmacs(
                    loading_mat[x,],
                    intercept_mat[x,],
                    pooled_item_sd = pooled_item_sd,
                    latent_mean = latent_mean,
                    latent_sd = latent_sd
                )
            },
            simplify = FALSE
        )
        ds <- do.call(rbind, ds)
        rownames(ds) <- utils::combn(gp_names,
                                     m = 2,
                                     FUN = paste,
                                     collapse = " vs ")
        colnames(ds) <- colnames(intercept_mat)
        ds
    }

# Function to extract specific parameters
getpt <- function(pt, type = c("load", "int", "thres", "uniq", "equality"),
                  ind_names) {
    # ind_names <- object@pta$vnames$ov.ind[[1]]
    rowids <- numeric(0)
    type <- match.arg(type, several.ok = TRUE)
    if ("load" %in% type) {
        rowids <- c(rowids, which(pt$op == "=~" & pt$rhs %in% ind_names))
    }
    if ("thres" %in% type) {
        rowids <- c(rowids, which(pt$op == "|" & pt$lhs %in% ind_names))
    }
    if ("int" %in% type) {
        rowids <- c(rowids, which(pt$op == "~1" & pt$lhs %in% ind_names))
    }
    if ("uniq" %in% type) {
        rowids <- c(rowids, which(pt$op == "~~" & pt$lhs %in% ind_names &
                                      pt$lhs == pt$rhs))
    }
    if ("equality" %in% type) {
        rowids <- c(rowids, which(pt$op == "=="))
    }
    pt[rowids, ]
}
# Compute pooled variance
implied_pooledvar_linear <- function(ns, ld, lat_sd, uniq) {
    ld <- as.matrix(ld)
    uniq <- as.matrix(uniq)
    vars <- ld^2 * lat_sd^2 + uniq
    pooledvar(vars, ns)
}

pooledsd_sampstat <- function(sampstat, ninv_ov, ns, ordered) {
    if (ordered) {
        vars <- lapply(sampstat,
                   function(x) vapply(
                       ninv_ov, function(x, y) {
                           th <- x$th[grep(paste0(y, "\\|t"), names(x$th))]
                           var_from_thres(th)
                       }, x = x, FUN.VALUE = numeric(1)
                   ))
    } else {
        vars <- lapply(sampstat, function(x) diag(x$cov)[ninv_ov])
    }
    vars <- do.call(rbind, vars)
    sqrt(pooledvar(vars, ns))
}

pooledvar <- function(vars, ns) {
    if (ncol(vars) == 0) {
        return(numeric(0))
    }
    w <- prop.table(as.matrix(ns), margin = 2)
    if (ncol(w) == 1) {
        w <- matrix(w, nrow = nrow(vars), ncol = ncol(vars))
    }
    # vars <- matrix(vars, ncol = length(ns))
    # apply(vars, MARGIN = 1, stats::weighted.mean, w = ns - 1)
    colSums(vars * w)
}

var_from_thres <- function(thres, mean = 0, sd = 1) {
    cum_ps <- stats::pnorm(thres, mean = mean, sd = sd)
    ps <- diff(c(0, cum_ps, 1))
    vals <- seq_len(length(thres) + 1) - 1
    sum(vals^2 * ps) - (sum(vals * ps))^2
}

check_inv <- function(ind, par_type, pt) {
    pt_par <- getpt(pt, type = par_type, ind_names = ind)
    if (any(pt_par$free == 0) &&
        (nrow(unique(pt_par[
            pt_par$free == 0,
            c("lhs", "op", "rhs")
        ])) <
            length(unique(pt_par$est[pt_par$free == 0])))) {
        return(FALSE)
    } else {
        npar_uniq <- nrow(unique(
            pt_par[pt_par$free != 0, c("lhs", "op", "rhs")]
        ))
        sum(pt$lhs %in% pt_par$plabel & pt$rhs %in% pt_par$plabel &
            pt$op == "==") >= nrow(pt_par) - npar_uniq -
            sum(pt_par$free == 0)
    }
}

suppress_zero_loadings <- function(x, loadings = NULL) {
    if (is.null(loadings)) {
        return(x)
    } else {
        return(x[, which(colSums(loadings != 0) > 0), drop = FALSE])
    }
}


to_mat_loadings <- function(pars, ninv_ov) {
    loadings <- lapply(
        pars,
        FUN = function(x) {
            sel_lambda <- x$lambda[ninv_ov, , drop = FALSE]
            dm <- dimnames(sel_lambda)
            out <- c(sel_lambda)
            names(out) <-
                outer(dm[[1]], dm[[2]], FUN = paste, sep = "-")
            out
        }
    )
    do.call(rbind, loadings)
}

to_mat_thresholds <- function(pars, ninv_ov, num_lvs) {
    thress <- lapply(
        pars,
        FUN = function(x) {
            th_names <- rownames(x$tau)
            out <- x$tau[
                grep(paste0(ninv_ov, "\\|t", collapse = "|"), th_names), ,
                drop = FALSE
            ]
            t(out)
        }
    )
    thres_mat <- do.call(rbind, thress)
    thres_names <- match(
        gsub("\\|t.$",
            replacement = "",
            x = colnames(thres_mat)
        ),
        table = ninv_ov
    )
    thres_mat <- t(rep(1, num_lvs)) %x% thres_mat
    colnames(thres_mat) <- t(rep(1, num_lvs)) %x% thres_names
    thres_mat
}

#' Item-level effect size for non-invariance
#' 
#' For two groups, the function uses [dmacs()] to compute
#'   \eqn{d_\text{MACS}}. For more than two groups, the function
#'   uses [fmacs()] to compute \eqn{f_\text{MACS}}, a
#'   generalisation of \eqn{d_\text{MACS}} similar to the 
#'   Cohen's \eqn{f} effect size.
#'
#' @param object A CFA model of class [`lavaan::lavaan-class`]
#'   fitted by [lavaan::cfa()]
#' @param ... Additional arguments passed to [dmacs()] or [fmacs()]
#' 
#' @return A matrix of 1 row showing the effect size values for
#'   each non-invariant item on each latent variable.
es_lavaan <- function(object, ...) {
    pt <- lavaan::parTable(object)
    ind_names <- object@pta$vnames$ov.ind[[1]]
    ordered <- length(lavaan::lavInspect(object, "ordered")) > 0
    if (ordered) {
        par_type <- c("load", "thres", "uniq")
    } else {
        par_type <- c("load", "int")
    }
    inv_ind <- vapply(ind_names, FUN = check_inv,
                      FUN.VALUE = logical(1),
                      par_type = par_type, pt = pt)
    pars <- lavaan::lavInspect(object, what = "est")
    if (is.null(list(...)$item_weights)) {
        ninv_ov <- ind_names[which(!inv_ind)]
        lmean <- rep(as.numeric(pars[[1]]$alpha),
                              each = length(ninv_ov))
        lsd <- rep(sqrt(as.numeric(pars[[1]]$psi)),
                            each = length(ninv_ov))
    } else {
        ninv_ov <- ind_names
        lmean <- as.numeric(pars[[1]]$alpha)
        lsd <- sqrt(as.numeric(pars[[1]]$psi))
    }
    # pt_par <- getpt(pt, type = par_type, ind_names = ind_names)
    # pt_eq <- getpt(pt, type = "equality", ind_names = ind_names)
    # ninv_par <-
    #     setdiff(pt_par$plabel, unique(c(pt_eq$lhs, pt_eq$rhs)))
    # # Remove pairs of values that are not free
    # free <- NULL
    # ninv_par <-
    #     pt_par$plabel[pt_par$plabel %in% ninv_par & pt_par$free != 0]
    # plabel <- NULL
    # ninv_ov <- unique(unlist(
    #     pt_par[pt_par$plabel %in% ninv_par, c("lhs", "rhs")]
    # ))
    # ninv_ov <- intersect(ninv_ov, ind_names)
    num_lvs <- length(object@pta$vnames$lv[[1]])
    loading_mat <- to_mat_loadings(pars, ninv_ov)
    sampstat <- lavaan::lavInspect(object, "sampstat")
    ns <- lavaan::lavInspect(object, "nobs")
    if (ordered) {
        # Need to think about input for thresholds to accommodate different
        # number of categories for different items
        thres_mat <- to_mat_thresholds(pars, ninv_ov, num_lvs)
        pooled_item_sd <- pooledsd_sampstat(
            sampstat, ninv_ov, ns, ordered = TRUE
        )
        if (lavaan::lavInspect(object, what = "ngroups") > 2) {
            es_fun <- function(...) fmacs_ordered(..., num_obs = ns)
        } else {
            es_fun <- dmacs_ordered
        }
        es_fun(
            thresholds = thres_mat,
            loadings = loading_mat,
            latent_mean = lmean,
            latent_sd = lsd,
            pooled_item_sd = rep(pooled_item_sd, num_lvs),
            ...
        )
    } else {
        intercept_mat <- lapply(
            pars,
            FUN = function(x) {
                x$nu[ninv_ov,]
            }
        )
        intercept_mat <- do.call(rbind, intercept_mat)
        intercept_mat <- t(rep(1, num_lvs)) %x% intercept_mat
        pooled_item_sd <- pooledsd_sampstat(
            sampstat, ninv_ov, ns, ordered = FALSE
        )
        if (lavaan::lavInspect(object, what = "ngroups") > 2) {
            es_fun <- function(...) fmacs(..., num_obs = ns)
        } else {
            es_fun <- dmacs
        }
        es_fun(
            intercepts = intercept_mat,
            loadings = loading_mat,
            pooled_item_sd = rep(pooled_item_sd, num_lvs),
            latent_mean = lmean,
            latent_sd = lsd,
            ...
        )
    }
}

#' @rdname es_lavaan
#' @export
pin_effsize <- es_lavaan
