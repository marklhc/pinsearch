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
#'
#' @return A 1 x p matrix of dMACS effect size.
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
#' @export
dmacs <- function(intercepts, loadings = NULL, pooled_item_sd,
                  latent_mean = 0, latent_sd = 1) {
    if (nrow(intercepts) != 2) {
        stop("Number of rows of loadings must be 2")
    }
    if (!is.null(loadings)) {
        dloading <- diff(loadings)
    } else {
        dloading <- 0
    }
    dintercept <- diff(intercepts)
    integral <- dintercept^2 + 2 * dintercept * dloading * latent_mean +
        dloading^2 * (latent_sd^2 + latent_mean^2)
    out <- sqrt(integral) / pooled_item_sd
    # if (!is.null(rownames(loadings))) {
    #     rownames(out) <- paste(rownames(loadings), collapse = " vs ")
    # }
    rownames(out) <- "dmacs"
    colnames(out) <- colnames(loadings)
    out
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
                          latent_mean = 0, latent_sd = 1) {
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
    # item_names <- colnames(loadings)
    thres_list <- split(seq_len(ncol(thresholds)),
                        as.numeric(colnames(thresholds)))
    integrals <- rep(NA, length(thres_list))
    for (j in seq_along(integrals)) {
        integrals[j] <- stats::integrate(
            function(x) {
                th <- thresholds[, thres_list[[j]], drop = FALSE]
                (expected_item(loadings[1, j], cs = th[1, ],
                               eta = x) -
                        expected_item(loadings[2, j], cs = th[2, ],
                                      eta = x)) ^ 2 *
                    stats::dnorm(x, mean = latent_mean, sd = latent_sd)
            },
            lower = -Inf, upper = Inf
        )$value
    }
    out <- matrix(sqrt(integrals) / pooled_item_sd, nrow = 1)
    # if (!is.null(rownames(loadings))) {
    #     rownames(out) <- paste(rownames(loadings), collapse = " vs ")
    # }
    rownames(out) <- "dmacs"
    colnames(out) <- colnames(loadings)
    out
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
# Compute pooled SD
pooledvar <- function(vars, ns) {
    vars <- matrix(vars, ncol = length(ns))
    apply(vars, MARGIN = 1, stats::weighted.mean, w = ns - 1)
}

var_from_thres <- function(thres, mean = 0, sd = 1) {
    cum_ps <- stats::pnorm(thres, mean = mean, sd = sd)
    ps <- diff(c(0, cum_ps, 1))
    vals <- seq_len(length(thres) + 1) - 1
    sum(vals^2 * ps) - (sum(vals * ps))^2
}

check_inv <- function(ind, par_type, pt) {
    pt_par <- getpt(pt, type = par_type, ind_names = ind)
    npar_uniq <- nrow(unique(
        pt_par[pt_par$free != 0, c("lhs", "op", "rhs")]
    ))
    sum(pt$lhs %in% pt_par$plabel & pt$rhs %in% pt_par$plabel &
            pt$op == "==") >= nrow(pt_par) - npar_uniq -
        sum(pt_par$free == 0)
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
#' 
#' @return A matrix of 1 row showing the effect size values for
#'   each non-invariant item on each latent variable.
es_lavaan <- function(object) {
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
    ninv_ov <- ind_names[which(!inv_ind)]
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
    pars <- lavaan::lavInspect(object, what = "est")
    num_lvs <- length(object@pta$vnames$lv[[1]])
    loading_mat <- lapply(
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
    loading_mat <- do.call(rbind, loading_mat)
    sampstat <- lavaan::lavInspect(object, "sampstat")
    ns <- lavaan::lavInspect(object, "nobs")
    if (ordered) {
        # Need to think about input for thresholds to accommodate different
        # number of categories for different items
        thres_mat <- lapply(
            pars,
            FUN = function(x) {
                th_names <- rownames(x$tau)
                x$tau[grep(paste0(ninv_ov, "\\|t", collapse = "|"), th_names),]
            }
        )
        thres_mat <- do.call(rbind, thres_mat)
        thres_names <- match(gsub("\\|t.$", replacement = "",
                                  x = colnames(thres_mat)),
                             table = ninv_ov)
        thres_mat <- t(rep(1, num_lvs)) %x% thres_mat
        colnames(thres_mat) <- t(rep(1, num_lvs)) %x% thres_names
        vars <- vapply(sampstat,
                       function(x) vapply(
                           ninv_ov, function(x, y) {
                               th <- x$th[grep(paste0(y, "\\|t"), names(x$th))]
                               var_from_thres(th)
                           }, x = x, FUN.VALUE = numeric(1)
                       ), FUN.VALUE = numeric(length(ninv_ov)))
        pooled_item_sd <- sqrt(pooledvar(vars, ns))
        if (lavaan::lavInspect(object, what = "ngroups") > 2) {
            es_fun <- function(...) fmacs_ordered(..., num_obs = ns)
        } else {
            es_fun <- dmacs_ordered
        }
        es_fun(
            thresholds = thres_mat,
            loadings = loading_mat,
            latent_mean = sqrt(as.numeric(pars[[1]]$alpha)),
            latent_sd = sqrt(as.numeric(pars[[1]]$psi)),
            pooled_item_sd = rep(pooled_item_sd, num_lvs)
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
        vars <- vapply(sampstat, function(x)
            diag(x$cov)[ninv_ov],
            FUN.VALUE = numeric(length(ninv_ov)))
        pooled_item_sd <- sqrt(pooledvar(vars, ns))
        if (lavaan::lavInspect(object, what = "ngroups") > 2) {
            es_fun <- function(...) fmacs(..., num_obs = ns)
        } else {
            es_fun <- dmacs
        }
        es_fun(
            intercepts = intercept_mat,
            loadings = loading_mat,
            pooled_item_sd = rep(pooled_item_sd, num_lvs)
        )
    }
}

#' @rdname es_lavaan
#' @export
pin_es <- es_lavaan
