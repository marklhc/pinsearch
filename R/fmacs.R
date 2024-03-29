#' Compute fMACS effect size for two or more groups.
#'
#' `fmacs` returns the fMACS effect size statistics given a set of loadings
#'   and intercepts.
#'
#' The \eqn{f_\text{MACS}} effect size is defined as
#' \deqn{f_{\text{MACS}, i} = \frac{1}{\mathit{SD}_{iP}}
#'   \sqrt{\int [(\nu_{ij} - \bar{\nu}_j) +
#'               (\lambda_{ij} - \bar{\lambda}_j) \eta]^2 f(\eta) d \eta}}
#'   where \eqn{\lambda} is the loading and \eqn{\nu} is the intercept,
#'   and j indexes group. The effect size reflects the square root of
#'   the ratio between the variance in observed item score due to
#'   measurement noninvariance and the variance of the observed item scores.
#'   \eqn{f_\text{MACS}} is analogous to the Cohen's f effect size. When there
#'   are two groups with equal sample sizes, \eqn{f_\text{MACS}} =
#'   \eqn{f_\text{MACS}} / 2
#'
#' @param loadings A \eqn{G \times p} matrix of factor loadings, where p
#'   is the number of items and G is the number of groups.
#' @param intercepts A \eqn{G \times p} matrix of measurement intercepts.
#' @param pooled_item_sd A numeric vector of length p of the pooled standard
#'   deviation (SD) of the items across groups.
#' @param num_obs A \eqn{G \times p} matrix of sample sizes. Default assumes
#'   sample sizes are equal across groups.
#' @param latent_mean latent factor mean for the reference group. Default to 0.
#' @param latent_sd latent factor SD for the reference group. Default to 1.
#'
#' @return A 1 x p matrix of fMACS effect size.
#' @examples
#' lambda <- rbind(c(.7, .8, .7, .9),
#'                 c(.7, .8, .7, .8),
#'                 c(.8, .7, .7, .5))
#' nu <- rbind(c(0, .5, 0, 1),
#'             c(0, .2, 0, 1.1),
#'             c(0, .3, 0, 1.2))
#' fmacs(lambda,
#'       loadings = nu,
#'       pooled_item_sd = c(1, 1, 1, 1),
#'       latent_mean = 0,
#'       latent_sd = 1)
#' @export
fmacs <- function(intercepts, loadings = NULL, pooled_item_sd,
                  num_obs = 0 * intercepts + 1,
                  latent_mean = 0, latent_sd = 1) {
    total_obs <- colSums(num_obs)
    if (!is.null(loadings)) {
        mean_loading <- colSums(num_obs * loadings) / total_obs
        dloading <- sweep(loadings, MARGIN = 2, STATS = mean_loading)
    } else {
        dloading <- 0
    }
    mean_intercept <- colSums(num_obs * intercepts) / total_obs
    dintercept <- sweep(intercepts, MARGIN = 2, STATS = mean_intercept)
    integral <- colSums(num_obs * (
        dintercept^2 + 2 * dintercept * dloading * latent_mean +
            dloading^2 * (latent_sd^2 + latent_mean^2)
    )) / total_obs
    out <- sqrt(integral) / pooled_item_sd
    matrix(out, nrow = 1, dimnames = list("fmacs", colnames(loadings)))
}

#' @rdname fmacs
#' @param thresholds A matrix with G rows for measurement thresholds. The
#'     matrix must have column names indicating to which item index each column
#'     corresponds.
#' @param link Link function for the model (probit or logit).
#' @param thetas Not currently used.
#' @examples
#' # Thresholds
#' lambda <- rbind(c(.8, .5, .7, .5),
#'                 c(.8, .5, .4, .6),
#'                 c(.8, .7, .7, .5))
#' tau <- rbind(c(-0.5, 0, 1, -0.3, 0.1, 0.5, -0.5, 1.5),
#'              c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1, 1.5),
#'              c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1, 0.5))
#' # three thresholds for items 1 and 2; one threshold for items 3 and 4
#' colnames(tau) <- c(1, 1, 1, 2, 2, 2, 3, 4)
#' fmacs_ordered(tau,
#'               loadings = lambda,
#'               pooled_item_sd = c(1, 1, 1, 1),
#'               latent_mean = 0,
#'               latent_sd = 1)
#' @export
fmacs_ordered <- function(thresholds, loadings,
                          thetas = 1,
                          num_obs = 0 * loadings + 1,
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
    total_obs <- colSums(num_obs)
    for (j in seq_along(integrals)) {
        integrals[j] <- stats::integrate(
            function(x) {
                th <- thresholds[, thres_list[[j]], drop = FALSE]
                exp_y <- lapply(seq_len(nrow(th)), function(i) {
                    expected_item(loadings[i, j], cs = th[i, ], eta = x)
                })
                exp_y <- do.call(rbind, exp_y)
                mean_exp_y <- crossprod(num_obs[, j], exp_y) / total_obs[j]
                crossprod(num_obs[, j],
                          sweep(exp_y, MARGIN = 2, STATS = mean_exp_y) ^ 2) /
                    total_obs[j] *
                    stats::dnorm(x, mean = latent_mean, sd = latent_sd)
            },
            lower = -Inf, upper = Inf
        )$value
    }
    out <- matrix(sqrt(integrals) / pooled_item_sd, nrow = 1)
    colnames(out) <- colnames(loadings)
    rownames(out) <- "fmacs"
    out
}
