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
#' @param num_obs A vector of length \eqn{G} of sample sizes. If not
#'   \code{NULL}, the weights will be proportional to sample sizes, assuming
#'   the same weights across items.
#' @param weights A \eqn{G \times p} matrix of weights. Default assumes
#'   equal weights across groups.
#' @param group_factor A vector of length \eqn{G} indicating grouping for 
#'   contrast. For example, `c(1, 1, 2)` means contrasting Group 1 & 2 vs. 
#'   Group 3. The default is to not combine any groups, meaning the 
#'   omnibus effect is computed.
#' @param contrast A \eqn{p \times k} contrast matrix where
#'   `colSums(contrast)` = 0. Default is `contr.sum(p)` if `group_factor` is
#'   not specified.
#' @param latent_mean latent factor mean for the reference group. Default to 0.
#' @param latent_sd latent factor SD for the reference group. Default to 1.
#'
#' @return A 1 x p matrix of fMACS effect size.
#' 
#' @importFrom stats contr.sum `contrasts<-` model.matrix
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
#' # With contrast (Group 1 & 2 vs. Group 3)
#' fmacs(lambda,
#'       loadings = nu,
#'       pooled_item_sd = c(1, 1, 1, 1),
#'       group_factor = c(1, 1, 2),
#'       latent_mean = 0,
#'       latent_sd = 1)
#' @export
fmacs <- function(intercepts, loadings = NULL, pooled_item_sd,
                  # Allow weighted vs. unweighted options
                  # Accept character and matrix input?
                  num_obs = NULL,
                  weights = 0 * intercepts + 1,
                  group_factor = NULL,
                  contrast = contr.sum(nrow(intercepts)),
                  latent_mean = 0, latent_sd = 1) {
    if (!is.null(num_obs)) {
        weights <- matrix(num_obs, nrow = nrow(intercepts),
                          ncol = ncol(intercepts))
    }
    # total_obs <- colSums(num_obs)
    weights <- sweep(weights, MARGIN = 2, STATS = colSums(weights), FUN = "/")
    if (!is.null(group_factor)) {
        stopifnot(nrow(intercepts) == nrow(loadings),
                  length(group_factor) == nrow(intercepts))
        fac <- as.factor(group_factor)
        contrasts(fac) <- contr.sum(nlevels(fac)) / as.numeric(table(fac))
        contrast <- model.matrix(~ fac)[, -1, drop = FALSE]
    }
    contrast <- as.matrix(contrast)
    if (any(colSums(contrast) != 0)) {
        warning("Contrast does not sum to 0, and results may not be correct.")
    }
    if (is.null(loadings)) {
        loadings <- 0 * intercepts
    }
    if (length(latent_mean) == 1) {
        latent_mean <- rep(latent_mean, ncol(loadings))
    }
    if (length(latent_sd) == 1) {
        latent_sd <- rep(latent_sd, ncol(loadings))
    }
    integral <- vapply(seq_len(ncol(weights)), FUN = function(j) {
        inv_ss_cont <- solve(
            crossprod(contrast, contrast / weights[, j])
        )
        cload <- crossprod(contrast, loadings[, j])
        cint <- crossprod(contrast, intercepts[, j])
        vload <- crossprod(cload, inv_ss_cont %*% cload)
        vint <- crossprod(cint, inv_ss_cont %*% cint)
        cp <- crossprod(cload, inv_ss_cont %*% cint)
        vint + 2 * cp * latent_mean[j] +
            vload * (latent_sd[j]^2 + latent_mean[j]^2)
    }, FUN.VALUE = numeric(1))
    out <- sqrt(integral) / pooled_item_sd
    out <- matrix(out, nrow = 1, dimnames = list("fmacs", colnames(loadings)))
    suppress_zero_loadings(out)
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
#' # With contrast (Group 1 & 2 vs. Group 3)
#' fmacs_ordered(tau,
#'               loadings = lambda,
#'               pooled_item_sd = c(1, 1, 1, 1),
#'               group_factor = c(1, 2, 1),
#'               latent_mean = 0,
#'               latent_sd = 1)
#' @export
fmacs_ordered <- function(thresholds, loadings,
                          thetas = 1,
                          num_obs = NULL,
                          weights = 0 * loadings + 1,
                          group_factor = NULL,
                          contrast = contr.sum(nrow(thresholds)),
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
    if (!is.null(group_factor)) {
        stopifnot(nrow(thresholds) == nrow(loadings),
                  length(group_factor) == nrow(thresholds))
        fac <- as.factor(group_factor)
        contrasts(fac) <- contr.sum(nlevels(fac)) / as.numeric(table(fac))
        contrast <- model.matrix(~ fac)[, -1, drop = FALSE]
    }
    contrast <- as.matrix(contrast)
    if (any(colSums(contrast) != 0)) {
        warning("Contrast does not sum to 0, and results may not be correct.")
    }
    if (!is.null(num_obs)) {
        weights <- matrix(num_obs, nrow = nrow(thresholds),
                          ncol = ncol(thresholds))
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
    weights <- sweep(weights, MARGIN = 2, STATS = colSums(weights), FUN = "/")
    # ww <- apply(weights, 2, function(v, g = group_factor) {
    #     ave(v, g, FUN = function(x) x / sum(x))
    # })
    for (j in seq_along(integrals)) {
        inv_ss_cont <- solve(
            crossprod(contrast, contrast / weights[, j])
        )
        integrals[j] <- stats::integrate(
            function(x) {
                th <- thresholds[, thres_list[[j]], drop = FALSE]
                exp_y <- lapply(seq_len(nrow(th)), function(i) {
                    expected_item(loadings[i, j], cs = th[i, ], eta = x)
                })
                exp_y <- do.call(rbind, exp_y)
                cexp_y <- crossprod(contrast, exp_y)
                # mean_exp_y <- crossprod(num_obs[, j], exp_y) / total_obs[j]
                # mean_exp_y <- crossprod(weights[, j], exp_y)
                # ww <- stats::ave(weights[, j], group_factor,
                #                  FUN = function(x) x / sum(x))
                # group_exp_y <-
                #     apply(ww * exp_y, 2, function(v, g = group_factor) {
                #         stats::ave(v, g, FUN = sum)
                #     })
                # crossprod(
                #     weights[, j],
                #     sweep(group_exp_y, MARGIN = 2, STATS = mean_exp_y) ^ 2) *
                #     stats::dnorm(x, mean = latent_mean, sd = latent_sd
                #     )
                colSums(cexp_y * (inv_ss_cont %*% cexp_y)) *
                    stats::dnorm(x, mean = latent_mean, sd = latent_sd)
            },
            lower = -Inf, upper = Inf
        )$value
    }
    out <- matrix(sqrt(integrals) / pooled_item_sd, nrow = 1)
    colnames(out) <- colnames(loadings)
    rownames(out) <- "fmacs"
    suppress_zero_loadings(out)
}
