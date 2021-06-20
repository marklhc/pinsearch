#' Compute dMACS effect size described in Nye & Drasgow (2011) for two groups.
#'
#' \code{dmacs} returns the dMACS effect size statistics given a set of loadings
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
#' @param loadings An \eqn{2 \times p} matrix of factor loadings, where p
#'   is the number of items.
#' @param intercepts An \eqn{2 \times p} matrix of measurement intercepts.
#' @param pooled_item_sd A numeric vector of length p of the pooled standard
#'   deviation (SD) of the items across groups.
#' @param latent_mean latent factor mean for the reference group. Default to 0.
#' @param latent_sd latent factor SD for the reference group. Default to 1.
#'
#' @return A numeric vector of length p of dMACS effect size.
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
  if (!is.null(rownames(loadings))) {
    rownames(out) <- paste(rownames(loadings), collapse = " vs ")
  }
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
    ds <- combn(
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
    rownames(ds) <- combn(gp_names,
                          m = 2,
                          FUN = paste,
                          collapse = " vs ")
    colnames(ds) <- colnames(intercept_mat)
    ds
  }

# Function to extract specific parameters
getpt <- function(pt, type = c("load", "int", "uniq", "equality"),
                  ind_names) {
  # ind_names <- object@pta$vnames$ov.ind[[1]]
  rowids <- numeric(0)
  type <- match.arg(type, several.ok = TRUE)
  if ("load" %in% type) {
    rowids <- c(rowids, which(pt$op == "=~" & pt$rhs %in% ind_names))
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
  apply(vars, MARGIN = 1, weighted.mean, w = ns - 1)
}

dmacs_lavaan <- function(object) {
  stopifnot(lavInspect(object, what = "ngroups") == 2)
  pt <- lavaan::parTable(object)
  ind_names <- object@pta$vnames$ov.ind[[1]]
  pt_par <- getpt(pt, type = c("load", "int"), ind_names = ind_names)
  pt_eq <- getpt(pt, type = "equality", ind_names = ind_names)
  ninv_par <- setdiff(pt_par$plabel, unique(c(pt_eq$lhs, pt_eq$rhs)))
  ninv_ov <-
    unique(unlist(subset(pt_par, plabel %in% ninv_par)[c("lhs", "rhs")]))
  ninv_ov <- intersect(ninv_ov, ind_names)
  pars <- lavaan::lavInspect(object, what = "est")
  sampstat <- lavInspect(object, "sampstat")
  vars <- vapply(sampstat, function(x) diag(x$cov)[ninv_ov],
                 FUN.VALUE = numeric(length(ninv_ov)))
  ns <- lavInspect(object, "nobs")
  pooled_item_sd <- sqrt(pooledvar(vars, ns))
  intercept_mat <- lapply(
    pars,
    FUN = function(x) {
      x$nu[ninv_ov, ]
    }
  )
  intercept_mat <- do.call(rbind, intercept_mat)
  loading_mat <- lapply(
    pars,
    FUN = function(x) {
      sel_lambda <- x$lambda[ninv_ov, , drop = FALSE]
      dm <- dimnames(sel_lambda)
      out <- c(sel_lambda)
      names(out) <- outer(dm[[1]], dm[[2]], FUN = paste, sep = "-")
      out
    }
  )
  loading_mat <- do.call(rbind, loading_mat)
  num_lvs <- length(object@pta$vnames$lv[[1]])
  intercept_mat <- t(rep(1, num_lvs)) %x% intercept_mat
  dmacs(intercepts = intercept_mat,
        loadings = loading_mat,
        pooled_item_sd = rep(pooled_item_sd, num_lvs))
}
