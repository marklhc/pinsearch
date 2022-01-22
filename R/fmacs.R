#' Compute fMACS effect size for two or more groups.
#'
#' \code{fmacs} returns the fMACS effect size statistics given a set of loadings
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
#' @param num_obs A \eqn{G \times p} matrix of sample sizes
#' @param latent_mean latent factor mean for the reference group. Default to 0.
#' @param latent_sd latent factor SD for the reference group. Default to 1.
#'
#' @return A numeric vector of length p of fMACS effect size.
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
    as.numeric(out)
}
