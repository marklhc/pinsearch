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
#' dmacs(rbind(lambdaf, lambdar),
#'       intercepts = rbind(nuf, nur),
#'       pooled_item_sd = c(1, 1, 1, 1),
#'       latent_mean = 0,
#'       latent_sd = 1)
dmacs <- function(loadings, intercepts, pooled_item_sd,
                  latent_mean = 0, latent_sd = 1) {
  dloading <- diff(loadings)
  dintercept <- diff(intercepts)
  integral <- dintercept^2 + 2 * dintercept * dloading * latent_mean +
    dloading^2 * (latent_sd^2 + latent_mean^2)
  c(sqrt(integral) / pooled_item_sd)
}

dmacs_pairwise <- function(loading_mat, intercept_mat, pooled_item_sd,
                           latent_mean = 0, latent_var = 1) {
  ngroups <- nrow(loading_mat)
  ds <- combn(
    ngroups,
    m = 2,
    FUN = function(x) {
      dmacs(
        loading_mat[x,],
        intercept_mat[x,],
        pooled_item_sd = pooled_item_sd,
        latent_mean = latent_mean,
        latent_var = latent_var
      )
    },
    simplify = FALSE
  )
  ds <- do.call(rbind, ds)
  rownames(ds) <- combn(ngroups,
                         m = 2,
                         FUN = paste,
                         collapse = " vs ")
  colnames(ds) <- colnames(loading_mat)
  ds
}
