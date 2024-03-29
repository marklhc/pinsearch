% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmacs.R
\name{dmacs}
\alias{dmacs}
\alias{dmacs_ordered}
\title{Compute dMACS effect size described in Nye & Drasgow (2011) for two groups.}
\usage{
dmacs(
  intercepts,
  loadings = NULL,
  pooled_item_sd,
  latent_mean = 0,
  latent_sd = 1
)

dmacs_ordered(
  thresholds,
  loadings,
  thetas = 1,
  link = c("probit", "logit"),
  pooled_item_sd = NULL,
  latent_mean = 0,
  latent_sd = 1
)
}
\arguments{
\item{intercepts}{A \eqn{2 \times p} matrix of measurement intercepts.}

\item{loadings}{A \eqn{2 \times p} matrix of factor loadings, where p
is the number of items.}

\item{pooled_item_sd}{A numeric vector of length p of the pooled standard
deviation (SD) of the items across groups.}

\item{latent_mean}{latent factor mean for the reference group. Default to 0.}

\item{latent_sd}{latent factor SD for the reference group. Default to 1.}

\item{thresholds}{A matrix with two rows for measurement thresholds. The
matrix must have column names indicating to which item index each column
corresponds.}

\item{thetas}{Not currently used.}

\item{link}{Link function for the model (probit or logit).}
}
\value{
A 1 x p matrix of dMACS effect size.
}
\description{
\code{dmacs} returns the dMACS effect size statistics given a set of loadings
and intercepts.
}
\details{
The \eqn{d_\text{MACS}} effect size is defined as
(Nye & Drasgow, 2011, p. 968)
\deqn{d_{\text{MACS}, i} = \frac{1}{\mathit{SD}_{iP}}
  \sqrt{\int [(\nu_{iR} - \nu{iF}) +
              (\lambda_{iR} - \lambda_{iF}) \eta]^2 f(\eta) d \eta}}
where \eqn{\lambda} is the loading and \eqn{\nu} is the intercept, F and R
denote the focal and the reference group. The effect size reflects the
standardized mean difference on an item due to measurement noninvariance,
and is analogous to the Cohen's d effect size.
}
\examples{
lambdaf <- c(.8, .5, .7, .5)
lambdar <- c(.8, .5, .4, .6)
nuf <- c(0.1, 0, 0.2, 0)
nur <- c(0.2, 0, 0, 0)
dmacs(rbind(nuf, nur),
      loadings = rbind(lambdaf, lambdar),
      pooled_item_sd = c(1, 1, 1, 1),
      latent_mean = 0,
      latent_sd = 1)
# Thresholds
lambda <- rbind(c(.8, .5, .7, .5),
                c(.8, .5, .4, .6))
tau <- rbind(c(-0.5, 0, 1, -0.3, 0.1, 0.5, -0.5, 1.5),
             c(-0.5, 0, 1, -0.5, 0.3, 0.5, -1, 1.5))
# three thresholds for items 1 and 2; one threshold for items 3 and 4
colnames(tau) <- c(1, 1, 1, 2, 2, 2, 3, 4)
dmacs_ordered(tau,
              loadings = lambda,
              pooled_item_sd = c(1, 1, 1, 1),
              latent_mean = 0,
              latent_sd = 1)
}
\references{
Nye, C. & Drasgow, F. (2011). Effect size indices for
analyses of measurement equivalence: Understanding the practical
importance of differences between groups.
Journal of Applied Psychology, 96(5), 966-980.
}
