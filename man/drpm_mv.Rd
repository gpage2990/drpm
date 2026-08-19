\name{drpm_mv}
\alias{drpm_mv}
\title{Fit the Hierarchical Multiview Partition Model}
\usage{
drpm_mv(
    Y,
    M0 = 1, Mv = 1,
    a_alpha = 1, b_alpha = 1, m0 = 0, s0_sq = 10,
    a0 = 2, b0 = 1, atau = 3, btau = 1,
    n_new = 1L, y_test = numeric(),
    independent = FALSE,
    c0_init = NULL, c_view_init = NULL,
    gamma_init = NULL, alpha_init = NULL,
    prior_pred_logdens_det_fn = NULL,
    do_prediction = TRUE,
    n_particles = 50L, n_pred_sweeps = 1L,
    exact_prediction = TRUE,
    marginalize_alpha_prediction = TRUE,
    n_iter = 1000L, burn = 0L, thin = 1L, )
}
\description{
Fits a Gaussian multiview restricted partition model with a latent consensus
partition. The sampler jointly updates each subject's consensus and view
labels, partially collapses the dependence parameters during the indicator
update, and restores those parameters for posterior inference.
}
\arguments{
  \item{Y}{Numeric matrix with subjects in rows and views in columns.}
  \item{M0}{Positive CRP mass parameter for the consensus partition.}
  \item{Mv}{Positive scalar or vector of view-specific CRP mass parameters.}
  \item{a_alpha, b_alpha}{Positive beta-prior shapes for the dependence parameters.}
  \item{m0, s0_sq}{Mean and variance of the prior for each view hypermean.}
  \item{a0, b0}{Positive inverse-gamma shape and scale for cluster variances.}
  \item{atau, btau}{Positive inverse-gamma shape and scale for between-cluster mean variances.}
  \item{n_new}{Number of predictive observations generated internally. Use one with the current returned \code{y_new} format.}
  \item{y_test}{Optional completely held-out vector of length \code{ncol(Y)}.}
  \item{independent}{If true, fit independent view partitions.}
  \item{c0_init}{Optional initial consensus labels.}
  \item{c_view_init}{Optional list of initial view-label vectors.}
  \item{gamma_init}{Optional list of binary dependence-indicator vectors.}
  \item{alpha_init}{Optional initial dependence probabilities.}
  \item{prior_pred_logdens_det_fn}{Optional callback replacing the built-in new-cluster prior predictive log density.}
  \item{do_prediction}{Whether to update predictive quantities.}
  \item{n_particles}{Number of auxiliary particles used for predictive simulation and diagnostic comparison.}
  \item{n_pred_sweeps}{Number of complete auxiliary-particle sweeps after each training update.}
  \item{exact_prediction}{Whether to calculate the exact joint predictive density for \code{y_test}.}
  \item{marginalize_alpha_prediction}{Whether to integrate alpha out of the predictive indicator probabilities. Posterior alpha draws are still retained.}
  \item{n_iter}{Total number of MCMC iterations.}
  \item{burn}{Number of initial iterations discarded.}
  \item{thin}{Positive thinning interval.}
}
\details{
Alpha is integrated out only while updating each training indicator vector and
is immediately restored from its beta full conditional. Exact joint prediction
sums over the held-out subject's consensus label, dependence indicators, and
view labels, avoiding predictive-particle mixing error.

For leave-one-out validation, omit one complete row, supply it as
\code{y_test}, and save \code{lp_test_exact_log_mean_exp}. Summing this scalar
over omitted rows gives the joint LOO score.
}
\value{
A list containing posterior draws and predictive summaries. The field
\code{lp_test_exact} contains one exact joint log predictive density per
retained MCMC draw. \code{lp_test_exact_log_mean_exp} is the scalar posterior
joint log predictive density for \code{y_test}. The older particle estimates
remain in \code{lp_test} and \code{lp_test_log_mean_exp} for diagnostics.
}
\examples{
set.seed(42)
z <- rep(1:2, each = 4)
Y <- cbind(
  rnorm(8, rep(c(-2, 2), each = 4), 0.4),
  rnorm(8, rep(c(-1, 3), each = 4), 0.4)
)
fit <- drpm_mv(
  Y[-1, ], n_iter = 12, burn = 4, thin = 2,
  c0_init = z[-1], c_view_init = list(z[-1], z[-1]),
  gamma_init = list(rep(1L, 7), rep(1L, 7)),
  y_test = Y[1, ], n_particles = 5
)
fit$lp_test_exact_log_mean_exp
}
