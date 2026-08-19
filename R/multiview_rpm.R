#' Fit the multiview restricted partition model
#'
#' Fits a Gaussian multiview restricted partition model with a latent consensus
#' partition. The sampler jointly updates each subject's consensus and view
#' labels, partially collapses the dependence parameters during the indicator
#' update, and restores those parameters for posterior inference. Joint
#' predictive densities can be evaluated exactly by summing over the latent
#' allocation of a completely held-out multivariate observation.
#'
#' @param Y Numeric matrix with subjects in rows and views in columns.
#' @param n_iter Total number of MCMC iterations.
#' @param burn Number of initial iterations discarded.
#' @param thin Positive thinning interval.
#' @param M0 Positive CRP mass parameter for the consensus partition.
#' @param Mv Positive scalar or vector of view-specific CRP mass parameters.
#' @param a_alpha,b_alpha Positive beta-prior shapes for the view dependence
#'   parameters.
#' @param m0,s0_sq Mean and variance of the prior for each view hypermean.
#' @param a0,b0 Positive inverse-gamma shape and scale for cluster variances.
#' @param atau,btau Positive inverse-gamma shape and scale for between-cluster
#'   mean variances.
#' @param n_new Number of posterior predictive observations generated internally.
#'   The current returned `y_new` object retains the first generated observation
#'   at each MCMC draw; use `n_new = 1` unless calling the lower-level routine.
#' @param y_test Optional numeric vector of length `ncol(Y)` containing one
#'   completely held-out multivariate observation.
#' @param independent If `TRUE`, fit independent view partitions by fixing all
#'   dependence indicators to zero.
#' @param c0_init Optional initial consensus labels of length `nrow(Y)`.
#' @param c_view_init Optional list of `ncol(Y)` initial view-label vectors.
#' @param gamma_init Optional list of binary initial dependence-indicator vectors.
#' @param alpha_init Optional initial dependence probabilities, one per view.
#' @param prior_pred_logdens_det_fn Optional advanced callback replacing the
#'   built-in new-cluster prior predictive log-density calculation.
#' @param do_prediction Whether to update predictive quantities.
#' @param n_particles Number of auxiliary particles retained for comparison and
#'   predictive simulation. Exact LOO scoring does not rely on their mixing.
#' @param n_pred_sweeps Number of complete auxiliary-particle sweeps after each
#'   training-state update.
#' @param exact_prediction Whether to compute the exact joint log predictive
#'   density for `y_test`.
#' @param marginalize_alpha_prediction Whether to integrate `alpha` out of the
#'   predictive indicator probabilities. Posterior draws of `alpha` are still
#'   retained.
#'
#' @return A list of posterior draws and predictive summaries. Of particular
#'   interest, `alpha` contains posterior dependence-parameter draws,
#'   `alpha_rb` contains their conditional posterior means,
#'   `lp_test_exact` contains one exact joint log predictive density per retained
#'   MCMC draw, and `lp_test_exact_log_mean_exp` is the scalar posterior joint
#'   log predictive density for `y_test`. The `lp_test` fields are the older
#'   auxiliary-particle estimates and are retained for diagnostics.
#'
#' @details During the update of each indicator vector, `alpha` is integrated
#'   out. It is immediately restored from its beta full conditional, so inference
#'   on dependence is retained. For leave-one-out validation, fit the function
#'   once for each omitted row and save `lp_test_exact_log_mean_exp`; summing
#'   these pointwise values gives the joint LOO score.
#'
#' @examples
#' set.seed(42)
#' z <- rep(1:2, each = 4)
#' Y <- cbind(
#'   rnorm(8, rep(c(-2, 2), each = 4), 0.4),
#'   rnorm(8, rep(c(-1, 3), each = 4), 0.4)
#' )
#' fit <- gibbs_mv_rpm(
#'   Y[-1, ], n_iter = 12, burn = 4, thin = 2,
#'   c0_init = z[-1],
#'   c_view_init = list(z[-1], z[-1]),
#'   gamma_init = list(rep(1L, 7), rep(1L, 7)),
#'   y_test = Y[1, ], n_particles = 5
#' )
#' fit$lp_test_exact_log_mean_exp
#'
#' @export
drpm_mv <- function(
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
    n_iter = 1000L, burn = 0L, thin = 1L) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (!length(Y) || nrow(Y) < 1L || ncol(Y) < 1L || any(!is.finite(Y)))
    stop("Y must be a non-empty finite numeric matrix", call. = FALSE)
  if (n_iter < 1L || burn < 0L || burn >= n_iter || thin < 1L)
    stop("require n_iter > burn >= 0 and thin >= 1", call. = FALSE)
  if (any(c(M0, Mv, a_alpha, b_alpha, s0_sq, a0, b0, atau, btau) <= 0))
    stop("all mass, variance, shape, and scale hyperparameters must be positive", call. = FALSE)
  if (!(length(Mv) %in% c(1L, ncol(Y))))
    stop("Mv must have length one or ncol(Y)", call. = FALSE)
  if (!(length(y_test) %in% c(0L, ncol(Y))) || any(!is.finite(y_test)))
    stop("y_test must be empty or a finite vector of length ncol(Y)", call. = FALSE)

  if (is.null(c_view_init)) {
    c_view_init <- lapply(seq_len(ncol(Y)), function(j) {
      k <- min(4L, nrow(Y), length(unique(Y[, j])))
      if (k <= 1L) rep.int(1L, nrow(Y))
      else stats::kmeans(Y[, j], centers = k)$cluster
    })
  }

  gibbs_mv_rpm_cpp(
    Y, as.integer(n_iter), as.integer(burn), as.integer(thin), M0,
    as.numeric(Mv), a_alpha, b_alpha, m0, s0_sq, a0, b0, atau,
    btau, as.integer(n_new), as.numeric(y_test), independent,
    c0_init, c_view_init, gamma_init, alpha_init,
    prior_pred_logdens_det_fn, do_prediction, as.integer(n_particles),
    as.integer(n_pred_sweeps), exact_prediction,
    marginalize_alpha_prediction
  )
}
