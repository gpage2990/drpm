#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
using namespace Rcpp;

// fixed11 implements the partition full conditionals in the supplied model:
// candidates are constructed in reduced canonical coordinates, compatibility
// is checked candidate by candidate (including pegged observations), and a
// selected new cluster receives an atom from its posterior given y.

// Numerical version of the R prior_pred_logdens_det() calculation.  The
// optional R callback remains supported, but fixed11 is self-contained.
struct PriorPredParams {
  double y;
  double m0;
  double s0_sq;
  double a0;
  double b0;
};

static void prior_pred_integrand_cpp(double* s2, int n, void* data) {
  const PriorPredParams* p = static_cast<PriorPredParams*>(data);
  const double log_const = p->a0 * std::log(p->b0) - R::lgammafn(p->a0);
  for (int i = 0; i < n; ++i) {
    const double x = s2[i];
    if (!(x > 0.0) || !R_finite(x)) {
      s2[i] = 0.0;
      continue;
    }
    const double log_value =
      R::dnorm4(p->y, p->m0, std::sqrt(p->s0_sq + x), 1) +
      log_const - (p->a0 + 1.0) * std::log(x) - p->b0 / x;
    s2[i] = std::exp(log_value);
  }
}

static double prior_pred_logdens_det_cpp(double y,
                                         double m0,
                                         double s0_sq,
                                         double a0,
                                         double b0) {
  if (!(s0_sq > 0.0) || !(a0 > 0.0) || !(b0 > 0.0)) return R_NegInf;
  PriorPredParams params = {y, m0, s0_sq, a0, b0};
  double bound = 0.0, epsabs = 1e-8, epsrel = 1e-8;
  double result = 0.0, abserr = 0.0;
  int inf = 1, neval = 0, ier = 0, limit = 100, lenw = 400, last = 0;
  std::vector<int> iwork(limit);
  std::vector<double> work(lenw);
  Rdqagi(prior_pred_integrand_cpp, &params, &bound, &inf,
         &epsabs, &epsrel, &result, &abserr, &neval, &ier,
         &limit, &lenw, &last, iwork.data(), work.data());
  if (ier != 0 || !(result > 0.0) || !R_finite(result)) {
    stop("prior_pred_logdens_det_cpp: numerical integration failed");
  }
  return std::log(result);
}

static double prior_pred_logdens_dispatch(const Nullable<Function>& fn,
                                          double y,
                                          double m0,
                                          double s0_sq,
                                          double a0,
                                          double b0) {
  if (fn.isNotNull()) {
    Function callback(fn.get());
    return as<double>(callback(y, m0, s0_sq, a0, b0));
  }
  return prior_pred_logdens_det_cpp(y, m0, s0_sq, a0, b0);
}

// Exact draw from p(mu,sigma2 | y) under the independent base measure
// mu ~ N(mu0,tau2), sigma2 ~ IG(a0,b0).  Rejection sampling sigma2 from
// its IG prior is valid because N(y | mu0,tau2+sigma2) is bounded above
// by N(mu0 | mu0,tau2); mu is then drawn from its Gaussian conditional.
static NumericVector draw_new_atom_posterior_cpp(double y,
                                                  double mu0,
                                                  double tau2,
                                                  double a0,
                                                  double b0) {
  if (!(tau2 > 0.0) || !(a0 > 0.0) || !(b0 > 0.0)) {
    stop("draw_new_atom_posterior_cpp: invalid hyperparameters");
  }
  const double log_bound = -0.5 * std::log(2.0 * M_PI * tau2);
  double s2 = NA_REAL;
  bool accepted = false;
  for (int attempt = 0; attempt < 1000000; ++attempt) {
    const double proposal = 1.0 / R::rgamma(a0, 1.0 / b0);
    const double log_accept =
      R::dnorm4(y, mu0, std::sqrt(tau2 + proposal), 1) - log_bound;
    if (std::log(R::runif(0.0, 1.0)) <= log_accept) {
      s2 = proposal;
      accepted = true;
      break;
    }
  }
  if (!accepted) stop("draw_new_atom_posterior_cpp: rejection sampler failed");

  const double v = 1.0 / (1.0 / tau2 + 1.0 / s2);
  const double m = v * (mu0 / tau2 + y / s2);
  return NumericVector::create(R::rnorm(m, std::sqrt(v)), s2);
}

// -----------------------------------------------------------------------------
// Shared helpers
// -----------------------------------------------------------------------------

static inline IntegerVector canon_cpp(const IntegerVector& z_in) {
  const int n = z_in.size();
  IntegerVector out(n);

  int next_label = 0;
  std::unordered_map<int, int> map;
  map.reserve(static_cast<size_t>(n) * 2 + 1);

  for (int i = 0; i < n; ++i) {
    const int lab = z_in[i];
    if (lab == NA_INTEGER) {
      stop("canon_cpp: NA labels are not allowed");
    }
    auto it = map.find(lab);
    if (it == map.end()) {
      ++next_label;
      map.emplace(lab, next_label);
      out[i] = next_label;
    } else {
      out[i] = it->second;
    }
  }
  return out;
}

static inline IntegerVector first_occurrence_labels_cpp(const IntegerVector& z_in) {
  const int n = z_in.size();
  std::unordered_set<int> seen;
  seen.reserve(static_cast<size_t>(n) * 2 + 1);
  std::vector<int> labs;
  labs.reserve(static_cast<size_t>(n));

  for (int i = 0; i < n; ++i) {
    const int lab = z_in[i];
    if (lab == NA_INTEGER) {
      stop("first_occurrence_labels_cpp: NA labels are not allowed");
    }
    if (seen.insert(lab).second) labs.push_back(lab);
  }

  IntegerVector out(labs.size());
  for (size_t k = 0; k < labs.size(); ++k) out[k] = labs[k];
  return out;
}

static inline IntegerVector subset_by_positions(const IntegerVector& x,
                                                const std::vector<int>& pos0) {
  IntegerVector out(pos0.size());
  for (size_t k = 0; k < pos0.size(); ++k) {
    out[k] = x[pos0[k]];
  }
  return out;
}

static inline IntegerVector insert_label_cpp(const IntegerVector& reduced,
                                             int position0,
                                             int label) {
  if (position0 < 0 || position0 > reduced.size()) {
    stop("insert_label_cpp: invalid insertion position");
  }
  IntegerVector out(reduced.size() + 1);
  for (int src = 0, dst = 0; dst < out.size(); ++dst) {
    if (dst == position0) out[dst] = label;
    else out[dst] = reduced[src++];
  }
  return out;
}

static inline int reduced_label_for_old_cluster_cpp(const IntegerVector& used_old_labels,
                                                     int old_label,
                                                     int new_label) {
  for (int k = 0; k < used_old_labels.size(); ++k) {
    if (used_old_labels[k] == old_label) return k + 1;
  }
  return new_label;
}

static inline IntegerVector counts_from_labels(const IntegerVector& labels) {
  if (labels.size() == 0) return IntegerVector(0);
  int maxlab = 0;
  for (int i = 0; i < labels.size(); ++i) {
    if (labels[i] > maxlab) maxlab = labels[i];
  }
  IntegerVector counts(maxlab);
  for (int i = 0; i < labels.size(); ++i) {
    counts[labels[i] - 1]++;
  }
  return counts;
}

bool same_partition_cpp(IntegerVector a, IntegerVector b) {
  const int n = a.size();
  if (n != b.size()) return false;
  if (n <= 1) return true;

  int max_a = 0, max_b = 0;
  for (int i = 0; i < n; ++i) {
    if (IntegerVector::is_na(a[i]) || IntegerVector::is_na(b[i])) return false;
    if (a[i] <= 0 || b[i] <= 0) return false;
    if (a[i] > max_a) max_a = a[i];
    if (b[i] > max_b) max_b = b[i];
  }

  std::vector<int> map_ab(max_a + 1, 0);
  std::vector<int> map_ba(max_b + 1, 0);

  for (int i = 0; i < n; ++i) {
    const int ak = a[i];
    const int bk = b[i];
    const int mab = map_ab[ak];
    const int mba = map_ba[bk];

    if (mab == 0 && mba == 0) {
      map_ab[ak] = bk;
      map_ba[bk] = ak;
    } else if (mab != bk || mba != ak) {
      return false;
    }
  }

  return true;
}

double log_crp_from_sizes_cpp(IntegerVector sizes, double M) {
  if (M <= 0.0) return R_NegInf;

  double n_d = 0.0;
  int K = 0;
  for (int i = 0; i < sizes.size(); ++i) {
    const int s = sizes[i];
    if (s > 0) {
      n_d += s;
      ++K;
    }
  }

  const int n = static_cast<int>(n_d);
  if (n == 0) return 0.0;

  double out = K * std::log(M);
  for (int i = 0; i < n; ++i) {
    out -= std::log(M + i);
  }
  for (int i = 0; i < sizes.size(); ++i) {
    const int s = sizes[i];
    if (s > 0) out += R::lgammafn(static_cast<double>(s));
  }
  return out;
}

IntegerVector sample_logw_cpp(NumericVector logw, int fallback = 1) {
  const int n = logw.size();
  if (n == 0) return IntegerVector::create(fallback);

  double mx = R_NegInf;
  bool any_finite = false;
  for (int i = 0; i < n; ++i) {
    if (R_finite(logw[i])) {
      any_finite = true;
      if (logw[i] > mx) mx = logw[i];
    }
  }
  if (!any_finite) return IntegerVector::create(fallback);

  NumericVector w(n);
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    if (R_finite(logw[i])) w[i] = std::exp(logw[i] - mx);
    else w[i] = 0.0;
    s += w[i];
  }
  if (!(s > 0.0) || !R_finite(s)) return IntegerVector::create(fallback);

  const double u = R::runif(0.0, s);
  double cum = 0.0;
  for (int i = 0; i < n; ++i) {
    cum += w[i];
    if (u <= cum) return IntegerVector::create(i + 1);
  }
  return IntegerVector::create(n);
}

static inline double log_sum_exp_core(const NumericVector& x) {
  double mx = R_NegInf;
  bool any_finite = false;
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) {
      any_finite = true;
      if (x[i] > mx) mx = x[i];
    }
  }
  if (!any_finite) return R_NegInf;
  double s = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) s += std::exp(x[i] - mx);
  }
  return mx + std::log(s);
}

static inline List relabel_partition_with_atoms_cpp(const IntegerVector& z_in,
                                                    const NumericVector& mu_in,
                                                    const NumericVector& sigma2_in) {
  IntegerVector oldlabs = first_occurrence_labels_cpp(z_in);
  IntegerVector z_new = canon_cpp(z_in);

  NumericVector mu_new(oldlabs.size());
  NumericVector sigma2_new(oldlabs.size());
  for (int k = 0; k < oldlabs.size(); ++k) {
    const int lab0 = oldlabs[k] - 1;
    if (lab0 < 0 || lab0 >= mu_in.size() || lab0 >= sigma2_in.size()) {
      stop("relabel_partition_with_atoms_cpp: label index out of bounds");
    }
    mu_new[k] = mu_in[lab0];
    sigma2_new[k] = sigma2_in[lab0];
  }

  return List::create(Named("z") = z_new,
                      Named("mu") = mu_new,
                      Named("sigma2") = sigma2_new);
}

// -----------------------------------------------------------------------------
// c0 update kernel
// -----------------------------------------------------------------------------

IntegerVector update_c0_cpp(IntegerVector c0,
                            List c_view,
                            List gamma,
                            double M0,
                            NumericVector Mv,
                            IntegerVector ord) {
  const int m = c0.size();
  const int J = c_view.size();
  if (gamma.size() != J) stop("update_c0_cpp: c_view and gamma must have same length");
  if (ord.size() == 0) stop("update_c0_cpp: ord must be non-empty");
  if (!(M0 > 0.0)) stop("update_c0_cpp: M0 must be positive");
  if (!(Mv.size() == 1 || Mv.size() == J)) stop("update_c0_cpp: invalid Mv length");

  for (int o = 0; o < ord.size(); ++o) {
    const int i0 = ord[o] - 1;
    if (i0 < 0 || i0 >= m) stop("update_c0_cpp: ord contains invalid index");
    const int old_lab = c0[i0];
    IntegerVector raw_minus(m - 1);
    for (int src = 0, dst = 0; src < m; ++src) if (src != i0) raw_minus[dst++] = c0[src];
    IntegerVector used = first_occurrence_labels_cpp(raw_minus);
    IntegerVector c0_red = canon_cpp(raw_minus);
    const int K = used.size();
    IntegerVector counts = counts_from_labels(c0_red);
    const int fallback_h = reduced_label_for_old_cluster_cpp(used, old_lab, K + 1);
    NumericVector logw(K + 1, R_NegInf);
    std::vector<IntegerVector> candidates;
    candidates.reserve(K + 1);
    for (int h = 1; h <= K + 1; ++h) {
      IntegerVector candidate = canon_cpp(insert_label_cpp(c0_red, i0, h));
      candidates.push_back(candidate);
      double score = (h <= K) ? std::log(static_cast<double>(counts[h - 1])) : std::log(M0);
      bool compatible = true;
      for (int j = 0; j < J; ++j) {
        IntegerVector cj = as<IntegerVector>(c_view[j]);
        IntegerVector gj = as<IntegerVector>(gamma[j]);
        if (cj.size() != m || gj.size() != m) stop("update_c0_cpp: invalid view state");
        std::vector<int> fixed;
        for (int ii = 0; ii < m; ++ii) if (gj[ii] == 1) fixed.push_back(ii);
        IntegerVector cj_fixed = subset_by_positions(cj, fixed);
        IntegerVector c0_fixed = subset_by_positions(candidate, fixed);
        if (!same_partition_cpp(cj_fixed, c0_fixed)) {
          compatible = false;
          break;
        }
        const double Mj = (Mv.size() == 1) ? Mv[0] : Mv[j];
        score += log_crp_from_sizes_cpp(counts_from_labels(cj), Mj) -
                 log_crp_from_sizes_cpp(counts_from_labels(cj_fixed), Mj);
      }
      if (compatible) logw[h - 1] = score;
    }
    const int chosen_h = sample_logw_cpp(logw, fallback_h)[0];
    c0 = clone(candidates[chosen_h - 1]);
  }
  return c0;
}

// -----------------------------------------------------------------------------
// c_view update kernel
// -----------------------------------------------------------------------------

List update_cview_cpp(NumericMatrix Y,
                      List c_view,
                      List mu,
                      List sigma2,
                      IntegerVector c0,
                      List gamma,
                      NumericVector mu0,
                      NumericVector tau2,
                      double a0,
                      double b0,
                      NumericVector Mv,
                      Nullable<Function> prior_pred_logdens_det_fn = R_NilValue) {
  const int m = Y.nrow();
  const int J = Y.ncol();

  if (c_view.size() != J || mu.size() != J || sigma2.size() != J || gamma.size() != J) {
    stop("update_cview_cpp: c_view, mu, sigma2, and gamma must all have length ncol(Y)");
  }
  if (mu0.size() != J || tau2.size() != J) {
    stop("update_cview_cpp: mu0 and tau2 must have length ncol(Y)");
  }
  if (!(Mv.size() == 1 || Mv.size() == J)) {
    stop("update_cview_cpp: Mv must have length 1 or ncol(Y)");
  }

  List c_view_out = clone(c_view);
  List mu_out = clone(mu);
  List sigma2_out = clone(sigma2);

  for (int j = 0; j < J; ++j) {
    IntegerVector c_view_j = as<IntegerVector>(c_view_out[j]);
    NumericVector mu_j = as<NumericVector>(mu_out[j]);
    NumericVector sigma2_j = as<NumericVector>(sigma2_out[j]);
    IntegerVector gamma_j = as<IntegerVector>(gamma[j]);

    if (c_view_j.size() != m || gamma_j.size() != m) {
      stop("update_cview_cpp: each c_view[[j]] and gamma[[j]] must have length nrow(Y)");
    }

    const double Mj = (Mv.size() == 1) ? Mv[0] : Mv[j];

    for (int i = 0; i < m; ++i) {
      const int old_lab = c_view_j[i];
      const double yij = Y(i, j);
      IntegerVector z_raw(m - 1);
      for (int src = 0, dst = 0; src < m; ++src) if (src != i) z_raw[dst++] = c_view_j[src];
      IntegerVector z_minus_i = canon_cpp(z_raw);
      IntegerVector used = first_occurrence_labels_cpp(z_raw);
      const int K = used.size();
      IntegerVector counts = counts_from_labels(z_minus_i);
      const int fallback_h = reduced_label_for_old_cluster_cpp(used, old_lab, K + 1);

      NumericVector mu_reduced(K), sigma2_reduced(K);
      for (int k = 0; k < K; ++k) {
        const int old_index = used[k] - 1;
        if (old_index < 0 || old_index >= mu_j.size() || old_index >= sigma2_j.size()) {
          stop("update_cview_cpp: cluster atom index out of bounds");
        }
        mu_reduced[k] = mu_j[old_index];
        sigma2_reduced[k] = sigma2_j[old_index];
      }

      std::vector<int> fixed;
      for (int ii = 0; ii < m; ++ii) if (gamma_j[ii] == 1) fixed.push_back(ii);
      NumericVector logw(K + 1, R_NegInf);
      std::vector<IntegerVector> candidates;
      candidates.reserve(K + 1);
      for (int h = 1; h <= K + 1; ++h) {
        IntegerVector candidate = canon_cpp(insert_label_cpp(z_minus_i, i, h));
        candidates.push_back(candidate);
        IntegerVector view_fixed = subset_by_positions(candidate, fixed);
        IntegerVector c0_fixed = subset_by_positions(c0, fixed);
        if (!same_partition_cpp(view_fixed, c0_fixed)) continue;
        if (h <= K) {
          logw[h - 1] = std::log(static_cast<double>(counts[h - 1])) +
                         R::dnorm4(yij, mu_reduced[h - 1],
                                   std::sqrt(sigma2_reduced[h - 1]), 1);
        } else {
          logw[h - 1] = std::log(Mj) +
                         prior_pred_logdens_dispatch(prior_pred_logdens_det_fn,
                                                     yij, mu0[j], tau2[j], a0, b0);
        }
      }

      const int chosen_h = sample_logw_cpp(logw, fallback_h)[0];
      NumericVector mu_candidate = clone(mu_reduced);
      NumericVector sigma2_candidate = clone(sigma2_reduced);
      if (chosen_h == K + 1) {
        NumericVector atom = draw_new_atom_posterior_cpp(yij, mu0[j], tau2[j], a0, b0);
        mu_candidate.push_back(atom[0]);
        sigma2_candidate.push_back(atom[1]);
      }
      // Relabel the candidate and its reduced-coordinate atoms together.  This
      // is essential when insertion occurs before a cluster's first member.
      IntegerVector raw_candidate = insert_label_cpp(z_minus_i, i, chosen_h);
      List relabs = relabel_partition_with_atoms_cpp(raw_candidate,
                                                     mu_candidate,
                                                     sigma2_candidate);
      c_view_j = as<IntegerVector>(relabs["z"]);
      mu_j = as<NumericVector>(relabs["mu"]);
      sigma2_j = as<NumericVector>(relabs["sigma2"]);
    }

    c_view_out[j] = c_view_j;
    mu_out[j] = mu_j;
    sigma2_out[j] = sigma2_j;
  }

  return List::create(Named("c_view") = c_view_out,
                      Named("mu") = mu_out,
                      Named("sigma2") = sigma2_out);
}

// -----------------------------------------------------------------------------
// Joint subject-wise update of (c0_i, c1_i, ..., cJ_i)
// -----------------------------------------------------------------------------

List update_partitions_joint_cpp(
    const NumericMatrix& Y,
    IntegerVector c0,
    List c_view,
    List mu,
    List sigma2,
    const List& gamma,
    const NumericVector& mu0,
    const NumericVector& tau2,
    double a0,
    double b0,
    double M0,
    const NumericVector& Mv,
    const IntegerVector& ord,
    Nullable<Function> prior_pred_logdens_det_fn = R_NilValue) {

  const int m = Y.nrow();
  const int J = Y.ncol();
  if (c0.size() != m || c_view.size() != J || mu.size() != J ||
      sigma2.size() != J || gamma.size() != J) {
    stop("update_partitions_joint_cpp: incompatible state dimensions");
  }
  if (ord.size() == 0) stop("update_partitions_joint_cpp: ord must be non-empty");
  if (!(M0 > 0.0) || !(Mv.size() == 1 || Mv.size() == J)) {
    stop("update_partitions_joint_cpp: invalid CRP parameters");
  }

  // Gamma does not change during this block, so its fixed positions can be
  // cached once for the entire sweep.
  std::vector<std::vector<int> > fixed_by_view(J);
  for (int j = 0; j < J; ++j) {
    IntegerVector gj = as<IntegerVector>(gamma[j]);
    if (gj.size() != m) stop("update_partitions_joint_cpp: invalid gamma length");
    for (int ii = 0; ii < m; ++ii) if (gj[ii] == 1) fixed_by_view[j].push_back(ii);
  }

  for (int oo = 0; oo < ord.size(); ++oo) {
    const int i0 = ord[oo] - 1;
    if (i0 < 0 || i0 >= m) stop("update_partitions_joint_cpp: invalid index in ord");

    // Consensus candidates are represented in the canonical coordinates of
    // the reduced partition, then expanded back to length m.
    const int old_c0_label = c0[i0];
    IntegerVector c0_raw_minus(m - 1);
    for (int src = 0, dst = 0; src < m; ++src) if (src != i0) c0_raw_minus[dst++] = c0[src];
    IntegerVector c0_used = first_occurrence_labels_cpp(c0_raw_minus);
    IntegerVector c0_red = canon_cpp(c0_raw_minus);
    const int K0 = c0_used.size();
    IntegerVector c0_counts = counts_from_labels(c0_red);
    const int c0_fallback = reduced_label_for_old_cluster_cpp(c0_used, old_c0_label, K0 + 1);
    std::vector<IntegerVector> c0_candidates;
    c0_candidates.reserve(K0 + 1);
    for (int h0 = 1; h0 <= K0 + 1; ++h0) {
      c0_candidates.push_back(canon_cpp(insert_label_cpp(c0_red, i0, h0)));
    }

    // Prepare every view once.  Given a consensus candidate, the view
    // candidates are conditionally independent, so their normalizing sums can
    // be multiplied instead of enumerating a Cartesian product.
    std::vector<int> Kj(J), fallback_j(J);
    std::vector<NumericVector> mu_red_j(J), s2_red_j(J), log_base_j(J), log_restricted_j(J);
    std::vector<std::vector<IntegerVector> > raw_candidates_j(J), fixed_candidates_j(J);

    for (int j = 0; j < J; ++j) {
      IntegerVector cj = as<IntegerVector>(c_view[j]);
      NumericVector muj = as<NumericVector>(mu[j]);
      NumericVector s2j = as<NumericVector>(sigma2[j]);
      if (cj.size() != m) stop("update_partitions_joint_cpp: invalid view length");

      const int old_label = cj[i0];
      IntegerVector raw_minus(m - 1);
      for (int src = 0, dst = 0; src < m; ++src) if (src != i0) raw_minus[dst++] = cj[src];
      IntegerVector used = first_occurrence_labels_cpp(raw_minus);
      IntegerVector red = canon_cpp(raw_minus);
      const int K = used.size();
      Kj[j] = K;
      fallback_j[j] = reduced_label_for_old_cluster_cpp(used, old_label, K + 1);
      IntegerVector counts = counts_from_labels(red);

      NumericVector mu_red(K), s2_red(K), log_base(K + 1, R_NegInf);
      for (int k = 0; k < K; ++k) {
        const int old_index = used[k] - 1;
        if (old_index < 0 || old_index >= muj.size() || old_index >= s2j.size()) {
          stop("update_partitions_joint_cpp: atom index out of bounds");
        }
        mu_red[k] = muj[old_index];
        s2_red[k] = s2j[old_index];
        log_base[k] = std::log(static_cast<double>(counts[k])) +
          R::dnorm4(Y(i0, j), mu_red[k], std::sqrt(s2_red[k]), 1);
      }
      const double Mj = (Mv.size() == 1) ? Mv[0] : Mv[j];
      if (!(Mj > 0.0)) stop("update_partitions_joint_cpp: Mv must be positive");
      log_base[K] = std::log(Mj) + prior_pred_logdens_dispatch(
        prior_pred_logdens_det_fn, Y(i0, j), mu0[j], tau2[j], a0, b0);

      mu_red_j[j] = mu_red;
      s2_red_j[j] = s2_red;
      log_base_j[j] = log_base;
      raw_candidates_j[j].reserve(K + 1);
      fixed_candidates_j[j].reserve(K + 1);
      NumericVector log_restricted(K + 1);
      for (int h = 1; h <= K + 1; ++h) {
        IntegerVector raw_candidate = insert_label_cpp(red, i0, h);
        IntegerVector candidate = canon_cpp(raw_candidate);
        IntegerVector fixed_candidate = subset_by_positions(candidate, fixed_by_view[j]);
        raw_candidates_j[j].push_back(raw_candidate);
        fixed_candidates_j[j].push_back(fixed_candidate);
        log_restricted[h - 1] = log_crp_from_sizes_cpp(counts_from_labels(fixed_candidate), Mj);
      }
      log_restricted_j[j] = log_restricted;
    }

    NumericVector logw0(K0 + 1, R_NegInf);
    for (int h0 = 1; h0 <= K0 + 1; ++h0) {
      double score = (h0 <= K0)
        ? std::log(static_cast<double>(c0_counts[h0 - 1]))
        : std::log(M0);
      bool possible = true;
      const IntegerVector& c0_candidate = c0_candidates[h0 - 1];

      for (int j = 0; j < J; ++j) {
        NumericVector logq(Kj[j] + 1, R_NegInf);
        IntegerVector c0_fixed = subset_by_positions(c0_candidate, fixed_by_view[j]);
        for (int h = 1; h <= Kj[j] + 1; ++h) {
          const IntegerVector& view_fixed = fixed_candidates_j[j][h - 1];
          if (!same_partition_cpp(view_fixed, c0_fixed)) continue;
          logq[h - 1] = log_base_j[j][h - 1] - log_restricted_j[j][h - 1];
        }
        const double view_norm = log_sum_exp_core(logq);
        if (!R_finite(view_norm)) {
          possible = false;
          break;
        }
        score += view_norm;
      }
      if (possible) logw0[h0 - 1] = score;
    }

    const int selected_h0 = sample_logw_cpp(logw0, c0_fallback)[0];
    if (selected_h0 < 1 || selected_h0 > K0 + 1 || !R_finite(logw0[selected_h0 - 1])) {
      stop("update_partitions_joint_cpp: no valid consensus candidate");
    }
    c0 = clone(c0_candidates[selected_h0 - 1]);

    // Conditional on the selected consensus candidate, sample each view and
    // carry its atoms through the same canonical relabeling operation.
    for (int j = 0; j < J; ++j) {
      NumericVector logq(Kj[j] + 1, R_NegInf);
      IntegerVector c0_fixed = subset_by_positions(c0, fixed_by_view[j]);
      for (int h = 1; h <= Kj[j] + 1; ++h) {
        const IntegerVector& view_fixed = fixed_candidates_j[j][h - 1];
        if (!same_partition_cpp(view_fixed, c0_fixed)) continue;
        logq[h - 1] = log_base_j[j][h - 1] - log_restricted_j[j][h - 1];
      }
      const int selected_h = sample_logw_cpp(logq, fallback_j[j])[0];
      if (selected_h < 1 || selected_h > Kj[j] + 1 || !R_finite(logq[selected_h - 1])) {
        stop("update_partitions_joint_cpp: no valid view candidate");
      }

      NumericVector mu_candidate = clone(mu_red_j[j]);
      NumericVector s2_candidate = clone(s2_red_j[j]);
      if (selected_h == Kj[j] + 1) {
        NumericVector atom = draw_new_atom_posterior_cpp(
          Y(i0, j), mu0[j], tau2[j], a0, b0);
        mu_candidate.push_back(atom[0]);
        s2_candidate.push_back(atom[1]);
      }
      List relabs = relabel_partition_with_atoms_cpp(
        raw_candidates_j[j][selected_h - 1], mu_candidate, s2_candidate);
      c_view[j] = relabs["z"];
      mu[j] = relabs["mu"];
      sigma2[j] = relabs["sigma2"];
    }
  }

  return List::create(Named("c0") = c0,
                      Named("c_view") = c_view,
                      Named("mu") = mu,
                      Named("sigma2") = sigma2);
}

// -----------------------------------------------------------------------------
// Gamma sweep kernel
// -----------------------------------------------------------------------------

List update_gamma_cpp(List c_view,
                      IntegerVector c0,
                      List gamma,
                      NumericVector alpha,
                      NumericVector Mv,
                      bool independent = true,
                      bool collapse_alpha = false,
                      double a_alpha = 1.0,
                      double b_alpha = 1.0) {
  const int J = c_view.size();
  if (gamma.size() != J || alpha.size() != J) {
    stop("update_gamma_cpp: c_view, gamma, and alpha must have the same length");
  }
  if (!(Mv.size() == 1 || Mv.size() == J)) {
    stop("update_gamma_cpp: Mv must have length 1 or length(c_view)");
  }
  if (collapse_alpha && (!(a_alpha > 0.0) || !(b_alpha > 0.0))) {
    stop("update_gamma_cpp: collapsed-alpha beta parameters must be positive");
  }

  const int m = c0.size();
  List gamma_out = clone(gamma);

  for (int j = 0; j < J; ++j) {
    IntegerVector c_view_j = as<IntegerVector>(c_view[j]);
    IntegerVector gamma_j = as<IntegerVector>(gamma_out[j]);

    if (c_view_j.size() != m || gamma_j.size() != m) {
      stop("update_gamma_cpp: each c_view[[j]] and gamma[[j]] must have length length(c0)");
    }

    const double Mj = (Mv.size() == 1) ? Mv[0] : Mv[j];
    IntegerVector full_counts = counts_from_labels(c_view_j);
    const double log_crp_full = log_crp_from_sizes_cpp(full_counts, Mj);
    IntegerVector fixed_counts(full_counts.size());
    int n_fixed = 0;
    std::vector<int> initial_fixed;
    for (int i = 0; i < m; ++i) {
      if (gamma_j[i] == 1) {
        fixed_counts[c_view_j[i] - 1]++;
        ++n_fixed;
        initial_fixed.push_back(i);
      }
    }
    if (!same_partition_cpp(subset_by_positions(c_view_j, initial_fixed),
                            subset_by_positions(c0, initial_fixed))) {
      stop("update_gamma_cpp: current state violates compatibility");
    }
    double log_crp_fixed = log_crp_from_sizes_cpp(fixed_counts, Mj);

    for (int i = 0; i < m; ++i) {
      const int gi = gamma_j[i];
      NumericVector logw(2, R_NegInf);
      const int lab0 = c_view_j[i] - 1;
      double log_crp_fixed0, log_crp_fixed1;
      if (gi == 0) {
        log_crp_fixed0 = log_crp_fixed;
        const int s = fixed_counts[lab0];
        const double log_add = (s > 0 ? std::log(static_cast<double>(s)) : std::log(Mj)) -
                               std::log(Mj + n_fixed);
        log_crp_fixed1 = log_crp_fixed0 + log_add;
      } else {
        log_crp_fixed1 = log_crp_fixed;
        const int remaining = fixed_counts[lab0] - 1;
        const int n0 = n_fixed - 1;
        const double log_add = (remaining > 0 ? std::log(static_cast<double>(remaining)) : std::log(Mj)) -
                               std::log(Mj + n0);
        log_crp_fixed0 = log_crp_fixed1 - log_add;
      }

      if (collapse_alpha) {
        const int s_minus_i = n_fixed - gi;
        logw[0] = log_crp_full - log_crp_fixed0 +
                  std::log(b_alpha + (m - 1 - s_minus_i));
        logw[1] = log_crp_full - log_crp_fixed1 +
                  std::log(a_alpha + s_minus_i);
      } else {
        logw[0] = log_crp_full - log_crp_fixed0 + std::log(1.0 - alpha[j]);
        logw[1] = log_crp_full - log_crp_fixed1 + std::log(alpha[j]);
      }

      // Given that the currently fixed restriction is compatible, adding i is
      // compatible iff its equality relation agrees with every fixed subject.
      bool compatible_add = true;
      for (int k = 0; k < m; ++k) {
        if (k == i || gamma_j[k] != 1) continue;
        if ((c_view_j[k] == c_view_j[i]) != (c0[k] == c0[i])) {
          compatible_add = false;
          break;
        }
      }
      if (!compatible_add) logw[1] = R_NegInf;

      IntegerVector sampled = sample_logw_cpp(logw, gi + 1);
      int new_gi = sampled[0] - 1;
      if (independent) new_gi = 0;
      if (new_gi != gi) {
        if (new_gi == 1) {
          fixed_counts[lab0]++;
          ++n_fixed;
          log_crp_fixed = log_crp_fixed1;
        } else {
          fixed_counts[lab0]--;
          --n_fixed;
          log_crp_fixed = log_crp_fixed0;
        }
      } else {
        log_crp_fixed = (gi == 1) ? log_crp_fixed1 : log_crp_fixed0;
      }
      gamma_j[i] = new_gi;
    }

    gamma_out[j] = gamma_j;
  }

  return List::create(Named("gamma") = gamma_out);
}


// -----------------------------------------------------------------------------
// Prediction helpers (synchronized predictive particles)
// -----------------------------------------------------------------------------

static inline NumericVector normalize_logw_pred(const NumericVector& logw) {
  const int K = logw.size();
  NumericVector p(K);
  if (K == 0) return p;

  double mx = R_NegInf;
  bool any_finite = false;
  for (int k = 0; k < K; ++k) {
    if (R_finite(logw[k])) {
      any_finite = true;
      if (logw[k] > mx) mx = logw[k];
    }
  }
  if (!any_finite) {
    for (int k = 0; k < K; ++k) p[k] = 1.0 / K;
    return p;
  }

  double s = 0.0;
  for (int k = 0; k < K; ++k) {
    p[k] = R_finite(logw[k]) ? std::exp(logw[k] - mx) : 0.0;
    s += p[k];
  }
  if (!(s > 0.0) || !R_finite(s)) {
    for (int k = 0; k < K; ++k) p[k] = 1.0 / K;
    return p;
  }

  for (int k = 0; k < K; ++k) p[k] /= s;
  return p;
}

static inline double log_mean_exp_pred(const NumericVector& x) {
  if (x.size() == 0) return R_NegInf;
  double mx = R_NegInf;
  bool any_finite = false;
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) {
      any_finite = true;
      if (x[i] > mx) mx = x[i];
    }
  }
  if (!any_finite) return R_NegInf;

  double s = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) {
      s += std::exp(x[i] - mx);
    }
  }
  return mx + std::log(s / static_cast<double>(x.size()));
}

static inline double log_sum_exp_pred(const NumericVector& x) {
  double mx = R_NegInf;
  bool any_finite = false;
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) {
      any_finite = true;
      if (x[i] > mx) mx = x[i];
    }
  }
  if (!any_finite) return R_NegInf;
  double s = 0.0;
  for (int i = 0; i < x.size(); ++i) if (R_finite(x[i])) s += std::exp(x[i] - mx);
  return mx + std::log(s);
}

static inline IntegerVector append_int_pred(const IntegerVector& x, int value) {
  IntegerVector out(x.size() + 1);
  for (int i = 0; i < x.size(); ++i) out[i] = x[i];
  out[x.size()] = value;
  return out;
}


static inline int max_label_pred(const IntegerVector& z) {
  int m = 0;
  for (int i = 0; i < z.size(); ++i) if (z[i] > m) m = z[i];
  return m;
}

static inline IntegerVector counts_from_labels_pred(const IntegerVector& z, int K) {
  IntegerVector counts(K);
  for (int i = 0; i < z.size(); ++i) {
    const int lab = z[i];
    if (lab >= 1 && lab <= K) counts[lab - 1]++;
  }
  return counts;
}

static inline std::vector<int> fixed_positions_from_gamma_pred(const IntegerVector& gaug_pred) {
  std::vector<int> fixed;
  fixed.reserve(gaug_pred.size());
  for (int i = 0; i < gaug_pred.size(); ++i) {
    if (gaug_pred[i] == 1L) fixed.push_back(i);
  }
  return fixed;
}

List update_g_pred_cpp(
    const List& c_view,
    const IntegerVector& c0,
    const List& gamma,
    const IntegerMatrix& g_pred,
    const IntegerVector& c0_pred,
    const IntegerMatrix& cview_pred,
    const NumericVector& alpha,
    const NumericVector& Mv,
    bool independent) {
  const int J = c_view.size();
  const int n_particles = g_pred.nrow();
  IntegerMatrix out = clone(g_pred);

  for (int j = 0; j < J; ++j) {
    const IntegerVector c_view_j = as<IntegerVector>(c_view[j]);
    const IntegerVector gamma_j = as<IntegerVector>(gamma[j]);

    for (int rr = 0; rr < n_particles; ++rr) {
      NumericVector logw(2, R_NegInf);

      IntegerVector g0_pred = append_int_pred(gamma_j, 0);
      IntegerVector g1_pred = append_int_pred(gamma_j, 1);
      IntegerVector cjaug_pred = append_int_pred(c_view_j, cview_pred(rr, j));
      IntegerVector c0aug_pred = append_int_pred(c0, c0_pred[rr]);

      const std::vector<int> fixed0 = fixed_positions_from_gamma_pred(g0_pred);
      const std::vector<int> fixed1 = fixed_positions_from_gamma_pred(g1_pred);

      IntegerVector cj_fixed0(fixed0.size()), c0_fixed0(fixed0.size());
      IntegerVector cj_fixed1(fixed1.size()), c0_fixed1(fixed1.size());
      for (int k = 0; k < (int)fixed0.size(); ++k) {
        cj_fixed0[k] = cjaug_pred[fixed0[k]];
        c0_fixed0[k] = c0aug_pred[fixed0[k]];
      }
      for (int k = 0; k < (int)fixed1.size(); ++k) {
        cj_fixed1[k] = cjaug_pred[fixed1[k]];
        c0_fixed1[k] = c0aug_pred[fixed1[k]];
      }

      const double Mj = (Mv.size() == 1) ? Mv[0] : Mv[j];
      const IntegerVector cj_sizes = counts_from_labels_pred(cjaug_pred, max_label_pred(cjaug_pred));
      logw[0] = log_crp_from_sizes_cpp(cj_sizes, Mj) -
                log_crp_from_sizes_cpp(counts_from_labels_pred(cj_fixed0, max_label_pred(cj_fixed0)), Mj) +
                std::log(1.0 - alpha[j]);

      if (same_partition_cpp(cj_fixed1, c0_fixed1)) {
        logw[1] = log_crp_from_sizes_cpp(cj_sizes, Mj) -
                  log_crp_from_sizes_cpp(counts_from_labels_pred(cj_fixed1, max_label_pred(cj_fixed1)), Mj) +
                  std::log(alpha[j]);
      }

      if (independent) {
        out(rr, j) = 0;
      } else {
        IntegerVector sampled = sample_logw_cpp(logw, out(rr, j) + 1);
        out(rr, j) = sampled[0] - 1;
      }
    }
  }

  return List::create(Named("g_pred") = out);
}

IntegerVector update_c0_pred_cpp(
    const IntegerVector& c0,
    const List& c_view,
    const List& gamma,
    const IntegerVector& c0_pred,
    const IntegerMatrix& g_pred,
    const IntegerMatrix& cview_pred,
    const NumericVector& Mv,
    double M0) {
  const int J = c_view.size();
  const int n_particles = g_pred.nrow();
  IntegerVector out = clone(c0_pred);

  for (int rr = 0; rr < n_particles; ++rr) {
    const int nclus = max_label_pred(c0);
    NumericVector logw(nclus + 1, R_NegInf);

    for (int h = 1; h <= nclus + 1; ++h) {
      double lw = 0.0;
      bool compatible = true;
      for (int j = 0; j < J; ++j) {
        const IntegerVector gamma_j = as<IntegerVector>(gamma[j]);
        IntegerVector gaug_pred = append_int_pred(gamma_j, g_pred(rr, j));
        IntegerVector c0aug_pred = append_int_pred(c0, h);
        IntegerVector cjaug_pred = append_int_pred(as<IntegerVector>(c_view[j]), cview_pred(rr, j));

        const std::vector<int> fixed = fixed_positions_from_gamma_pred(gaug_pred);
        IntegerVector cj_fixed(fixed.size()), c0_fixed(fixed.size());
        for (int k = 0; k < (int)fixed.size(); ++k) {
          cj_fixed[k] = cjaug_pred[fixed[k]];
          c0_fixed[k] = c0aug_pred[fixed[k]];
        }

        if (!same_partition_cpp(cj_fixed, c0_fixed)) {
          compatible = false;
          break;
        }
        IntegerVector cj_sizes = counts_from_labels_pred(cjaug_pred, max_label_pred(cjaug_pred));
        IntegerVector cj_fixed_sizes = counts_from_labels_pred(cj_fixed, max_label_pred(cj_fixed));
        lw += log_crp_from_sizes_cpp(cj_sizes, (Mv.size() == 1) ? Mv[0] : Mv[j]) -
              log_crp_from_sizes_cpp(cj_fixed_sizes, (Mv.size() == 1) ? Mv[0] : Mv[j]);
      }
      if (!compatible) continue;
      if (h <= nclus) {
        int count_h = 0;
        for (int i = 0; i < c0.size(); ++i) if (c0[i] == h) ++count_h;
        lw += std::log(static_cast<double>(count_h));
      } else {
        lw += std::log(M0);
      }
      logw[h - 1] = lw;
    }

    IntegerVector sampled = sample_logw_cpp(logw, out[rr]);
    out[rr] = sampled[0];
  }

  return out;
}

List update_cview_pred_cpp(
    const NumericMatrix& Y,
    const List& c_view,
    const List& mu,
    const List& sigma2,
    const IntegerVector& c0,
    const List& gamma,
    const IntegerVector& c0_pred,
    const IntegerMatrix& g_pred,
    const IntegerMatrix& cview_pred,
    const NumericVector& mu0,
    const NumericVector& tau2,
    double a0,
    double b0,
    const NumericVector& Mv,
    Nullable<Function> prior_pred_logdens_det_fn = R_NilValue) {
  const int J = Y.ncol();
  const int n_particles = g_pred.nrow();
  IntegerMatrix cview_out = clone(cview_pred);
  List cprob_pred(J);

  for (int j = 0; j < J; ++j) {
    const IntegerVector c_view_j = as<IntegerVector>(c_view[j]);
    const IntegerVector gamma_j = as<IntegerVector>(gamma[j]);
    const int nclus = max_label_pred(c_view_j);
    List cprob_j(n_particles);

    for (int rr = 0; rr < n_particles; ++rr) {
      IntegerVector gaug_pred = append_int_pred(gamma_j, g_pred(rr, j));
      IntegerVector c0aug_pred = append_int_pred(c0, c0_pred[rr]);
      const std::vector<int> fixed = fixed_positions_from_gamma_pred(gaug_pred);
      const double Mj = (Mv.size() == 1) ? Mv[0] : Mv[j];

      NumericVector logw(nclus + 1, R_NegInf);
      for (int h = 1; h <= nclus + 1; ++h) {
        IntegerVector cjaug_pred = append_int_pred(c_view_j, h);
        IntegerVector cj_fixed(fixed.size()), c0_fixed(fixed.size());
        for (int k = 0; k < (int)fixed.size(); ++k) {
          cj_fixed[k] = cjaug_pred[fixed[k]];
          c0_fixed[k] = c0aug_pred[fixed[k]];
        }
        if (same_partition_cpp(cj_fixed, c0_fixed)) {
          // Exact predictive restricted-CRP allocation weight:
          //   Pr(rho_j^+) / Pr((rho_j^+)^{gamma_j^+}) * compatibility.
          // The restricted denominator is constant over compatible h once
          // c0_pred and g_pred are fixed, but retaining it here makes the
          // predictive calculation match the model density explicitly.
          const IntegerVector cj_sizes = counts_from_labels_pred(
            cjaug_pred, max_label_pred(cjaug_pred));
          const IntegerVector cj_fixed_sizes = counts_from_labels_pred(
            cj_fixed, max_label_pred(cj_fixed));
          logw[h - 1] = log_crp_from_sizes_cpp(cj_sizes, Mj) -
                        log_crp_from_sizes_cpp(cj_fixed_sizes, Mj);
        }
      }

      NumericVector p = normalize_logw_pred(logw);
      cprob_j[rr] = p;
      IntegerVector sampled = sample_logw_cpp(logw, cview_out(rr, j));
      cview_out(rr, j) = sampled[0];
    }

    cprob_pred[j] = cprob_j;
  }

  return List::create(Named("cview_pred") = cview_out,
                      Named("cprob_pred") = cprob_pred);
}

// Exact joint predictive density for one completely held-out J-variate
// observation. Conditional on a consensus candidate, the view sums factor.
static double exact_joint_predict_logdens_cpp(
    const NumericVector& y_test,
    const IntegerVector& c0,
    const List& c_view,
    const List& gamma,
    const NumericVector& alpha,
    const List& mu,
    const List& sigma2,
    const NumericVector& mu0,
    const NumericVector& tau2,
    double M0,
    const NumericVector& Mv,
    double a_alpha,
    double b_alpha,
    double a0,
    double b0,
    bool independent,
    bool marginalize_alpha_prediction,
    Nullable<Function> prior_pred_logdens_det_fn) {
  const int m = c0.size();
  const int J = c_view.size();
  if (y_test.size() != J || gamma.size() != J || alpha.size() != J ||
      mu.size() != J || sigma2.size() != J || mu0.size() != J ||
      tau2.size() != J || !(Mv.size() == 1 || Mv.size() == J))
    stop("exact_joint_predict_logdens_cpp: incompatible dimensions");
  if (!(M0 > 0.0) || !(a_alpha > 0.0) || !(b_alpha > 0.0) ||
      !(a0 > 0.0) || !(b0 > 0.0))
    stop("exact_joint_predict_logdens_cpp: hyperparameters must be positive");

  const int K0 = max_label_pred(c0);
  IntegerVector n0 = counts_from_labels_pred(c0, K0);
  NumericVector log_den_h0(K0 + 1, R_NegInf), log_num_h0(K0 + 1, R_NegInf);

  for (int h0 = 1; h0 <= K0 + 1; ++h0) {
    const double log_a0 = h0 <= K0 ? std::log(static_cast<double>(n0[h0 - 1])) : std::log(M0);
    IntegerVector c0_aug = append_int_pred(c0, h0);
    double log_den = log_a0, log_num = log_a0;
    bool possible = true;

    for (int j = 0; j < J; ++j) {
      const IntegerVector cj = as<IntegerVector>(c_view[j]);
      const IntegerVector gj = as<IntegerVector>(gamma[j]);
      const NumericVector muj = as<NumericVector>(mu[j]);
      const NumericVector s2j = as<NumericVector>(sigma2[j]);
      if (cj.size() != m || gj.size() != m || muj.size() != s2j.size())
        stop("exact_joint_predict_logdens_cpp: invalid view state");
      const int Kj = max_label_pred(cj);
      if (muj.size() != Kj)
        stop("exact_joint_predict_logdens_cpp: atoms do not match occupied clusters");
      const double Mj = Mv.size() == 1 ? Mv[0] : Mv[j];
      if (!(Mj > 0.0)) stop("exact_joint_predict_logdens_cpp: Mv must be positive");

      int sj = 0;
      for (int i = 0; i < m; ++i) sj += gj[i];
      double q1 = marginalize_alpha_prediction
        ? (a_alpha + sj) / (a_alpha + b_alpha + m) : alpha[j];
      if (independent) q1 = 0.0;
      if (!(q1 >= 0.0 && q1 <= 1.0))
        stop("exact_joint_predict_logdens_cpp: invalid predictive gamma probability");

      NumericVector ltd(2 * (Kj + 1), R_NegInf), ltn(2 * (Kj + 1), R_NegInf);
      int term = 0;
      for (int g = 0; g <= 1; ++g) {
        if (independent && g == 1) { term += Kj + 1; continue; }
        const double pg = g == 1 ? q1 : 1.0 - q1;
        if (!(pg > 0.0)) { term += Kj + 1; continue; }
        IntegerVector g_aug = append_int_pred(gj, g);
        const std::vector<int> fixed = fixed_positions_from_gamma_pred(g_aug);
        for (int h = 1; h <= Kj + 1; ++h, ++term) {
          IntegerVector cj_aug = append_int_pred(cj, h);
          IntegerVector cj_fixed(fixed.size()), c0_fixed(fixed.size());
          for (int k = 0; k < static_cast<int>(fixed.size()); ++k) {
            cj_fixed[k] = cj_aug[fixed[k]];
            c0_fixed[k] = c0_aug[fixed[k]];
          }
          if (!same_partition_cpp(cj_fixed, c0_fixed)) continue;
          const double log_rpm =
            log_crp_from_sizes_cpp(counts_from_labels_pred(cj_aug, max_label_pred(cj_aug)), Mj) -
            log_crp_from_sizes_cpp(counts_from_labels_pred(cj_fixed, max_label_pred(cj_fixed)), Mj);
          const double log_latent = std::log(pg) + log_rpm;
          const double log_lik = h <= Kj
            ? R::dnorm4(y_test[j], muj[h - 1], std::sqrt(s2j[h - 1]), 1)
            : prior_pred_logdens_dispatch(prior_pred_logdens_det_fn,
                y_test[j], mu0[j], tau2[j], a0, b0);
          ltd[term] = log_latent;
          ltn[term] = log_latent + log_lik;
        }
      }
      const double vd = log_sum_exp_pred(ltd), vn = log_sum_exp_pred(ltn);
      if (!R_finite(vd) || !R_finite(vn)) { possible = false; break; }
      log_den += vd;
      log_num += vn;
    }
    if (possible) {
      log_den_h0[h0 - 1] = log_den;
      log_num_h0[h0 - 1] = log_num;
    }
  }
  const double log_den = log_sum_exp_pred(log_den_h0);
  const double log_num = log_sum_exp_pred(log_num_h0);
  if (!R_finite(log_den) || !R_finite(log_num))
    stop("exact_joint_predict_logdens_cpp: no compatible predictive state");
  return log_num - log_den;
}

List generate_prediction_outputs_cpp(
    const NumericMatrix& Y,
    const NumericVector& y_test,
    const List& c_view,
    const List& mu,
    const List& sigma2,
    const NumericVector& mu0,
    const NumericVector& tau2,
    const IntegerMatrix& g_pred,
    const IntegerVector& c0_pred,
    const IntegerMatrix& cview_pred,
    const List& cprob_pred,
    int n_new,
    double a0,
    double b0,
    Nullable<Function> prior_pred_logdens_det_fn = R_NilValue) {
  const int J = Y.ncol();
  const int n_particles = g_pred.nrow();
  NumericMatrix y_new(n_new, J);
  NumericMatrix y_new_rb(n_new, J);
  NumericVector lp_rr(n_particles);
  if (!(y_test.size() == 0 || y_test.size() == J)) {
    stop("generate_prediction_outputs_cpp: y_test must be empty or have length ncol(Y)");
  }
  const bool compute_lp = (y_test.size() == J);
  double lp_test = NA_REAL;

  if (compute_lp) {
    for (int rr = 0; rr < n_particles; ++rr) {
      double lp_this_rr = 0.0;
      for (int j = 0; j < J; ++j) {
        List cprob_j = as<List>(cprob_pred[j]);
        NumericVector cprob = as<NumericVector>(cprob_j[rr]);
        NumericVector mu_j = as<NumericVector>(mu[j]);
        NumericVector s2_j = as<NumericVector>(sigma2[j]);
        NumericVector log_comp(mu_j.size() + 1);
        for (int k = 0; k < mu_j.size(); ++k) {
          log_comp[k] = R::dnorm4(y_test[j], mu_j[k], std::sqrt(s2_j[k]), 1);
        }
        log_comp[mu_j.size()] = prior_pred_logdens_dispatch(
          prior_pred_logdens_det_fn, y_test[j], mu0[j], tau2[j], a0, b0);
        lp_this_rr += log_sum_exp_pred(log(cprob) + log_comp);
      }
      lp_rr[rr] = lp_this_rr;
    }
    lp_test = log_mean_exp_pred(lp_rr);
  }

  // One shared predictive particle per new subject, as in the R algorithm.
  for (int i = 0; i < n_new; ++i) {
    int rr_star = static_cast<int>(std::floor(R::runif(0.0, n_particles)));
    if (rr_star >= n_particles) rr_star = n_particles - 1;
    for (int j = 0; j < J; ++j) {
      List cprob_j = as<List>(cprob_pred[j]);
      NumericVector cprob = as<NumericVector>(cprob_j[rr_star]);
      NumericVector mu_j = as<NumericVector>(mu[j]);
      NumericVector s2_j = as<NumericVector>(sigma2[j]);
      NumericVector mu_tmp(mu_j.size() + 1);
      for (int k = 0; k < mu_j.size(); ++k) mu_tmp[k] = mu_j[k];
      mu_tmp[mu_j.size()] = mu0[j];
      double rb = 0.0;
      for (int k = 0; k < cprob.size(); ++k) rb += cprob[k] * mu_tmp[k];
      y_new_rb(i, j) = rb;

      const int chosen = sample_logw_cpp(log(cprob), 1)[0];
      if (chosen <= mu_j.size()) {
        y_new(i, j) = R::rnorm(mu_j[chosen - 1], std::sqrt(s2_j[chosen - 1]));
      } else {
        const double mu_star = R::rnorm(mu0[j], std::sqrt(tau2[j]));
        const double s2_star = 1.0 / R::rgamma(a0, 1.0 / b0);
        y_new(i, j) = R::rnorm(mu_star, std::sqrt(s2_star));
      }
    }
  }

  return List::create(Named("y_new") = y_new,
                      Named("y_new_rb") = y_new_rb,
                      Named("lp_test") = lp_test);
}

// -----------------------------------------------------------------------------
// Outer Gibbs driver (initial scaffold)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List gibbs_mv_rpm_cpp(
    NumericMatrix Y,
    int n_iter = 1000,
    int burn = 0,
    int thin = 1,
    double M0 = 1.0,
    NumericVector Mv = NumericVector::create(1.0),
    double a_alpha = 1.0,
    double b_alpha = 1.0,
    double m0 = 0.0,
    double s0_sq = 10.0,
    double a0 = 2.0,
    double b0 = 1.0,
    double atau = 3.0,
    double btau = 1.0,
    int n_new = 1,
    NumericVector y_test = NumericVector::create(),
    bool independent = false,
    Nullable<IntegerVector> c0_init = R_NilValue,
    Nullable<List> c_view_init = R_NilValue,
    Nullable<List> gamma_init = R_NilValue,
    Nullable<NumericVector> alpha_init = R_NilValue,
    Nullable<Function> prior_pred_logdens_det_fn = R_NilValue,
    bool do_prediction = true,
    int n_particles = 50,
    int n_pred_sweeps = 1,
    bool exact_prediction = true,
    bool marginalize_alpha_prediction = true) {

  const int m = Y.nrow();
  const int J = Y.ncol();
  if (m < 1 || J < 1) stop("gibbs_mv_rpm_cpp: Y must be non-empty");
  if (thin < 1 || burn < 0 || n_iter < burn) stop("gibbs_mv_rpm_cpp: invalid MCMC controls");
  if (n_new < 1) stop("gibbs_mv_rpm_cpp: n_new must be positive");
  if (do_prediction && n_particles < 1) stop("gibbs_mv_rpm_cpp: n_particles must be positive");
  if (do_prediction && n_pred_sweeps < 1) stop("gibbs_mv_rpm_cpp: n_pred_sweeps must be positive");
  if (!(y_test.size() == 0 || y_test.size() == J)) {
    stop("gibbs_mv_rpm_cpp: y_test must be omitted or have length ncol(Y)");
  }
  if (Mv.size() == 1) Mv = NumericVector(J, Mv[0]);
  if (Mv.size() != J) stop("gibbs_mv_rpm_cpp: Mv must have length 1 or ncol(Y)");
  if (!(M0 > 0.0)) stop("gibbs_mv_rpm_cpp: M0 must be positive");
  if (!(a_alpha > 0.0) || !(b_alpha > 0.0)) {
    stop("gibbs_mv_rpm_cpp: alpha beta-prior parameters must be positive");
  }
  for (int j = 0; j < J; ++j) if (!(Mv[j] > 0.0)) stop("gibbs_mv_rpm_cpp: Mv must be positive");

  Function sample_fn("sample");
  Function kmeans_fn("kmeans");

  IntegerVector c0;
  if (c0_init.isNotNull()) {
    c0 = canon_cpp(as<IntegerVector>(c0_init));
  } else {
    IntegerVector pool(2); pool[0] = 1; pool[1] = 2;
    IntegerVector tmp = as<IntegerVector>(sample_fn(pool, Named("size") = m, Named("replace") = true));
    c0 = canon_cpp(tmp);
  }

  List c_view(J);
  if (c_view_init.isNotNull()) {
    List init = as<List>(c_view_init);
    if (init.size() != J) stop("gibbs_mv_rpm_cpp: c_view_init has wrong length");
    for (int j = 0; j < J; ++j) c_view[j] = canon_cpp(as<IntegerVector>(init[j]));
  } else {
    for (int j = 0; j < J; ++j) {
      NumericVector col = Y(_, j);
      List km = kmeans_fn(col, Named("centers") = 4);
      IntegerVector cl = as<IntegerVector>(km["cluster"]);
      c_view[j] = canon_cpp(cl);
    }
  }

  List gamma(J);
  if (gamma_init.isNotNull()) {
    List init = as<List>(gamma_init);
    for (int j = 0; j < J; ++j) gamma[j] = as<IntegerVector>(init[j]);
  } else {
    for (int j = 0; j < J; ++j) gamma[j] = IntegerVector(m, 0);
  }
  if (c0.size() != m) stop("gibbs_mv_rpm_cpp: c0_init has wrong length");
  for (int j = 0; j < J; ++j) {
    IntegerVector cj = as<IntegerVector>(c_view[j]);
    IntegerVector gj = as<IntegerVector>(gamma[j]);
    if (cj.size() != m || gj.size() != m) stop("gibbs_mv_rpm_cpp: invalid initial state length");
    std::vector<int> fixed;
    for (int i = 0; i < m; ++i) {
      if (gj[i] != 0 && gj[i] != 1) stop("gibbs_mv_rpm_cpp: gamma_init must be binary");
      if (gj[i] == 1) fixed.push_back(i);
    }
    if (!same_partition_cpp(subset_by_positions(cj, fixed), subset_by_positions(c0, fixed))) {
      stop("gibbs_mv_rpm_cpp: initial view and consensus partitions are incompatible with gamma_init");
    }
  }

  NumericVector alpha(J);
  if (alpha_init.isNotNull()) {
    alpha = as<NumericVector>(alpha_init);
    if (alpha.size() != J) stop("gibbs_mv_rpm_cpp: alpha_init has wrong length");
  } else {
    for (int j = 0; j < J; ++j) alpha[j] = 0.5;
  }

  List mu(J), sigma2(J);
  for (int j = 0; j < J; ++j) {
    IntegerVector cj = as<IntegerVector>(c_view[j]);
    int nclus = 0;
    for (int i = 0; i < cj.size(); ++i) if (cj[i] > nclus) nclus = cj[i];
    NumericVector mu_j(nclus);
    NumericVector s2_j(nclus, 0.5);
    for (int h = 1; h <= nclus; ++h) {
      double sumy = 0.0;
      int nh = 0;
      for (int i = 0; i < m; ++i) if (cj[i] == h) {
        sumy += Y(i, j);
        ++nh;
      }
      mu_j[h - 1] = (nh > 0) ? (sumy / nh) : 0.0;
      s2_j[h - 1] = 0.5;
    }
    mu[j] = mu_j;
    sigma2[j] = s2_j;
  }

  NumericVector mu0(J), tau2(J);
  for (int j = 0; j < J; ++j) { mu0[j] = 0.0; tau2[j] = 1.0; }

  const int particle_count = do_prediction ? n_particles : 0;
  IntegerMatrix g_pred(particle_count, J);
  IntegerVector c0_pred(particle_count);
  IntegerMatrix cview_pred(particle_count, J);
  List cprob_pred(J);
  for (int j = 0; j < J; ++j) cprob_pred[j] = List(particle_count);
  for (int rr = 0; rr < particle_count; ++rr) {
    c0_pred[rr] = 1;
    for (int j = 0; j < J; ++j) {
      g_pred(rr, j) = 0;
      cview_pred(rr, j) = 1;
    }
  }

  const int n_keep = std::max(0, (n_iter - burn) / thin);
  IntegerMatrix c0_samp(n_keep, m);
  IntegerVector c_view_samp(n_keep * m * J);
  c_view_samp.attr("dim") = IntegerVector::create(n_keep, m, J);
  IntegerVector gamma_samp(n_keep * m * J);
  gamma_samp.attr("dim") = IntegerVector::create(n_keep, m, J);
  NumericMatrix alpha_samp(n_keep, J);
  NumericMatrix alpha_rb_samp(n_keep, J);
  List mu_samp(n_keep), sigma2_samp(n_keep);
  NumericMatrix mu0_samp(n_keep, J), tau2_samp(n_keep, J);
  NumericMatrix y_new_samp(n_keep, J), y_new_rb_samp(n_keep, J);
  NumericVector y_in_samp(n_keep * m * J);
  y_in_samp.attr("dim") = IntegerVector::create(n_keep, m, J);
  NumericVector lp_test_samp(n_keep);
  NumericVector lp_test_exact_samp(n_keep);
  NumericVector gamma_sum_samp(n_keep);
  List prediction_state_samp(n_keep);
  std::fill(y_new_samp.begin(), y_new_samp.end(), NA_REAL);
  std::fill(y_new_rb_samp.begin(), y_new_rb_samp.end(), NA_REAL);
  std::fill(y_in_samp.begin(), y_in_samp.end(), NA_REAL);
  std::fill(lp_test_samp.begin(), lp_test_samp.end(), NA_REAL);
  std::fill(lp_test_exact_samp.begin(), lp_test_exact_samp.end(), NA_REAL);

  int keep_ctr = 0;

  for (int iter = 1; iter <= n_iter; ++iter) {
    if (iter % 500 == 0) {
      Rcout << "iter = " << iter << "\n";
    }

    // Partially collapsed update: integrate alpha out only while sampling
    // gamma, then immediately restore alpha with a draw from alpha | gamma.
    gamma = as<List>(update_gamma_cpp(
      c_view, c0, gamma, alpha, Mv, independent,
      true, a_alpha, b_alpha)["gamma"]);
    for (int j = 0; j < J; ++j) {
      int gsum = 0;
      IntegerVector gj = as<IntegerVector>(gamma[j]);
      for (int i = 0; i < m; ++i) gsum += gj[i];
      alpha[j] = independent
        ? 0.0
        : R::rbeta(a_alpha + gsum, b_alpha + m - gsum);
    }

    // Jointly update each subject's consensus and all view labels.  This block
    // can cross compatibility barriers that trap separate c0/cj sweeps.
    IntegerVector pool(m); for (int i = 0; i < m; ++i) pool[i] = i + 1;
    IntegerVector ord = as<IntegerVector>(sample_fn(pool, Named("size") = m, Named("replace") = false));
    List joint_res = update_partitions_joint_cpp(
      Y, c0, c_view, mu, sigma2, gamma, mu0, tau2,
      a0, b0, M0, Mv, ord, prior_pred_logdens_det_fn);
    c0 = as<IntegerVector>(joint_res["c0"]);
    c_view = as<List>(joint_res["c_view"]);
    mu = as<List>(joint_res["mu"]);
    sigma2 = as<List>(joint_res["sigma2"]);

    // --------------------------------
    // 4) Update mu and sigma2
    // --------------------------------
    for (int j = 0; j < J; ++j) {
      IntegerVector cj = as<IntegerVector>(c_view[j]);
      NumericVector mu_j = as<NumericVector>(mu[j]);
      NumericVector s2_j = as<NumericVector>(sigma2[j]);

      int nclus = 0;
      for (int i = 0; i < cj.size(); ++i) if (cj[i] > nclus) nclus = cj[i];
      IntegerVector cluster_n(nclus);
      NumericVector cluster_sum(nclus), cluster_sumsq(nclus);
      for (int i = 0; i < m; ++i) {
        const int h0 = cj[i] - 1;
        const double y = Y(i, j);
        cluster_n[h0]++;
        cluster_sum[h0] += y;
        cluster_sumsq[h0] += y * y;
      }

      for (int h = 1; h <= nclus; ++h) {
        const int idx = h - 1;
        const int n_h = cluster_n[idx];
        const double sumy = cluster_sum[idx];
        const double s2_curr = s2_j[h - 1];
        const double v_mu = 1.0 / (1.0 / tau2[j] + n_h / s2_curr);
        const double m_mu = v_mu * (mu0[j] / tau2[j] + sumy / s2_curr);
        mu_j[h - 1] = R::rnorm(m_mu, std::sqrt(v_mu));

        const double a_post = a0 + n_h / 2.0;
        double ss = cluster_sumsq[idx] - 2.0 * mu_j[idx] * sumy +
                    n_h * mu_j[idx] * mu_j[idx];
        if (ss < 0.0 && ss > -1e-10) ss = 0.0;
        const double b_post = b0 + 0.5 * ss;
        s2_j[h - 1] = 1.0 / R::rgamma(a_post, 1.0 / b_post);
      }

      const double v_mu0 = 1.0 / (1.0 / s0_sq + nclus / tau2[j]);
      const double m_mu0 = v_mu0 * (m0 / s0_sq + sum(mu_j) / tau2[j]);
      mu0[j] = R::rnorm(m_mu0, std::sqrt(v_mu0));

      const double a_post_tau = atau + nclus / 2.0;
      const double b_post_tau = btau + 0.5 * sum((mu_j - mu0[j]) * (mu_j - mu0[j]));
      tau2[j] = 1.0 / R::rgamma(a_post_tau, 1.0 / b_post_tau);

      mu[j] = mu_j;
      sigma2[j] = s2_j;
    }

    if (do_prediction) {
      for (int ps = 0; ps < n_pred_sweeps; ++ps) {
        g_pred = as<IntegerMatrix>(update_g_pred_cpp(
          c_view, c0, gamma, g_pred, c0_pred, cview_pred,
          alpha, Mv, independent)["g_pred"]);
        c0_pred = update_c0_pred_cpp(c0, c_view, gamma, c0_pred,
                                     g_pred, cview_pred, Mv, M0);
        List cvpred_res = update_cview_pred_cpp(
          Y, c_view, mu, sigma2, c0, gamma, c0_pred, g_pred,
          cview_pred, mu0, tau2, a0, b0, Mv,
          prior_pred_logdens_det_fn);
        cview_pred = as<IntegerMatrix>(cvpred_res["cview_pred"]);
        cprob_pred = as<List>(cvpred_res["cprob_pred"]);
      }
    }

    if (iter > burn && ((iter - burn) % thin == 0)) {
      ++keep_ctr;
      NumericMatrix y_insample(m, J);
      NumericMatrix y_new(n_new, J), y_new_rb(n_new, J);
      double lp_test = NA_REAL;
      double lp_test_exact = NA_REAL;
      if (do_prediction) {
        for (int i = 0; i < m; ++i) {
          for (int j = 0; j < J; ++j) {
            IntegerVector cj = as<IntegerVector>(c_view[j]);
            NumericVector mu_j = as<NumericVector>(mu[j]);
            NumericVector s2_j = as<NumericVector>(sigma2[j]);
            const int h0 = cj[i] - 1;
            y_insample(i, j) = R::rnorm(mu_j[h0], std::sqrt(s2_j[h0]));
          }
        }
        List pred_res = generate_prediction_outputs_cpp(
          Y, y_test, c_view, mu, sigma2, mu0, tau2, g_pred,
          c0_pred, cview_pred, cprob_pred, n_new, a0, b0,
          prior_pred_logdens_det_fn);
        y_new = as<NumericMatrix>(pred_res["y_new"]);
        y_new_rb = as<NumericMatrix>(pred_res["y_new_rb"]);
        lp_test = as<double>(pred_res["lp_test"]);
        if (exact_prediction && y_test.size() == J) {
          lp_test_exact = exact_joint_predict_logdens_cpp(
            y_test, c0, c_view, gamma, alpha, mu, sigma2, mu0, tau2,
            M0, Mv, a_alpha, b_alpha, a0, b0, independent,
            marginalize_alpha_prediction, prior_pred_logdens_det_fn);
        }
      }
      for (int i = 0; i < m; ++i) c0_samp(keep_ctr - 1, i) = c0[i];
      for (int j = 0; j < J; ++j) {
        alpha_samp(keep_ctr - 1, j) = alpha[j];
        IntegerVector gj_for_alpha = as<IntegerVector>(gamma[j]);
        int gsum_j = 0;
        for (int i = 0; i < m; ++i) gsum_j += gj_for_alpha[i];
        alpha_rb_samp(keep_ctr - 1, j) =
          (a_alpha + gsum_j) / (a_alpha + b_alpha + m);
        mu0_samp(keep_ctr - 1, j) = mu0[j];
        tau2_samp(keep_ctr - 1, j) = tau2[j];
        if (do_prediction) {
          y_new_samp(keep_ctr - 1, j) = y_new(0, j);
          y_new_rb_samp(keep_ctr - 1, j) = y_new_rb(0, j);
        }
      }
      c_view_samp.attr("dim");
      gamma_samp.attr("dim");
      // Store c_view and gamma as arrays with dimensions (n_keep, m, J).
      // In R's column-major layout, the linear index is:
      //   iter + n_keep * i + n_keep * m * j
      for (int j = 0; j < J; ++j) {
        IntegerVector cj = as<IntegerVector>(c_view[j]);
        IntegerVector gj = as<IntegerVector>(gamma[j]);
        for (int i = 0; i < m; ++i) {
          const int idx = keep_ctr - 1 + n_keep * i + n_keep * m * j;
          c_view_samp[idx] = cj[i];
          gamma_samp[idx] = gj[i];
          if (do_prediction) y_in_samp[idx] = y_insample(i, j);
        }
      }
      // Deep copies are required: later atom updates mutate list elements.
      mu_samp[keep_ctr - 1] = clone(mu);
      sigma2_samp[keep_ctr - 1] = clone(sigma2);
      if (do_prediction) lp_test_samp[keep_ctr - 1] = lp_test;
      if (do_prediction && exact_prediction) lp_test_exact_samp[keep_ctr - 1] = lp_test_exact;

      int gsum_all = 0;
      for (int jj = 0; jj < J; ++jj) {
        IntegerVector gj = as<IntegerVector>(gamma[jj]);
        for (int ii = 0; ii < m; ++ii) gsum_all += gj[ii];
      }
      gamma_sum_samp[keep_ctr - 1] = gsum_all;

      if (do_prediction) {
        prediction_state_samp[keep_ctr - 1] = List::create(
          Named("g_pred") = clone(g_pred),
          Named("c0_pred") = clone(c0_pred),
          Named("cview_pred") = clone(cview_pred),
          Named("cprob_pred") = clone(cprob_pred)
        );
      }
    }
  }

  // Each lp_test element is conditional on one retained posterior draw and
  // has already integrated over the synchronized predictive particles.  LOO
  // posterior prediction must average these densities, not their logarithms.
  double lp_test_log_mean_exp = NA_REAL;
  if (do_prediction && y_test.size() == J && keep_ctr > 0) {
    NumericVector lp_used(keep_ctr);
    for (int t = 0; t < keep_ctr; ++t) lp_used[t] = lp_test_samp[t];
    lp_test_log_mean_exp = log_mean_exp_pred(lp_used);
  }
  double lp_test_exact_log_mean_exp = NA_REAL;
  if (do_prediction && exact_prediction && y_test.size() == J && keep_ctr > 0) {
    NumericVector lp_used(keep_ctr);
    for (int t = 0; t < keep_ctr; ++t) lp_used[t] = lp_test_exact_samp[t];
    lp_test_exact_log_mean_exp = log_mean_exp_pred(lp_used);
  }

  return List::create(
    Named("c0") = c0_samp,
    Named("c_view") = c_view_samp,
    Named("gamma") = gamma_samp,
    Named("alpha") = alpha_samp,
    Named("alpha_rb") = alpha_rb_samp,
    Named("mu") = mu_samp,
    Named("sigma2") = sigma2_samp,
    Named("mu0") = mu0_samp,
    Named("tau2") = tau2_samp,
    Named("y_new") = y_new_samp,
    Named("y_new_rb") = y_new_rb_samp,
    Named("y_in") = y_in_samp,
    Named("lp_test") = lp_test_samp,
    Named("lp_test_log_mean_exp") = lp_test_log_mean_exp,
    Named("lp_test_exact") = lp_test_exact_samp,
    Named("lp_test_exact_log_mean_exp") = lp_test_exact_log_mean_exp,
    Named("gamma_sum") = gamma_sum_samp,
    Named("prediction_state") = prediction_state_samp,
    Named("prediction_state_final") = List::create(
      Named("g_pred") = g_pred,
      Named("c0_pred") = c0_pred,
      Named("cview_pred") = cview_pred,
      Named("cprob_pred") = cprob_pred
    ),
    Named("final_state") = List::create(
      Named("c0") = c0,
      Named("c_view") = c_view,
      Named("gamma") = gamma,
      Named("alpha") = alpha,
      Named("mu") = mu,
      Named("sigma2") = sigma2
    )
  );
}
