#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

static double quantile_cpp(std::vector<double>& v, double prob) {
  int n = v.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return v[0];
  std::sort(v.begin(), v.end());
  double index = 1.0 + (n - 1.0) * prob;
  int lo = (int)floor(index) - 1, hi = (int)ceil(index) - 1;
  if (lo < 0) lo = 0; if (hi >= n) hi = n - 1;
  double gamma = index - floor(index);
  if (gamma < 1e-15) gamma = 0.0;
  return (1.0 - gamma) * v[lo] + gamma * v[hi];
}

// [[Rcpp::export]]
List trajectory_kernel_cpp(
    NumericMatrix date_mat, IntegerMatrix event_mat,
    List base_idx_list, List paths, int n_minimum
) {
  int n_traj = paths.size(), D = date_mat.ncol();

  std::vector<std::vector<int>> base_cpp(D);
  for (int c = 0; c < D; c++) {
    SEXP elem = base_idx_list[c];
    if (elem == R_NilValue) continue;
    IntegerVector r_idx(elem);
    int sz = r_idx.size();
    std::vector<int>& v = base_cpp[c];
    v.reserve(sz);
    for (int i = 0; i < sz; i++) v.push_back(r_idx[i] - 1);
  }

  IntegerMatrix n_mat(n_traj, 6);
  NumericMatrix interval_stats(n_traj, 30);
  IntegerVector n_valid(n_traj);
  List raw_intervals(n_traj);
  std::fill(n_mat.begin(), n_mat.end(), NA_INTEGER);
  std::fill(interval_stats.begin(), interval_stats.end(), NA_REAL);
  std::fill(n_valid.begin(), n_valid.end(), 0);

  std::vector<int> idx_all[6];
  std::vector<double> intervals[5];

  for (int t = 0; t < n_traj; t++) {
    IntegerVector path_int = paths[t];
    int L = path_int.size();
    int d_idx[6];
    for (int k = 0; k < L; k++) d_idx[k] = path_int[k];

    const std::vector<int>& base = base_cpp[d_idx[0]];
    if (base.empty()) continue;

    for (int k = 0; k < L; k++) idx_all[k].clear();
    std::vector<int>& n1_vec = idx_all[0];
    int bsz = base.size();
    n1_vec.reserve(bsz);
    for (int bi = 0; bi < bsz; bi++) {
      int i = base[bi];
      bool ok = true;
      for (int j = 1; j < L; j++) {
        int dj = d_idx[j];
        if (event_mat(i, dj) == 1) {
          double dj_date = date_mat(i, dj);
          if (!std::isnan(dj_date) && dj_date <= date_mat(i, d_idx[0])) { ok = false; break; }
        }
      }
      if (ok) n1_vec.push_back(i);
    }
    if ((int)n1_vec.size() < n_minimum) continue;

    bool skip = false;
    for (int k = 1; k < L; k++) {
      int dk = d_idx[k], dk_1 = d_idx[k - 1];
      const std::vector<int>& prev = idx_all[k - 1];
      std::vector<int>& next = idx_all[k];
      next.clear(); next.reserve(prev.size() / 3);
      int psz = prev.size();
      for (int pi = 0; pi < psz; pi++) {
        int i = prev[pi];
        if (event_mat(i, dk) != 1) continue;
        double dk_date = date_mat(i, dk);
        if (std::isnan(dk_date)) continue;
        if (dk_date <= date_mat(i, dk_1)) continue;
        bool ok = true;
        for (int j = k + 1; j < L; j++) {
          int dj = d_idx[j];
          if (event_mat(i, dj) == 1) {
            double dj_date = date_mat(i, dj);
            if (!std::isnan(dj_date) && dj_date <= dk_date) { ok = false; break; }
          }
        }
        if (ok) next.push_back(i);
      }
      if ((int)next.size() < n_minimum) { skip = true; break; }
    }
    if (skip) continue;

    for (int s = 0; s < L - 1; s++) intervals[s].clear();
    for (int s = 0; s < L - 1; s++) {
      const std::vector<int>& ids = idx_all[s + 1];
      int n_ids = ids.size();
      std::vector<double>& intv = intervals[s];
      intv.reserve(n_ids);
      int ds_now = d_idx[s], ds_next = d_idx[s + 1];
      for (int pi = 0; pi < n_ids; pi++) {
        int i = ids[pi];
        intv.push_back((date_mat(i, ds_next) - date_mat(i, ds_now)) / 365.25);
      }
      if (n_ids > 0) {
        int cb = s * 5;
        interval_stats(t, cb + 0) = quantile_cpp(intv, 0.50);
        interval_stats(t, cb + 1) = quantile_cpp(intv, 0.25);
        interval_stats(t, cb + 2) = quantile_cpp(intv, 0.75);
        interval_stats(t, cb + 3) = quantile_cpp(intv, 0.025);
        interval_stats(t, cb + 4) = quantile_cpp(intv, 0.975);
      }
    }

    for (int k = 0; k < L; k++) n_mat(t, k) = idx_all[k].size();
    n_valid[t] = idx_all[L - 1].size();

    int total_sz = 0;
    for (int s = 0; s < L - 1; s++) total_sz += intervals[s].size();
    NumericVector flat_intv(total_sz);
    IntegerVector sizes_vec(L - 1);
    int offset = 0;
    for (int s = 0; s < L - 1; s++) {
      sizes_vec[s] = intervals[s].size();
      for (size_t kk = 0; kk < intervals[s].size(); kk++)
        flat_intv[offset + kk] = intervals[s][kk];
      offset += intervals[s].size();
    }
    raw_intervals[t] = List::create(_["v"] = flat_intv, _["sz"] = sizes_vec);
  }

  return List::create(
    _["n_mat"] = n_mat, _["interval_stats"] = interval_stats,
    _["raw_intervals"] = raw_intervals, _["n_valid"] = n_valid
  );
}
