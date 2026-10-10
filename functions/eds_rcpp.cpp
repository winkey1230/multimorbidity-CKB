// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <random>
#include <algorithm>
using namespace Rcpp;

//' EDS 核心匹配循环 (C++ 原生 RNG 版)
//'
//' 与 eds_match_rcpp 逻辑完全相同，优化点：
//' 1. CharacterVector 预分配 + 下标写入，消除 push_back 扩容（越跑越慢）
//' 2. valid_ctrl 复用 clear() + reserve()，消除每次迭代的堆分配
//' 3. C++ <random> 替代 Rcpp::sample()，消除跨语言 RNG 调用
//'
//' 注意：C++ <random> 与 R 的 RNG 路径不同，即使同 seed 也会产生
//' 不同的随机样本，但统计推断等价。
//'
//' @keywords internal
// [[Rcpp::export]]
List eds_match_rcpp_v4(
    List ctrl_time_list,
    List ctrl_row_list,
    List ctrl_Ctime_list,
    NumericVector unique_A_time,
    IntegerVector unique_ctrl_group,
    List exp_rows_list,
    CharacterVector id_vec,
    int ratio,
    double window_time,
    int c_condition,
    bool trace,
    int seed
) {
  int n_groups = unique_A_time.size();
  bool has_C = (c_condition != 0);

  // 估算总输出行数
  int est_rows = 0;
  for (int i = 0; i < n_groups; i++) {
    int n_exp = as<IntegerVector>(exp_rows_list[i]).size();
    est_rows += n_exp * (ratio + 1);
  }

  // CharacterVector 预分配：消除 push_back 的反复扩容
  CharacterVector all_ids(est_rows);
  std::vector<double> all_starts;
  std::vector<int> all_exp;
  std::vector<int> group_sizes;
  all_starts.reserve(est_rows);
  all_exp.reserve(est_rows);
  group_sizes.reserve(n_groups);

  int pos = 0;  // 当前写入下标

  // C++ 原生 RNG：用 seed 初始化，保证可复现
  std::mt19937 gen(seed);

  // valid_ctrl 复用：提至循环外，避免每次迭代堆分配
  std::vector<int> valid_ctrl;
  valid_ctrl.reserve(65536);

  for (int i = 0; i < n_groups; i++) {
    double start_time_i = unique_A_time[i];
    int ctrl_group_i = unique_ctrl_group[i];
    IntegerVector exp_rows = as<IntegerVector>(exp_rows_list[i]);
    int n_exp = exp_rows.size();

    // 写入暴露组（下标写入，无扩容）
    for (int j = 0; j < n_exp; j++) {
      int ridx = exp_rows[j] - 1;
      all_ids[pos] = id_vec[ridx];
      all_starts.push_back(start_time_i);
      all_exp.push_back(1);
      pos++;
    }

    int n_sample = 0;

    bool has_ctrl_group = !IntegerVector::is_na(ctrl_group_i) && ctrl_group_i != 0;
    if (has_ctrl_group) {
      NumericVector tvec = ctrl_time_list[ctrl_group_i - 1];
      IntegerVector row_vec = ctrl_row_list[ctrl_group_i - 1];

      double threshold = start_time_i + window_time;
      int lo = 0, hi = tvec.size();
      while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (tvec[mid] <= threshold) lo = mid + 1; else hi = mid;
      }

      if (lo < tvec.size()) {
        // 复用 valid_ctrl，只清空不释放内存
        valid_ctrl.clear();
        if (has_C) {
          IntegerVector ctime_vec = ctrl_Ctime_list[ctrl_group_i - 1];
          for (int j = lo; j < row_vec.size(); j++) {
            bool ok = (c_condition == 1) ? (ctime_vec[j] <= start_time_i)
                                         : (ctime_vec[j] > start_time_i);
            if (ok) valid_ctrl.push_back(row_vec[j]);
          }
        } else {
          for (int j = lo; j < row_vec.size(); j++) {
            valid_ctrl.push_back(row_vec[j]);
          }
        }

        int n_ctrl = (int)valid_ctrl.size();
        n_sample = std::min(n_exp * ratio, n_ctrl);

        if (n_sample > 0) {
          if (n_ctrl <= 1) {
            int ridx = valid_ctrl[0] - 1;
            all_ids[pos] = id_vec[ridx];
            all_starts.push_back(start_time_i);
            all_exp.push_back(0);
            pos++;
          } else {
            // C++ 原生随机采样：shuffle 后取前 n_sample 个
            std::shuffle(valid_ctrl.begin(), valid_ctrl.end(), gen);
            for (int j = 0; j < n_sample; j++) {
              int ridx = valid_ctrl[j] - 1;
              all_ids[pos] = id_vec[ridx];
              all_starts.push_back(start_time_i);
              all_exp.push_back(0);
              pos++;
            }
          }
        }
      }
    }

    if (n_sample == 0) {
      all_ids[pos] = NA_STRING;
      all_starts.push_back(start_time_i);
      all_exp.push_back(0);
      n_sample = 1;
      pos++;
    }

    group_sizes.push_back(n_exp + n_sample);

    if (trace && (i + 1) % 100 == 0)
      Rcpp::Rcout << round((i + 1.0) / n_groups * 100) << "% ";
  }

  // 截断到实际长度（预分配可能偏大）
  int n_total = pos;
  CharacterVector id_out(n_total);
  for (int i = 0; i < n_total; i++) id_out[i] = all_ids[i];

  NumericVector start_out(n_total);
  IntegerVector exp_out(n_total);
  for (int i = 0; i < n_total; i++) {
    start_out[i] = all_starts[i];
    exp_out[i]   = all_exp[i];
  }

  return List::create(
    _["id_out"]      = id_out,
    _["start_out"]   = start_out,
    _["exp_out"]     = exp_out,
    _["n_total"]     = n_total,
    _["group_sizes"] = group_sizes
  );
}
