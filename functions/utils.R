###############################################################################-#
## 共享工具函数库
## 用于暴露密度抽样(EDS)下的疾病对效应估算
##
## 本文件被 comorbidity-fun.R 和 trajectory-driven-fun.R 共同 source
## 用户不应直接加载本文件
##
## 依赖: tidyverse, gtools, doSNOW, data.table, survival, glmnet
##
## ==== 方法学概述 ====
##
## 暴露密度抽样 (Exposure Density Sampling, EDS) 是一种嵌套病例-对照设计，
## 用于时依暴露的时间-to-event 数据分析。
##
## 核心流程:
## 1. 倾向性评分 (PS) 估计: 通过分层 Cox 模型估计"发生暴露A"的倾向性
## 2. 风险集抽样: 在每个暴露事件时间点，从风险集中按匹配变量抽样对照
## 3. 匹配后处理: 定义随访区间、构建时依协变量 (tmerge)
## 4. 效应估计: LASSO-Cox 惩罚回归 + 后选择未惩罚 Cox 模型
##
## 五种模型变体:
## - Unadjusted: 仅包含暴露变量，不调整其他疾病
## - AIC-lambda: LASSO 选择 lambda 使 AIC 最小
## - BIC-lambda: LASSO 选择 lambda 使 BIC 最小（更保守）
## - CV-min:     交叉验证选择使偏差最小的 lambda
## - CV-1se:     交叉验证选择 1 标准误范围内的最正则化 lambda
##
## LASSO 策略: 暴露变量始终不惩罚 (penalty.factor = 0)，其他疾病被 L1 惩罚，
## 实现自动变量选择。
###############################################################################-#


# ---- 连续变量分箱 ----

#' 按固定带宽将连续变量划分为类别
#'
#' @param x   一个数值向量，待离散化的连续变量
#' @param bandwidth 每个区间的宽度
#' @return 一个整数向量，表示每个观测值所属的类别编号
#' @keywords internal
cut_cat <- function(x, bandwidth) {
  breaks_seq <- seq(from = min(x), to = max(x), by = bandwidth)
  breaks_seq[which.max(breaks_seq)] <- max(x)
  cat_value <- cut(x, breaks = breaks_seq, include.lowest = TRUE) %>% as.integer
  cat_value
}


# ---- LASSO 信息准则 ----

#' 从 glmnet 拟合结果中按 AIC 和 BIC 选择最优 lambda
#'
#' @param fit 一个 glmnet 或 cv.glmnet 拟合对象
#' @param n   整数，有效样本量（事件数），用于 BIC 惩罚项
#' @return 一个列表，包含 lambda.AIC 和 lambda.BIC
#' @keywords internal
extractAICBIC.glmnet <- function(fit, n) {
  df_seq <- fit$df
  deviance <- fit$nulldev * (1 - fit$dev.ratio)
  aic <- deviance + 2 * df_seq
  bic <- deviance + df_seq * log(n)
  return(list(lambda.AIC = fit$lambda[which.min(aic)],
              lambda.BIC = fit$lambda[which.min(bic)]))
}


# ---- LASSO 超参数自适应选择 ----

#' 根据匹配后样本量自适应选择 LASSO 超参数
#'
#' 在计算速度与精度之间取得平衡。样本量越大，加速越激进
#' （更少的 lambda 值、更高的最小 lambda 比率、更少的 CV 折数）。
#'
#' @param n_total 整数，匹配后数据集的总行数
#' @return 一个列表，包含: nlambda_val, lambda_min_ratio_val,
#'   nfolds_val, control_list（传给 glmnet 的 control 参数列表）
#' @keywords internal
select_lasso_parameters <- function(n_total) {
  if (n_total > 100000) {
    # 超大样本加速（10万+）
    nlambda_val <- 18
    lambda_min_ratio_val <- 0.12
    nfolds_val <- 5
    mxitnr_val <- 12
    epsnr_val <- 1e-4
  } else if (n_total > 50000) {
    # 极限加速（5万+样本）
    nlambda_val <- 22
    lambda_min_ratio_val <- 0.1
    nfolds_val <- 5
    mxitnr_val <- 15
    epsnr_val <- 1e-4
  } else if (n_total > 20000) {
    # 激进加速（2-5万样本）
    nlambda_val <- 30
    lambda_min_ratio_val <- 0.08
    nfolds_val <- 5
    mxitnr_val <- 15
    epsnr_val <- 1e-4
  } else if (n_total > 5000) {
    # 中度加速（5千-2万样本）
    nlambda_val <- 35
    lambda_min_ratio_val <- 0.05
    nfolds_val <- 8
    mxitnr_val <- 18
    epsnr_val <- 1e-4
  } else if (n_total > 1000) {
    # 轻度加速（1千-5千样本）
    nlambda_val <- 40
    lambda_min_ratio_val <- 0.03
    nfolds_val <- 10
    mxitnr_val <- 20
    epsnr_val <- 5e-5
  } else {
    # 小样本保守（<1000样本）
    nlambda_val <- 50
    lambda_min_ratio_val <- 0.01
    nfolds_val <- 10
    mxitnr_val <- 25
    epsnr_val <- 1e-5
  }
  control_list <- list(
    maxit   = 2000,
    mxitnr  = mxitnr_val,
    epsnr   = epsnr_val,
    fdev    = 1e-4,
    thresh  = 1e-5,   # relax coordinate-descent convergence (default 1e-7)
    pmax    = 60,     # limit max nonzero variables (default: no limit)
    itrace  = 0
  )
  return(list(
    nlambda_val = nlambda_val,
    lambda_min_ratio_val = lambda_min_ratio_val,
    nfolds_val = nfolds_val,
    control_list = control_list
  ))
}


# ---- 匹配变量处理 ----

#' 对模糊匹配变量进行分箱处理
#'
#' 对每个模糊匹配变量，使用 cut_cat 按指定的 caliper 带宽生成分类版本。
#'
#' @param data           数据框，包含待处理的匹配变量
#' @param matchvar_fuzzy 字符向量，模糊匹配变量名
#' @param caliper_fuzzy  数值向量，与 matchvar_fuzzy 一一对应的 caliper 带宽
#' @return 添加了 _cat 列的数据框
#' @keywords internal
build_matchvar_categories <- function(data, matchvar_fuzzy, caliper_fuzzy) {
  for (i in seq_along(matchvar_fuzzy)) {
    data[, paste0(matchvar_fuzzy[i], "_cat")] <- cut_cat(
      data[, matchvar_fuzzy[i]], caliper_fuzzy[i]
    )
  }
  return(data)
}


# ---- 倾向性评分 ----

#' 计算倾向性评分（Propensity Score）
#'
#' 通过分层 Cox 模型估计暴露 A 的倾向性评分。
#' 模型：Surv(ps_time_var, A_event) ~ covariates_ps + strata(matchvar_exact) + strata(matchvar_fuzzy_cat)
#'
#' @param data              数据框，包含所有变量
#' @param covariates_ps     字符向量，PS 模型的协变量名
#' @param matchvar_exact    字符向量，精确匹配变量名
#' @param matchvar_fuzzy_cat 字符向量，模糊匹配变量的分类列名（_cat 后缀）
#' @param caliper_ps_sd     数值标量，PS 匹配的 caliper（以 PS 的 SD 为单位）
#' @param ps_time_var       字符，PS 模型中 Surv 的时间变量名（"A_time" 或 "A_time_o"）
#' @param ps_filter_expr    字符，用于筛选 PS 拟合数据的表达式（如 "B_time >= A_time"）
#' @return 添加了 ps 和 ps_cat 列的数据框
#' @keywords internal
calculate_propensity_score <- function(data, covariates_ps, matchvar_exact,
                                        matchvar_fuzzy_cat, caliper_ps_sd,
                                        ps_time_var = "A_time",
                                        ps_filter_expr = NULL) {
  # 筛选用于拟合 PS 的数据
  data_fit_ps <- if (!is.null(ps_filter_expr)) {
    subset(data, eval(parse(text = ps_filter_expr)))
  } else {
    data
  }

  # 构建 PS 公式（条件拼接，支持 matchvar_exact 和 matchvar_fuzzy_cat 为空）
  cov_part <- paste0(covariates_ps, collapse = "+")
  strata_vars <- c(matchvar_exact, matchvar_fuzzy_cat)
  if (length(strata_vars) > 0) {
    strata_part <- paste0(" + ", paste0(paste0("strata(", strata_vars, ")"), collapse = "+"))
  } else {
    strata_part <- ""
  }
  rhs <- if (nchar(cov_part) > 0) paste0(cov_part, strata_part) else "1"
  ps_formula <- as.formula(paste0("Surv(", ps_time_var, ", A_event) ~ ", rhs))

  fit_ps <- coxph(ps_formula, data = data_fit_ps)
  data$ps <- predict(fit_ps, newdata = data, type = "lp", reference = "sample")
  caliper_val <- caliper_ps_sd * sd(data$ps)
  names(caliper_val) <- "ps"
  data$ps_cat <- cut_cat(data$ps, caliper_val)
  return(data)
}


# ---- 编译/加载 Rcpp 加速模块 ----
# fork 模式（mclapply）下子进程继承父进程的已编译函数，跳过重复编译
# socket 模式下每个 worker 独立编译一次
local({
  if (exists("eds_match_rcpp_v4", envir = .GlobalEnv, inherits = FALSE)) {
    return(invisible(NULL))
  }

  src_file <- sys.frame(1)$ofile
  rcpp_path <- if (!is.null(src_file) && nchar(src_file) > 0L) {
    file.path(dirname(normalizePath(src_file)), "eds_rcpp.cpp")
  } else {
    "eds_rcpp.cpp"
  }
  # Fallback: also try version4/eds_rcpp.cpp
  if (!file.exists(rcpp_path)) rcpp_path <- "version4/eds_rcpp.cpp"
  if (file.exists(rcpp_path)) Rcpp::sourceCpp(rcpp_path)
})

# ---- 暴露密度抽样 (EDS) ----

#' 执行暴露密度抽样（EDS）匹配
#'
#' 对每个唯一的 (A_time, 匹配变量) 组合，从该时间点的风险集中按匹配变量
#' 抽样对照。这是 EDS 方法的核心步骤。
#'
#' ## EDS 方法学
#'
#' EDS 是一种嵌套病例-对照抽样设计，用于时依暴露的时间-to-event 数据：
#' 1. **风险集形成**：对于时间 t 的每个暴露事件，风险集由在时间 t
#'    仍处于结局风险中（即尚未发生结局或被删失）的所有个体组成。
#' 2. **匹配**：在每个风险集中，对照按以下条件与暴露个体匹配：
#'    - 精确匹配变量（如性别、地区）
#'    - 模糊匹配变量（如年龄、BMI），通过 caliper 分箱
#'    - 倾向性评分分箱（由时间-to-暴露的 Cox 模型估计）
#'
#' @param data                  data.table，包含 A_time, B_time, matchvar_names 等列
#' @param matchvar_names        字符向量，匹配变量列名（含 _cat 和 ps_cat）
#' @param ratio                 整数，每个暴露个体的对照数
#' @param window_time           数值（天），A 和 B 之间的最小时间间隔
#' @param seed                  整数，随机种子
#' @param exposure_filter_expr  字符，筛选暴露组的表达式（如 "A_time + window_time < B_time & A_event == 1"）
#' @param control_filter_expr   字符，筛选对照候选的表达式（使用 start_time_i 表示当前索引时间）
#' @param id                    字符，个体标识符列名
#' @param include_C_time        逻辑，control_dt 是否包含 C_time 列（用于 Ctrue/Cfalse 的对照筛选）
#' @param trace                 逻辑，每100组打印进度
#' @return 一个列表，包含: data_matched_list（匹配数据列表）和 exposure_dt（暴露data.table）
#' @keywords internal
run_eds_sampling <- function(data, matchvar_names, ratio, window_time, seed,
                              exposure_filter_expr, control_filter_expr,
                              id, include_C_time = FALSE, trace = FALSE) {
  set.seed(seed)
  data <- as.data.table(data)
  data[, row_number := .I]

  # ---- 预解析筛选表达式（避免在循环内反复 parse） ----
  exposure_filter_call <- parse(text = exposure_filter_expr)[[1]]

  # ---- 自动探测 C 条件方向（用于窄集上的向量化过滤） ----
  c_condition <- "none"
  if (include_C_time) {
    if (grepl("C_time\\s*<=\\s*start_time_i", control_filter_expr)) {
      c_condition <- "lte"
    } else if (grepl("C_time\\s*>\\s*start_time_i", control_filter_expr)) {
      c_condition <- "gt"
    }
  }

  # ---- 预计算列名向量（避免 .. 作用域问题） ----
  exp_cols  <- c("A_time", matchvar_names, "row_number")
  uniq_cols <- c("A_time", matchvar_names)
  map_cols  <- c(matchvar_names, "ctrl_group")
  ctrl_cols <- if (include_C_time) {
    c("A_time", "B_time", "C_time", matchvar_names, "row_number")
  } else {
    c("A_time", "B_time", matchvar_names, "row_number")
  }

  # ---- 构建暴露池 ----
  exposure_dt <- data[eval(exposure_filter_call),
    exp_cols, with = FALSE
  ][, group_id := .GRP, by = c("A_time", matchvar_names)]

  # ---- 构建对照池，按匹配变量分组 + 按 pmin(A_time, B_time) 排序 ----
  control_dt <- data[, ctrl_cols, with = FALSE]
  # pmin(A_time, B_time) > t  ⇔  A_time > t & B_time > t
  # 预排序后可用二分查找定位阈值，避免每轮全表扫描
  control_dt[, `:=`(
    time_min  = pmin(A_time, B_time),
    ctrl_group = .GRP
  ), by = matchvar_names]
  setkeyv(control_dt, c("ctrl_group", "time_min"))

  # ---- 按 ctrl_group 拆分对照池（每个 group 内按 time_min 升序） ----
  ctrl_split_cols <- if (include_C_time) {
    c("row_number", "time_min", "C_time")
  } else {
    c("row_number", "time_min")
  }
  ctrl_dt_list <- split(
    control_dt[, ctrl_split_cols, with = FALSE],
    control_dt$ctrl_group
  )
  # 预提取 time_min 向量列表（加速 findInterval）
  ctrl_time_list <- lapply(ctrl_dt_list, `[[`, "time_min")

  # ---- 暴露组 → 对照组的映射（一次性 join） ----
  unique_exposure_dt <- unique(exposure_dt[, uniq_cols, with = FALSE])
  ctrl_group_map <- unique(
    control_dt[, map_cols, with = FALSE]
  )
  unique_exposure_dt <- merge(
    unique_exposure_dt, ctrl_group_map,
    by = matchvar_names, all.x = TRUE, sort = FALSE
  )
  n_exposure_group <- nrow(unique_exposure_dt)

  # ---- 预提取（同上一版优化） ----
  exposure_rows_by_group <- split(exposure_dt$row_number, exposure_dt$group_id)
  id_vec <- data[[id]]

  data_matched_list <- vector("list", n_exposure_group)

  for (idx in seq_len(n_exposure_group)) {
    tryCatch({
      start_time_i <- unique_exposure_dt[idx, A_time]
      ctrl_group_i <- unique_exposure_dt[idx, ctrl_group]
      exposure_row_i <- exposure_rows_by_group[[as.character(idx)]]

      # ---- 对照筛选：二分查找 + 向量化 C 条件 ----
      if (is.na(ctrl_group_i)) {
        control_row_i <- NA_integer_
      } else {
        tvec <- ctrl_time_list[[ctrl_group_i]]
        # findInterval: 最后一个 <= threshold 的位置，+1 = 第一个 > threshold
        lo <- findInterval(start_time_i + window_time, tvec) + 1L
        if (lo <= length(tvec)) {
          grp <- ctrl_dt_list[[ctrl_group_i]]
          cand <- grp[lo:.N]  # 仅取时间达标行（已是 matchvar 匹配的）
          if (c_condition == "lte") {
            cand <- cand[C_time <= start_time_i]
          } else if (c_condition == "gt") {
            cand <- cand[C_time > start_time_i]
          }
          control_row_i <- sort(cand$row_number)
          # 保持与原版一致：无匹配时产生 NA（供 caller anyNA 剔除）
          if (length(control_row_i) == 0L) control_row_i <- NA_integer_
        } else {
          control_row_i <- NA_integer_
        }
      }

      control_row_sampling_i <- if (length(control_row_i) <= 1L) {
        control_row_i
      } else {
        sample(control_row_i,
               size = min(length(exposure_row_i) * ratio, length(control_row_i)))
      }

      out <- data.table(
        id_val     = c(id_vec[exposure_row_i], id_vec[control_row_sampling_i]),
        start_time = start_time_i,
        exposure   = c(rep(1L, length(exposure_row_i)),
                       rep(0L, length(control_row_sampling_i)))
      )
      setnames(out, "id_val", id)
      data_matched_list[[idx]] <- out
    }, error = function(e) {
      message(id, " ", paste(id_vec[exposure_row_i], collapse = ","), " ", e$message)
    })
    if (trace && idx %% 100 == 0) {
      cat(round(idx / n_exposure_group * 100), "%", " ")
    }
  }

  return(list(
    data_matched_list = data_matched_list,
    exposure_dt = exposure_dt
  ))
}


# ---- 疾病对解析 ----

#' 解析疾病对：从 D_all 自动生成排列对，或使用显式指定的 As/Bs
#'
#' 这是所有 batch 函数的公共入口逻辑。
#'
#' @param data       数据框，用于校验列名
#' @param D_all      字符向量或缺失。如果提供，生成所有有序对（排列）
#' @param As         字符向量，每个疾病对的第一个疾病
#' @param Bs         字符向量，每个疾病对的第二个疾病
#' @param adjusted_D 字符向量，待调整的疾病名
#' @param C          字符或 NULL，条件疾病（从 D_all 中排除）
#' @return 一个列表，包含: As, Bs, adjusted_D, n_pair,
#'   A_dates, B_dates, A_events, B_events, first_Ds, second_Ds
#' @keywords internal
resolve_disease_pairs <- function(data, D_all, As, Bs, adjusted_D, C = NULL) {
  if (!missing(D_all)) {
    if (!missing(As) || !missing(Bs)) {
      warning("The parameters As and Bs have been omitted")
    }
    if (any(duplicated(D_all))) stop("D_all have duplicated values")
    if (!is.null(C)) D_all <- setdiff(D_all, C)
    if (missing(adjusted_D)) adjusted_D <- D_all
    d_pairs <- gtools::permutations(length(D_all), 2, D_all)
    As <- d_pairs[, 1]
    Bs <- d_pairs[, 2]
  } else {
    if (missing(As) || missing(Bs)) {
      stop("These variables As and Bs cannot be missing")
    }
    if (length(As) != length(Bs)) {
      stop("These variables As and Bs must have the same length")
    }
    if (any(As == Bs)) stop("Disease A and B have the same value")
  }

  n_pair <- length(As)
  if (!all(c(paste0(As, "_date"), paste0(Bs, "_event")) %in% names(data))) {
    stop("Some D_event or D_date is not in data")
  }
  if (!is.null(C) && C %in% c(As, Bs)) {
    stop("C cannot be included in As or Bs")
  }

  return(list(
    As = As, Bs = Bs, adjusted_D = adjusted_D, n_pair = n_pair,
    A_dates = paste0(As, "_date"), B_dates = paste0(Bs, "_date"),
    A_events = paste0(As, "_event"), B_events = paste0(Bs, "_event"),
    first_Ds = As, second_Ds = Bs
  ))
}

# ---- Rcpp 加速版 EDS 抽样 ----
run_eds_sampling_rcpp <- function(data, matchvar_names, ratio, window_time,
                                   seed, exposure_filter_expr,
                                   control_filter_expr, id,
                                   include_C_time = FALSE, trace = FALSE) {
  set.seed(seed)
  data <- as.data.table(data)
  data[, row_number := .I]

  exposure_filter_call <- parse(text = exposure_filter_expr)[[1]]
  c_condition <- "none"
  if (include_C_time) {
    if (grepl("C_time\\s*<=\\s*start_time_i", control_filter_expr)) {
      c_condition <- "lte"
    } else if (grepl("C_time\\s*>\\s*start_time_i", control_filter_expr)) {
      c_condition <- "gt"
    }
  }
  c_cond_int <- switch(c_condition, none = 0L, lte = 1L, gt = 2L)

  exp_cols  <- c("A_time", matchvar_names, "row_number")
  uniq_cols <- c("A_time", matchvar_names)
  map_cols  <- c(matchvar_names, "ctrl_group")
  ctrl_cols <- if (include_C_time) {
    c("A_time", "B_time", "C_time", matchvar_names, "row_number")
  } else {
    c("A_time", "B_time", matchvar_names, "row_number")
  }

  exposure_dt <- data[eval(exposure_filter_call),
    exp_cols, with = FALSE
  ][, group_id := .GRP, by = c("A_time", matchvar_names)]

  control_dt <- data[, ctrl_cols, with = FALSE]
  control_dt[, `:=`(
    time_min  = pmin(A_time, B_time),
    ctrl_group = .GRP
  ), by = matchvar_names]
  setkeyv(control_dt, c("ctrl_group", "time_min"))

  ctrl_split_cols <- if (include_C_time) {
    c("row_number", "time_min", "C_time")
  } else {
    c("row_number", "time_min")
  }
  ctrl_dt_list <- split(
    control_dt[, ctrl_split_cols, with = FALSE],
    control_dt$ctrl_group
  )
  ctrl_time_list <- lapply(ctrl_dt_list, `[[`, "time_min")
  ctrl_row_list  <- lapply(ctrl_dt_list, `[[`, "row_number")
  ctrl_C_list    <- if (include_C_time) lapply(ctrl_dt_list, `[[`, "C_time") else list()

  unique_exposure_dt <- unique(exposure_dt[, uniq_cols, with = FALSE])
  ctrl_group_map <- unique(control_dt[, map_cols, with = FALSE])
  unique_exposure_dt <- merge(
    unique_exposure_dt, ctrl_group_map,
    by = matchvar_names, all.x = TRUE, sort = FALSE
  )

  exposure_rows_by_group <- split(exposure_dt$row_number, exposure_dt$group_id)
  id_vec <- data[[id]]                        # 保持原始类型，不改 RNG 状态
  id_is_int <- is.integer(id_vec)
  id_is_num <- is.numeric(id_vec) && !id_is_int
  ctrl_group_use <- unique_exposure_dt$ctrl_group
  ctrl_group_use[is.na(ctrl_group_use)] <- 0L

  # ---- C++ 加速（v4: 原生 RNG + 预分配 + 复用，比 v3 Rcpp::sample 快 ~2×） ----
  rcpp_res <- eds_match_rcpp_v4(
    ctrl_time_list    = ctrl_time_list,
    ctrl_row_list     = ctrl_row_list,
    ctrl_Ctime_list   = ctrl_C_list,
    unique_A_time     = unique_exposure_dt$A_time,
    unique_ctrl_group = ctrl_group_use,
    exp_rows_list     = exposure_rows_by_group,
    id_vec            = id_vec,
    ratio             = as.integer(ratio),
    window_time       = as.numeric(window_time),
    c_condition       = c_cond_int,
    trace             = isTRUE(trace),
    seed              = seed
  )

  # ---- 按 group_sizes 拆分回 data_matched_list ----
  id_out    <- rcpp_res$id_out
  start_out <- rcpp_res$start_out
  exp_out   <- rcpp_res$exp_out
  gsizes    <- rcpp_res$group_sizes
  if (id_is_int) { id_out <- as.integer(id_out); na_id <- NA_integer_ }
  else if (id_is_num) { id_out <- as.numeric(id_out); na_id <- NA_real_ }
  else na_id <- NA_character_

  n_groups <- nrow(unique_exposure_dt)
  data_matched_list <- vector("list", n_groups)
  pos <- 1L
  for (idx in seq_len(n_groups)) {
    exp_rows <- exposure_rows_by_group[[as.character(idx)]]
    n_exp <- length(exp_rows)
    gs <- if (idx <= length(gsizes)) gsizes[idx] else NA_integer_

    if (!is.na(gs) && gs > 0L && pos + gs - 1L <= length(id_out)) {
      end_pos <- pos + gs - 1L
      out <- data.table(id_val = id_out[pos:end_pos],
                        start_time = start_out[pos:end_pos],
                        exposure = exp_out[pos:end_pos])
      pos <- end_pos + 1L
    } else {
      # 回退：该组无对照或 gs 不可用
      out <- data.table(
        id_val     = c(id_vec[exp_rows], na_id),
        start_time = unique_exposure_dt[idx, A_time],
        exposure   = c(rep(1L, n_exp), 0L)
      )
    }
    setnames(out, "id_val", id)
    data_matched_list[[idx]] <- out
  }

  list(data_matched_list = data_matched_list, exposure_dt = exposure_dt)
}
