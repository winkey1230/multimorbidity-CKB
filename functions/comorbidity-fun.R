##############################################################################-#
## 该函数用于实现暴露密度抽样下的疾病对效应估算
## 校正基线变量和其它conditions
###############################################################################-#


library(tidyverse)
library(gtools)
library(doSNOW)
library(data.table)
library(survival)
library(glmnet)
library(Rcpp)
library(pbapply)

# Resolve utils.R relative to this script file, not the working directory
# Search through all frames for 'ofile' (set by source()), robust to nested sourcing
this_file <- NULL
for (i in sys.nframe():1) {
  env <- sys.frame(i)
  if (exists("ofile", envir = env, inherits = FALSE)) {
    this_file <- normalizePath(get("ofile", envir = env, inherits = FALSE))
    break
  }
}
if (!is.null(this_file)) {
  SCRIPT_DIR <- dirname(this_file)
} else {
  # Fallback: if run interactively (e.g., line-by-line), assume cwd is the project root
  SCRIPT_DIR <- file.path(getwd(), "version2")
}
source(file.path(SCRIPT_DIR, "utils.R"))

# ---- 系数提取辅助函数 ----

#' Extract coefficient row(s) from coxph fit without full summary()
#'
#' Avoids summary.coxph() overhead (robust variance, LR/Wald/score tests, R-squared)
#' by directly extracting coef, vcov diagonal, and computing z and p manually.
#'
#' @param fit  A coxph object
#' @param vars Character vector of coefficient names to extract
#' @return A matrix with columns: coef, exp(coef), se(coef), z, Pr(>|z|)
#' @keywords internal
extract_coxph_coef <- function(fit, vars) {
  b <- coef(fit)[vars]
  se <- sqrt(diag(vcov(fit))[vars])
  z <- b / se
  p <- 2 * pnorm(-abs(z))
  mat <- cbind(b, exp(b), se, z, p)
  colnames(mat) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  rownames(mat) <- vars
  mat
}

#' Estimate the causal effect of disease A on disease B
#'
#' Uses exposure density sampling (EDS) with propensity score matching
#' and LASSO-Cox regression to estimate whether disease A increases
#' the risk of disease B, adjusting for other conditions.
#'
#' @param data A data frame containing patient-level data. Must include
#'   an \code{id} column, and for each disease columns named
#'   \code{{Disease}_date} (Date) and \code{{Disease}_event} (0/1).
#' @param id Character. The column name of the unique patient identifier.
#' @param A Character. The name of the first disease (exposure).
#' @param B Character. The name of the second disease (outcome).
#' @param start_date Character. The column name in \code{data} containing
#'   the study start date (Date). Must be <= disease dates for all rows.
#' @param matchvar_exact Character vector. Variable names for exact
#'   matching (categorical variables). Use \code{NULL} to skip.
#' @param matchvar_fuzzy Character vector. Variable names for caliper
#'   matching (numeric variables only). Use \code{NULL} to skip.
#' @param caliper_fuzzy Numeric vector. Caliper bandwidths corresponding
#'   to \code{matchvar_fuzzy}. Must have the same length.
#' @param covariates_ps Character vector. Covariates for the propensity
#'   score Cox model. Use \code{NULL} to skip PS (all PS = 1).
#' @param caliper_ps_sd Numeric scalar. Caliper for PS matching, in units
#'   of SD of the propensity score. Default: 0.2.
#' @param ratio Integer. Number of controls to match per exposed case.
#'   Default: 1.
#' @param window_time Numeric. Minimum interval (days) between A and B
#'   onset. Only participants with A_time + window_time < B_time are
#'   considered exposed. Default: 0.
#' @param adjusted_D Character vector. Disease names to include as
#'   candidate confounders in LASSO-Cox selection.
#' @param output_match Logical. If \code{TRUE}, include the matched
#'   dataset in the returned list. Default: \code{FALSE}.
#' @param trace Logical. If \code{TRUE}, print EDS sampling progress
#'   every 100 risk sets. Default: \code{FALSE}.
#' @param cv_for_lambda Logical. If \code{TRUE}, perform cross-validation
#'   to compute CV-min and CV-1se lambda variants (slower).
#'   Default: \code{FALSE}.
#' @param seed Integer. Random seed for EDS sampling. Default: 1.
#' @return A list with components:
#'   \item{res_coef}{A 5-row matrix with rownames
#'     (coef_unadjusted_D, coef_AIC_lambda, coef_BIC_lambda,
#'     coef_cv_lambda.min, coef_cv_lambda.1se) and columns
#'     (coef, exp(coef), se(coef), z, Pr(>|z|)).}
#'   \item{n_exposure}{Number of exposed cases.}
#'   \item{n_control}{Number of matched controls.}
#'   \item{n_unmatched}{Number of exposed cases with no matched control.}
#'   \item{n_event}{Number of outcome events.}
#'   \item{data_matched}{Matched dataset (only when output_match = TRUE).}
#'
getAB <- function(data,id,A,B,adjusted_D,start_date,matchvar_exact = NULL,
                  matchvar_fuzzy = NULL,caliper_fuzzy = NULL,covariates_ps = NULL,
                  caliper_ps_sd = 0.2,
                  ratio = 1,window_time = 0,
                  output_match = FALSE,trace = FALSE, cv_for_lambda = F,seed = 1){
  # 参数标准化：NULL → 零长度向量
  matchvar_exact  <- if (is.null(matchvar_exact))  character(0) else matchvar_exact
  matchvar_fuzzy  <- if (is.null(matchvar_fuzzy))  character(0) else matchvar_fuzzy
  caliper_fuzzy   <- if (is.null(caliper_fuzzy))   numeric(0)   else caliper_fuzzy
  covariates_ps   <- if (is.null(covariates_ps))   character(0) else covariates_ps
  if(length(matchvar_fuzzy) != length(caliper_fuzzy)) stop("matchvar_fuzzy and caliper_fuzzy have different length")
  data <- as.data.frame(data)
  if(any(sapply(c(id,A,B,start_date,ratio,caliper_ps_sd),length) != 1))
    stop("The length of id,A,B,start_date,ratio,caliper_ps_sd must be 1")
  A_date <- paste0(A,"_date")
  B_date <- paste0(B,"_date")
  A_event <- paste0(A,"_event")
  B_event <- paste0(B,"_event")
  
  adjusted_D <- setdiff(adjusted_D,c(A,B))
  adjusted_D_date <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_date") else character(0)
  adjusted_D_event <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_event") else character(0)
  if (!start_date %in% names(data)) stop("start_date column not found in data")
  if (class(data[[start_date]]) != "Date") stop("start_date must be Date type")
  if (any(data[[start_date]] > data[[A_date]] | data[[start_date]] > data[[B_date]], na.rm = TRUE))
    stop("start_date must be <= A_date and B_date for all rows")
  
  keepvars <- c(id,A_date,A_event,B_date,B_event,start_date,matchvar_exact,
                matchvar_fuzzy,covariates_ps,adjusted_D_date,adjusted_D_event)
  if(!all(keepvars %in% names(data))) 
    stop("There were variables not being included in the data")
  if(any(sapply(data[,c(A_date,B_date,start_date,adjusted_D_date)], class) != "Date")) stop("There are date variables not being Date type")
  if(class(data[,A_event]) == "logical") data[,A_event] <- as.integer(data[,A_event])
  if(class(data[,B_event]) == "logical") data[,B_event] <- as.integer(data[,B_event])
  if(A == B) stop("A is the same as B")
  if(!all(sapply(data[,matchvar_fuzzy,drop = FALSE], class) %in% c("integer","numeric"))) stop("matchvar_fuzzy must be a numeric or integer variable")
  
  origin_data <- na.omit(data[,c(id,A_date,A_event,B_date,B_event,start_date,matchvar_exact,
                                 matchvar_fuzzy,covariates_ps,adjusted_D_date,adjusted_D_event)])  
  for (i in adjusted_D) {
    origin_data[,paste0(i,"_time")] <- difftime(origin_data[,paste0(i,"_date")],origin_data[,start_date],units = "days") %>% as.integer
  }
  
  data <- origin_data[,c(id,A_date,A_event,B_date,B_event,start_date,matchvar_exact,matchvar_fuzzy,covariates_ps)]  
  data[,"A_time"] <- difftime(data[,A_date],data[,start_date],units = "days") %>% as.integer
  data[,"B_time"] <- difftime(data[,B_date],data[,start_date],units = "days") %>% as.integer
  data[,"A_event"] <- data[,A_event]
  data[,"B_event"] <- data[,B_event]
  data <- subset(data,data$A_time > 0 & data$B_time > 0)

  # 根据matchvar_fuzzy划分类别并计算倾向性评分（支持空参数）
  if (length(matchvar_fuzzy) > 0) {
    data <- build_matchvar_categories(data, matchvar_fuzzy, caliper_fuzzy)
  }
  # 占位符：确保 matchvar_names 永不为空
  if (length(matchvar_exact) == 0) {
    data[["_placeholder_exact"]] <- 1L
  }
  if (length(matchvar_fuzzy) == 0) {
    data[["_placeholder_fuzzy_cat"]] <- 1L
  }
  matchvar_fuzzy_cat <- if (length(matchvar_fuzzy) > 0) paste0(matchvar_fuzzy, "_cat") else character(0)
  if (length(covariates_ps) > 0) {
    data <- calculate_propensity_score(data, covariates_ps, matchvar_exact,
                                        matchvar_fuzzy_cat, caliper_ps_sd,
                                        ps_time_var = "A_time",
                                        ps_filter_expr = "A_time > 0 & B_time > 0 & B_time >= A_time")
  } else {
    data[["ps"]] <- 1
    data[["ps_cat"]] <- 1L
  }


  # 用ccwc进行EDS抽样， 因为ccwc函数无法限制B_time和A_time的关系，以及无法灵活改变window_time参数，所以被弃用
  # set.seed(seed)
  # matchvar_names <- paste(c(matchvar_exact,paste0(matchvar_fuzzy,"_cat") ,"ps_cat"), collapse = ", ")
  # include_names <- paste0(c(id,matchvar_fuzzy,"ps","A_time","A_event","B_time","B_event"), collapse = ", ")
  # expr <- paste0(
  #   "ccwc(entry = 0, exit = A_time, fail = A_event, origin = 0, ",
  #   "controls = ratio, match = list(", matchvar_names, "), ",
  #   "include = list(", include_names, "), ",
  #   "data = data, silent = TRUE)"
  # )
  # matched_data_ccwc <- eval(parse(text = expr))
  
  ######## EDS抽样
  matchvar_names <- c(
    if (length(matchvar_exact) > 0) matchvar_exact else "_placeholder_exact",
    if (length(matchvar_fuzzy) > 0) paste0(matchvar_fuzzy, "_cat") else "_placeholder_fuzzy_cat",
    "ps_cat"
  )
  eds_result <- run_eds_sampling_rcpp(
    data = data, matchvar_names = matchvar_names,
    ratio = ratio, window_time = window_time, seed = seed,
    exposure_filter_expr = "A_time + window_time < B_time & A_event == 1",
    control_filter_expr = "A_time > start_time_i + window_time & B_time > start_time_i + window_time",
    id = id, include_C_time = FALSE, trace = trace
  )
  data_matched_list <- eds_result$data_matched_list
  exposure_dt <- eds_result$exposure_dt
  
  # 剔除没有匹配到对照的list
  NA_index <- sapply(data_matched_list, function(x) anyNA(x[[id]])) %>% which()
  if(length(NA_index) >= 1) data_matched_list <- data_matched_list[-NA_index]


  adjusted_D_time <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_time") else character(0)
  join_cols <- c(id, matchvar_names, matchvar_fuzzy, "ps", "A_time", "B_time", "A_event", "B_event")
  data_matched <- rbindlist(data_matched_list) %>%
    left_join(data[, join_cols], by = id) %>%
    left_join(origin_data[,c(id,adjusted_D_time,adjusted_D_event),drop = FALSE],by = id) %>%
    group_by(!!!syms(matchvar_names)) %>%
    mutate(strata = cur_group_id()) %>%
    ungroup()
  
  # 处理对照组中B发生在A之前的：只保留区间1，因为区间2已经在exposure组了 
  data_matched <- data_matched %>% mutate(
    end_time = fifelse(
      exposure == 0 & A_time < B_time, A_time,
      B_time),
    status = fifelse(
      exposure == 0 & A_time < B_time, 0L,
      B_event),
  )
  
  # 去掉strata内部的重复个体，保留start_time最早的
  data_matched <- data_matched %>% 
    group_by(strata,!!sym(id),exposure) %>%
    slice_min(start_time, n = 1, with_ties = FALSE) %>%  # 保留 start_time 最小的行
    ungroup()
  
  # 去掉没有事件发生的strata
  data_matched_model <- data_matched %>%
    group_by(strata) %>%
    filter(any(status == 1)) %>%  # 保留至少有一个事件的 strata
    ungroup()
  
  # 添加其它condition的时依协变量
  data_matched_model$unique_row <- 1:nrow(data_matched_model)
  tdc_part <- if (length(adjusted_D) > 0) {
    paste0(",", paste0(paste0("`", adjusted_D, "` = tdc(`", adjusted_D, "_time`)"), collapse = ","))
  } else {
    ""
  }
  tdc_expr <- paste0(
    "tmerge(data_matched_model, data_matched_model, id=unique_row,tstart = start_time,tstop = end_time,status = event(end_time,status)",
    tdc_part, ")"
  )
  td_matched_data <- eval(parse(text = tdc_expr)) 
  
  n_unmatched <- nrow(exposure_dt) - sum(data_matched$exposure == 1)
  if(n_unmatched > 0) warning(sprintf("There are %.0f participants having not been matched for %.0f participants",n_unmatched,nrow(exposure_dt)))
  n_exposure <- nrow(exposure_dt)
  n_control <- sum(data_matched$exposure==0)
  n_event <- sum(data_matched$status==1)
  #### fit model: cox + lasso

  # not adjusting for other conditions----------
  fit_unadjusted_D <- coxph(Surv(tstart,tstop,status) ~ exposure + strata(strata),data = td_matched_data)
  coef_unadjusted_D <- extract_coxph_coef(fit_unadjusted_D, "exposure")[1L, ]
  
  # adjusting for other conditions---------
  # ===== 根据样本量分层选择参数 
  n_total <- nrow(td_matched_data)
  params <- select_lasso_parameters(n_total)
  nlambda_val <- params$nlambda_val
  lambda_min_ratio_val <- params$lambda_min_ratio_val
  nfolds_val <- params$nfolds_val
  control_list <- params$control_list
  # browser()
  if (length(adjusted_D) > 0) {
  x0 <- td_matched_data[,c("exposure",adjusted_D)] %>% as.matrix()
  y0 <- stratifySurv(Surv(td_matched_data$tstart,td_matched_data$tstop,td_matched_data$status),strata = td_matched_data$strata)
  
  ## 因为glmnet只能从global environment调取参数，所以在嵌套函数中无法识别control_list等参数，
  ## 故采用do.call方法
  # fit_lasso <- glmnet(x = x0,y = y0,family = "cox",penalty.factor = c(0,rep(1, length(adjusted_D))),
  #                     relax = T,nlambda = nlambda_val,cox.ties = "breslow",
  #                     lambda.min.ratio = lambda_min_ratio_val,control = control_list)
  args_list <- list(
    x = x0,y = y0,family = "cox",penalty.factor = c(0,rep(1, length(adjusted_D))),
    relax = T,nlambda = nlambda_val,cox.ties = "breslow",
    lambda.min.ratio = lambda_min_ratio_val,control = control_list,
    standardize = FALSE
  )
  fit_lasso <- do.call(glmnet,args_list)
  
  
  lambda_aic_bic <- extractAICBIC.glmnet(fit_lasso,n_event)
  
  coef_lasso <- coef(fit_lasso,s = lambda_aic_bic$lambda.AIC)
  vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
  fit_cox <- coxph(cox_formula,data = td_matched_data)
  coef_AIC_lambda <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  
  coef_lasso <- coef(fit_lasso,s = lambda_aic_bic$lambda.BIC)
  vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
  fit_cox <- coxph(cox_formula,data = td_matched_data)
  coef_BIC_lambda <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  
  if(cv_for_lambda){
    # fit_lasso_cv <- cv.glmnet(x = x0,y = y0,family = "cox",penalty.factor = c(0,rep(1, length(adjusted_D))),
    #                           relax = T,nlambda = nlambda_val,cox.ties = "breslow",nfolds = nfolds_val,gamma = 0,
    #                           lambda.min.ratio = lambda_min_ratio_val,control = control_list)
    fit_lasso_cv <- do.call(cv.glmnet,modifyList(args_list,val = list(gamma = 0,nfolds = nfolds_val)))
    
    coef_lasso <- coef(fit_lasso_cv,s = fit_lasso_cv$lambda.min)
    vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
    fit_cox <- coxph(cox_formula,data = td_matched_data)
    coef_cv_lambda.min <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
    
    coef_lasso <- coef(fit_lasso_cv,s = fit_lasso_cv$lambda.1se)
    vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
    fit_cox <- coxph(cox_formula,data = td_matched_data)
    coef_cv_lambda.1se <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  } else{
    coef_cv_lambda.1se <- coef_cv_lambda.min <- NA
  }
    
  } else {
    coef_AIC_lambda <- coef_BIC_lambda <- coef_cv_lambda.min <- coef_cv_lambda.1se <- NA
  }
  res_coef <- rbind(coef_unadjusted_D,coef_AIC_lambda,coef_BIC_lambda,coef_cv_lambda.min,coef_cv_lambda.1se)
  
  if(output_match) {
    return(list(res_coef = res_coef,n_exposure = n_exposure,n_control = n_control,
                 n_unmatched = n_unmatched,n_event = n_event,data_matched = data_matched))
  } else {
    return(list(res_coef = res_coef,n_exposure = n_exposure,n_control = n_control,
                 n_unmatched = n_unmatched,n_event = n_event))
  }
}



#' Batch estimation of disease-pair causal effects
#'
#' Runs \code{getAB()} for multiple disease pairs in parallel. Disease pairs
#' are either generated as all ordered permutations of \code{D_all}, or
#' specified explicitly via \code{As} and \code{Bs}.
#'
#' @param data A data frame containing patient-level data.
#' @param id Character. The column name of the unique patient identifier.
#' @param D_all Character vector. All disease names. If provided (and
#'   \code{As}/\code{Bs} are missing), all ordered pairs are generated
#'   via \code{gtools::permutations()}.
#' @param As Character vector. First disease of each pair. Required if
#'   \code{D_all} is missing. Must have the same length as \code{Bs}.
#' @param Bs Character vector. Second disease of each pair. Required if
#'   \code{D_all} is missing.
#' @param start_date Character. Column name of the study start date (Date).
#' @param matchvar_exact Character vector. Exact matching variables.
#'   Use \code{NULL} to skip.
#' @param matchvar_fuzzy Character vector. Fuzzy matching variables
#'   (numeric only). Use \code{NULL} to skip.
#' @param caliper_fuzzy Numeric vector. Caliper bandwidths for
#'   \code{matchvar_fuzzy}. Must have the same length.
#' @param covariates_ps Character vector. Covariates for the PS model.
#'   Use \code{NULL} to skip.
#' @param caliper_ps_sd Numeric scalar. PS caliper in SD units.
#'   Default: 0.2.
#' @param ratio Integer. Controls per exposed case. Default: 1.
#' @param window_time Numeric. Minimum interval (days) between A and B.
#'   Default: 0.
#' @param adjusted_D Character vector. Disease names as candidate
#'   confounders for LASSO-Cox.
#' @param ncore Integer. Number of CPU cores for parallel execution.
#'   Default: 1.
#' @param cv_for_lambda Logical. If \code{TRUE}, compute CV-based lambda
#'   variants. Default: \code{FALSE}.
#' @param seed Integer. Base random seed. Each pair receives \code{seed + i}
#'   for reproducibility. Default: 1.
#' @param parallel_backend Character. Parallel backend to use:
#'   \code{"auto"} (default, detects OS), \code{"fork"} (Linux/macOS,
#'   shared memory), or \code{"socket"} (Windows).
#' @param stop_on_error Logical. If \code{FALSE} (default), failing pairs
#'   return \code{NA} rows instead of aborting the entire batch.
#' @param ... Additional arguments (currently unused).
#' @return A data frame with columns: \code{coef_type}, \code{first_D},
#'   \code{second_D}, \code{n_exposure}, \code{n_unmatched},
#'   \code{n_control}, \code{n_event}, \code{coef}, \code{exp(coef)},
#'   \code{se(coef)}, \code{z}, \code{Pr(>|z|)}.
getAB_batch <- function(data,id,D_all,As,Bs,start_date,matchvar_exact = NULL,
                               matchvar_fuzzy = NULL,caliper_fuzzy = NULL,covariates_ps = NULL,caliper_ps_sd = 0.2,
                               ratio = 1,window_time = 0,adjusted_D,ncore = 1,
                               cv_for_lambda = F,seed = 1,
                               parallel_backend = c("auto", "fork", "socket"),
                               stop_on_error = FALSE,...){
  # 参数标准化
  matchvar_exact  <- if (is.null(matchvar_exact))  character(0) else matchvar_exact
  matchvar_fuzzy  <- if (is.null(matchvar_fuzzy))  character(0) else matchvar_fuzzy
  caliper_fuzzy   <- if (is.null(caliper_fuzzy))   numeric(0)   else caliper_fuzzy
  covariates_ps   <- if (is.null(covariates_ps))   character(0) else covariates_ps
  if(length(matchvar_fuzzy) != length(caliper_fuzzy)) stop("matchvar_fuzzy and caliper_fuzzy have different length")

  # 解析疾病对
  pairs <- resolve_disease_pairs(data, D_all, As, Bs, adjusted_D)
  As <- pairs$As; Bs <- pairs$Bs; adjusted_D <- pairs$adjusted_D
  n_pair <- pairs$n_pair
  A_dates <- pairs$A_dates; B_dates <- pairs$B_dates
  A_events <- pairs$A_events; B_events <- pairs$B_events
  first_Ds <- pairs$first_Ds; second_Ds <- pairs$second_Ds

  # 选择并行后端
  parallel_backend <- match.arg(parallel_backend)
  if (parallel_backend == "auto") {
    parallel_backend <- if (.Platform$OS.type == "windows") "socket" else "fork"
  }

  if (parallel_backend == "fork") {
    # Linux/macOS: fork 模式，共享内存，无需序列化 data
    ncore <- min(ncore, parallel::detectCores() - 1L, n_pair)
    results <- pbapply::pblapply(seq_len(n_pair), function(i) {
      tryCatch({
        fit <- getAB(data, id, As[i], Bs[i],
          adjusted_D,
          start_date, matchvar_exact, matchvar_fuzzy, caliper_fuzzy,
          covariates_ps, caliper_ps_sd, ratio, window_time,
          output_match = FALSE, trace = FALSE, cv_for_lambda,
          seed = seed + i)
        
        if (i %% max(1, floor(n_pair / 500)) == 0 || i == n_pair)
          cat(sprintf("[%s] %d/%d (%.0f%%) %s -> %s\n", 
                      Sys.time(), i, n_pair, round(i/n_pair*100), first_Ds[i], second_Ds[i]),
              file = "logs/progress.txt", append = TRUE)
        
        tibble(first_D = first_Ds[i], second_D = second_Ds[i],
               n_exposure = fit$n_exposure, n_unmatched = fit$n_unmatched,
               n_control = fit$n_control, n_event = fit$n_event) %>%
          cbind(fit$res_coef) %>%
          { cbind(coef_type = rownames(.), ., error = NA_character_) }
      }, error = function(e) {
        if (stop_on_error) stop(e)
        dA <- tryCatch(as.Date(data[[A_dates[i]]]), error = function(.) NULL)
        dB <- tryCatch(as.Date(data[[B_dates[i]]]), error = function(.) NULL)
        valid <- if (!is.null(dA) && !is.null(dB)) {
          idx <- as.integer(difftime(dB, dA, units = "days")) > window_time
          idx[is.na(idx)] <- FALSE; idx
        } else rep(FALSE, nrow(data))
        tibble(
          coef_type = NA_character_, first_D = first_Ds[i],
          second_D = second_Ds[i],
          n_exposure = sum(data[[A_events[i]]][valid], na.rm = TRUE),
          n_unmatched = NA_integer_, n_control = NA_integer_,
          n_event = sum(data[[B_events[i]]][valid], na.rm = TRUE) * (ratio + 1),
          coef = NA_real_, `exp(coef)` = NA_real_, `se(coef)` = NA_real_,
          z = NA_real_, `Pr(>|z|)` = NA_real_, error = conditionMessage(e))
      })
    }, cl = ncore)
    results <- rbindlist(results[!vapply(results, is.null, logical(1L))])

  } else {
    # Windows: socket 模式，clusterExport 导出 data（仅一次）
    ncore <- min(ncore, parallel::detectCores() - 1L, n_pair)
    cl <- snow::makeCluster(ncore, type = "SOCK")
    doSNOW::registerDoSNOW(cl)

    # 导出函数和数据（data 通过 clusterExport 而非词法作用域，避免重复序列化）
    snow::clusterExport(cl, c("getAB", "run_eds_sampling",
                              "calculate_propensity_score",
                              "build_matchvar_categories",
                              "resolve_disease_pairs",
                              "cut_cat", "extractAICBIC.glmnet",
                              "select_lasso_parameters", "run_eds_sampling_rcpp",
                              "extract_coxph_coef",
                              "SCRIPT_DIR"), envir = environment())
    snow::clusterEvalQ(cl, {
      library(tidyverse)
      library(data.table)
      library(survival)
      library(glmnet)
      library(Rcpp)
      Rcpp::sourceCpp(file.path(SCRIPT_DIR, "eds_rcpp.cpp"))
    })

    pb <- txtProgressBar(max = n_pair, style = 3)
    progress <- function(nn) setTxtProgressBar(pb, nn)
    opts <- list(progress = progress)

    results <- foreach(i = seq_len(n_pair),
                       .combine = function(x, y) rbindlist(list(x, y)),
                       .options.snow = opts) %dopar% {
      tryCatch({
        fit <- getAB(data, id, As[i], Bs[i],
          adjusted_D,
          start_date, matchvar_exact, matchvar_fuzzy, caliper_fuzzy,
          covariates_ps, caliper_ps_sd, ratio, window_time,
          output_match = FALSE, trace = FALSE, cv_for_lambda,
          seed = seed + i)
        tibble(first_D = first_Ds[i], second_D = second_Ds[i],
               n_exposure = fit$n_exposure, n_unmatched = fit$n_unmatched,
               n_control = fit$n_control, n_event = fit$n_event) %>%
          cbind(fit$res_coef) %>%
          { cbind(coef_type = rownames(.), ., error = NA_character_) }
      }, error = function(e) {
        if (stop_on_error) stop(e)
        dA <- tryCatch(as.Date(data[[A_dates[i]]]), error = function(.) NULL)
        dB <- tryCatch(as.Date(data[[B_dates[i]]]), error = function(.) NULL)
        valid <- if (!is.null(dA) && !is.null(dB)) {
          idx <- as.integer(difftime(dB, dA, units = "days")) > window_time
          idx[is.na(idx)] <- FALSE; idx
        } else rep(FALSE, nrow(data))
        tibble(
          coef_type = NA_character_, first_D = first_Ds[i],
          second_D = second_Ds[i],
          n_exposure = sum(data[[A_events[i]]][valid], na.rm = TRUE),
          n_unmatched = NA_integer_, n_control = NA_integer_,
          n_event = sum(data[[B_events[i]]][valid], na.rm = TRUE) * (ratio + 1),
          coef = NA_real_, `exp(coef)` = NA_real_, `se(coef)` = NA_real_,
          z = NA_real_, `Pr(>|z|)` = NA_real_, error = conditionMessage(e))
      })
    }
    close(pb)
    snow::stopCluster(cl)
  }

  rownames(results) <- NULL
  as.data.frame(results)
}






#' Estimate the effect of A on B given that C has already occurred
#'
#' Uses exposure density sampling (EDS) to estimate the causal effect of
#' disease A on disease B, conditional on disease C having occurred
#' before or at the same time as A. The exposure start time is redefined
#' as \code{max(A_time, C_time)}. C is excluded from \code{adjusted_D}.
#'
#' @param data A data frame containing patient-level data.
#' @param id Character. The column name of the unique patient identifier.
#' @param A Character. The name of the first disease (exposure).
#' @param B Character. The name of the second disease (outcome).
#' @param C Character. The name of the conditional disease. Only patients
#'   with \code{C_event == 1} and \code{C_time < B_time} are included.
#' @param start_date Character. Column name of the study start date (Date).
#' @param matchvar_exact Character vector. Exact matching variables.
#'   Use \code{NULL} to skip.
#' @param matchvar_fuzzy Character vector. Caliper matching variables
#'   (numeric only). Use \code{NULL} to skip.
#' @param caliper_fuzzy Numeric vector. Caliper bandwidths for
#'   \code{matchvar_fuzzy}. Must have the same length.
#' @param covariates_ps Character vector. Covariates for the PS model.
#'   Use \code{NULL} to skip.
#' @param caliper_ps_sd Numeric scalar. PS caliper in SD units.
#'   Default: 0.2.
#' @param ratio Integer. Controls per exposed case. Default: 1.
#' @param window_time Numeric. Minimum interval (days) between A and B.
#'   Default: 0.
#' @param adjusted_D Character vector. Disease names as candidate
#'   confounders (C is automatically excluded).
#' @param output_match Logical. If \code{TRUE}, include matched dataset
#'   in output. Default: \code{FALSE}.
#' @param trace Logical. Print EDS progress. Default: \code{FALSE}.
#' @param cv_for_lambda Logical. Compute CV-based lambda variants.
#'   Default: \code{FALSE}.
#' @param seed Integer. Random seed. Default: 1.
#' @return A list with the same structure as \code{getAB()}.
getAB_Ctrue <- function(data,id,A,B,C,adjusted_D,
                        start_date,matchvar_exact = NULL,matchvar_fuzzy = NULL,
                        caliper_fuzzy = NULL,covariates_ps = NULL,caliper_ps_sd = 0.2,ratio = 1,
                        window_time = 0,output_match = FALSE,
                        trace = FALSE, cv_for_lambda = F,seed = 1){
  # 参数标准化：NULL → 零长度向量
  matchvar_exact  <- if (is.null(matchvar_exact))  character(0) else matchvar_exact
  matchvar_fuzzy  <- if (is.null(matchvar_fuzzy))  character(0) else matchvar_fuzzy
  caliper_fuzzy   <- if (is.null(caliper_fuzzy))   numeric(0)   else caliper_fuzzy
  covariates_ps   <- if (is.null(covariates_ps))   character(0) else covariates_ps
  data <- as.data.frame(data)
  if(any(sapply(c(id,A,B,C,start_date,ratio,caliper_ps_sd),length) != 1))
    stop("The length of id,A,B,C,start_date,ratio,caliper_ps_sd must be 1")
  if(length(matchvar_fuzzy) != length(caliper_fuzzy)) stop("matchvar_fuzzy and caliper_fuzzy have different length")
  
  A_date <- paste0(A,"_date")
  B_date <- paste0(B,"_date")
  C_date <- paste0(C,"_date")
  A_event <- paste0(A,"_event")
  B_event <- paste0(B,"_event")
  C_event <- paste0(C,"_event")
  adjusted_D <- setdiff(adjusted_D, c(A, B, C))
  adjusted_D_date <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_date") else character(0)
  adjusted_D_event <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_event") else character(0)
  
  if (!start_date %in% names(data)) stop("start_date column not found in data")
  if (class(data[[start_date]]) != "Date") stop("start_date must be Date type")
  if (any(data[[start_date]] > data[[A_date]] | data[[start_date]] > data[[B_date]], na.rm = TRUE))
    stop("start_date must be <= A_date and B_date for all rows")
  
  keepvars <- c(id,A_date,A_event,B_date,B_event,C_date,C_event,start_date,matchvar_exact,
                matchvar_fuzzy,covariates_ps,adjusted_D_date,adjusted_D_event)
  if(!all(keepvars %in% names(data))) 
    stop("There were variables not being included in the data")
  if(any(sapply(data[,c(A_date,B_date,C_date,start_date,adjusted_D_date)], class) != "Date")) stop("There are date variables not being Date type")
  if(class(data[,A_event]) == "logical") data[,A_event] <- as.integer(data[,A_event])
  if(class(data[,B_event]) == "logical") data[,B_event] <- as.integer(data[,B_event])
  if(class(data[,C_event]) == "logical") data[,C_event] <- as.integer(data[,C_event])
  if(any(duplicated(c(A,B,C)))) stop("There is duplicated variable for A, B and C")
  if(!all(sapply(data[,matchvar_fuzzy,drop = FALSE], class) %in% c("integer","numeric"))) stop("matchvar_fuzzy must be a numeric or integer variable")
  
  origin_data <- na.omit(data[,c(id,A_date,A_event,B_date,B_event,C_date, C_event,start_date,matchvar_exact,
                                 matchvar_fuzzy,covariates_ps,adjusted_D_date,adjusted_D_event)])  
  for (i in adjusted_D) {
    origin_data[,paste0(i,"_time")] <- difftime(origin_data[,paste0(i,"_date")],origin_data[,start_date],units = "days") %>% as.integer
  }
  
  data <- origin_data[,c(id,A_date,A_event,B_date,B_event,C_date,C_event,start_date,matchvar_exact,matchvar_fuzzy,covariates_ps)]  
  data[,"A_time"] <- difftime(data[,A_date],data[,start_date],units = "days") %>% as.integer
  data[,"B_time"] <- difftime(data[,B_date],data[,start_date],units = "days") %>% as.integer
  data[,"C_time"] <- difftime(data[,C_date],data[,start_date],units = "days") %>% as.integer
  data[,"A_event"] <- data[,A_event]
  data[,"B_event"] <- data[,B_event]
  data[,"C_event"] <- data[,C_event]
  
  # 提取C_event == 1，并重定义A_date
  data <- as.data.frame(data) %>% 
    filter(C_event == 1 & C_time < B_time) %>%
    mutate(
      A_time_o = A_time,
      A_time = fifelse(C_time > A_time, C_time, A_time)
    )

  # 过滤负时间值（A_date或B_date在start_date之前的情况）
  data <- subset(data, A_time > 0 & B_time > 0)

  # 根据matchvar_fuzzy划分类别并计算倾向性评分（支持空参数）
  # 注：PS模型使用原始A_time_o估计发生A的概率，EDS使用重定义的A_time
  if (length(matchvar_fuzzy) > 0) {
    data <- build_matchvar_categories(data, matchvar_fuzzy, caliper_fuzzy)
  }
  # 占位符：确保 matchvar_names 永不为空
  if (length(matchvar_exact) == 0) {
    data[["_placeholder_exact"]] <- 1L
  }
  if (length(matchvar_fuzzy) == 0) {
    data[["_placeholder_fuzzy_cat"]] <- 1L
  }
  matchvar_fuzzy_cat <- if (length(matchvar_fuzzy) > 0) paste0(matchvar_fuzzy, "_cat") else character(0)
  if (length(covariates_ps) > 0) {
    data <- calculate_propensity_score(data, covariates_ps, matchvar_exact,
                                        matchvar_fuzzy_cat, caliper_ps_sd,
                                        ps_time_var = "A_time_o",
                                        ps_filter_expr = "B_time >= A_time")
  } else {
    data[["ps"]] <- 1
    data[["ps_cat"]] <- 1L
  }


  ######## EDS抽样
  matchvar_names <- c(
    if (length(matchvar_exact) > 0) matchvar_exact else "_placeholder_exact",
    if (length(matchvar_fuzzy) > 0) paste0(matchvar_fuzzy, "_cat") else "_placeholder_fuzzy_cat",
    "ps_cat"
  )
  eds_result <- run_eds_sampling_rcpp(
    data = data, matchvar_names = matchvar_names,
    ratio = ratio, window_time = window_time, seed = seed,
    exposure_filter_expr = "C_time <= A_time & A_time + window_time < B_time & A_event == 1",
    control_filter_expr = "C_time <= start_time_i & A_time > start_time_i + window_time & B_time > start_time_i + window_time",
    id = id, include_C_time = TRUE, trace = trace
  )
  data_matched_list <- eds_result$data_matched_list
  exposure_dt <- eds_result$exposure_dt
  
  # 剔除没有匹配到对照的list
  NA_index <- sapply(data_matched_list, function(x) anyNA(x[[id]])) %>% which
  if(length(NA_index) >= 1) data_matched_list <- data_matched_list[-NA_index]
  
  
  adjusted_D_time <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_time") else character(0)
  data_matched <- rbindlist(data_matched_list) %>%
    left_join(data[,c(id,matchvar_names,matchvar_fuzzy,"ps","A_time","B_time","C_time","A_event","B_event","C_event")],by = id) %>%
    left_join(origin_data[,c(id,adjusted_D_time,adjusted_D_event),drop = FALSE],by = id) %>%
    group_by(!!!syms(matchvar_names)) %>%
    mutate(strata = cur_group_id()) %>%
    ungroup()
  
  # 处理对照组
  # if B_time ≤ A_time, start_time = A_time_e，end_time = B_time，status = B_event
  # if B_time > A_time, start_time =A_time_e，end_time = A_time，status = 0
  data_matched <- data_matched %>% mutate(
    end_time = fifelse(
      exposure == 0 & A_time < B_time, A_time,
      B_time),
    status = fifelse(
      exposure == 0 & A_time < B_time, 0L,
      B_event),
  )
  
  # 去掉strata内部的重复个体，保留start_time最早的
  data_matched <- data_matched %>% 
    group_by(strata,!!sym(id),exposure) %>%
    slice_min(start_time, n = 1, with_ties = FALSE) %>%  # 保留 start_time 最小的行
    ungroup()
  
  # 去掉没有事件发生的strata
  data_matched_model <- data_matched %>%
    group_by(strata) %>%
    filter(any(status == 1)) %>%  # 保留至少有一个事件的 strata
    ungroup()
  
  # 添加其它condition的时依协变量
  data_matched_model$unique_row <- 1:nrow(data_matched_model)
  tdc_part <- if (length(adjusted_D) > 0) {
    paste0(",", paste0(paste0("`", adjusted_D, "` = tdc(`", adjusted_D, "_time`)"), collapse = ","))
  } else {
    ""
  }
  tdc_expr <- paste0(
    "tmerge(data_matched_model, data_matched_model, id=unique_row,tstart = start_time,tstop = end_time,status = event(end_time,status)",
    tdc_part, ")"
  )
  td_matched_data <- eval(parse(text = tdc_expr)) 
  
  n_unmatched <- nrow(exposure_dt) - sum(data_matched$exposure == 1)
  if(n_unmatched > 0) warning(sprintf("There are %.0f participants having not been matched for %.0f participants",n_unmatched,nrow(exposure_dt)))
  n_exposure <- nrow(exposure_dt)
  n_control <- sum(data_matched$exposure==0)
  n_event <- sum(data_matched$status==1)
  #### fit model: cox + lasso

  # not adjusting for other conditions
  fit_unadjusted_D <- coxph(Surv(tstart,tstop,status) ~ exposure + strata(strata),data = td_matched_data)
  coef_unadjusted_D <- extract_coxph_coef(fit_unadjusted_D, "exposure")[1L, ]
  
  # adjusting for other conditions
  n_total <- nrow(td_matched_data)
  params <- select_lasso_parameters(n_total)
  if (length(adjusted_D) > 0) {
  nlambda_val <- params$nlambda_val
  lambda_min_ratio_val <- params$lambda_min_ratio_val
  nfolds_val <- params$nfolds_val
  control_list <- params$control_list
  x0 <- td_matched_data[,c("exposure",adjusted_D)] %>% as.matrix()
  y0 <- stratifySurv(Surv(td_matched_data$tstart,td_matched_data$tstop,td_matched_data$status),strata = td_matched_data$strata)
  args_list <- list(
    x = x0,y = y0,family = "cox",penalty.factor = c(0,rep(1, length(adjusted_D))),
    relax = T,nlambda = nlambda_val,cox.ties = "breslow",
    lambda.min.ratio = lambda_min_ratio_val,control = control_list,
    standardize = FALSE
  )
  fit_lasso <- do.call(glmnet,args_list)
  
  lambda_aic_bic <- extractAICBIC.glmnet(fit_lasso,n_event)
  
  coef_lasso <- coef(fit_lasso,s = lambda_aic_bic$lambda.AIC)
  vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
  fit_cox <- coxph(cox_formula,data = td_matched_data)
  coef_AIC_lambda <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  
  coef_lasso <- coef(fit_lasso,s = lambda_aic_bic$lambda.BIC)
  vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
  fit_cox <- coxph(cox_formula,data = td_matched_data)
  coef_BIC_lambda <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  
  if(cv_for_lambda){
    fit_lasso_cv <- do.call(cv.glmnet,modifyList(args_list,val = list(gamma = 0,nfolds = nfolds_val)))
    coef_lasso <- coef(fit_lasso_cv,s = fit_lasso_cv$lambda.min)
    vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
    fit_cox <- coxph(cox_formula,data = td_matched_data)
    coef_cv_lambda.min <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
    
    coef_lasso <- coef(fit_lasso_cv,s = fit_lasso_cv$lambda.1se)
    vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
    fit_cox <- coxph(cox_formula,data = td_matched_data)
    coef_cv_lambda.1se <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  } else{
    coef_cv_lambda.1se <- coef_cv_lambda.min <- NA
  }
    
  } else {
    coef_AIC_lambda <- coef_BIC_lambda <- coef_cv_lambda.min <- coef_cv_lambda.1se <- NA
  }
  res_coef <- rbind(coef_unadjusted_D,coef_AIC_lambda,coef_BIC_lambda,coef_cv_lambda.min,coef_cv_lambda.1se)
  
  if(output_match) {
    return(list(res_coef = res_coef,n_exposure = n_exposure,n_control = n_control,
                 n_unmatched = n_unmatched,n_event = n_event,data_matched = data_matched))
  } else {
    return(list(res_coef = res_coef,n_exposure = n_exposure,n_control = n_control,
                 n_unmatched = n_unmatched,n_event = n_event))
  }
}








#' Batch estimation for C-true conditional effects
#'
#' Runs \code{getAB_Ctrue()} for multiple disease pairs in parallel,
#' estimating the effect of A on B given that C has already occurred.
#'
#' @param data A data frame containing patient-level data.
#' @param id Character. The column name of the unique patient identifier.
#' @param D_all Character vector. All disease names for permutation.
#' @param As Character vector. First disease of each pair.
#' @param Bs Character vector. Second disease of each pair.
#' @param C Character. The conditional disease name.
#' @param start_date Character. Column name of the study start date (Date).
#' @param matchvar_exact Character vector. Exact matching variables.
#'   Use \code{NULL} to skip.
#' @param matchvar_fuzzy Character vector. Caliper matching variables
#'   (numeric only). Use \code{NULL} to skip.
#' @param caliper_fuzzy Numeric vector. Caliper bandwidths for
#'   \code{matchvar_fuzzy}. Must have the same length.
#' @param covariates_ps Character vector. Covariates for the PS model.
#'   Use \code{NULL} to skip.
#' @param caliper_ps_sd Numeric scalar. PS caliper in SD units.
#'   Default: 0.2.
#' @param ratio Integer. Controls per exposed case. Default: 1.
#' @param window_time Numeric. Minimum interval (days) between A and B.
#'   Default: 0.
#' @param adjusted_D Character vector. Disease names as candidate
#'   confounders (C is automatically excluded).
#' @param ncore Integer. Number of CPU cores for parallel execution.
#'   Default: 1.
#' @param cv_for_lambda Logical. Compute CV-based lambda variants.
#'   Default: \code{FALSE}.
#' @param seed Integer. Base random seed (\code{seed + i} per pair).
#'   Default: 1.
#' @param parallel_backend Character. Parallel backend:
#'   \code{"auto"} (default), \code{"fork"}, or \code{"socket"}.
#' @param stop_on_error Logical. If \code{FALSE} (default), return
#'   \code{NA} rows for failing pairs.
#' @param ... Additional arguments (currently unused).
#' @return A data frame with the same structure as \code{getAB_batch()}.
getAB_Ctrue_batch <- function(data,id,D_all,As,Bs,C,start_date,matchvar_exact = NULL,
                               matchvar_fuzzy = NULL,caliper_fuzzy = NULL,covariates_ps = NULL,caliper_ps_sd = 0.2,
                               ratio = 1,window_time = 0,adjusted_D,ncore = 1,
                               cv_for_lambda = F,seed = 1,
                               parallel_backend = c("auto", "fork", "socket"),
                               stop_on_error = FALSE,...){
  # 参数标准化
  matchvar_exact  <- if (is.null(matchvar_exact))  character(0) else matchvar_exact
  matchvar_fuzzy  <- if (is.null(matchvar_fuzzy))  character(0) else matchvar_fuzzy
  caliper_fuzzy   <- if (is.null(caliper_fuzzy))   numeric(0)   else caliper_fuzzy
  covariates_ps   <- if (is.null(covariates_ps))   character(0) else covariates_ps
  if(length(matchvar_fuzzy) != length(caliper_fuzzy)) stop("matchvar_fuzzy and caliper_fuzzy have different length")

  # 解析疾病对
  pairs <- resolve_disease_pairs(data, D_all, As, Bs, adjusted_D, C = C)
  As <- pairs$As; Bs <- pairs$Bs; adjusted_D <- pairs$adjusted_D
  n_pair <- pairs$n_pair
  A_dates <- pairs$A_dates; B_dates <- pairs$B_dates
  A_events <- pairs$A_events; B_events <- pairs$B_events
  first_Ds <- pairs$first_Ds; second_Ds <- pairs$second_Ds

  # 选择并行后端
  parallel_backend <- match.arg(parallel_backend)
  if (parallel_backend == "auto") {
    parallel_backend <- if (.Platform$OS.type == "windows") "socket" else "fork"
  }

  if (parallel_backend == "fork") {
    # Linux/macOS: fork 模式，共享内存，无需序列化 data
    ncore <- min(ncore, parallel::detectCores() - 1L, n_pair)
    results <- pbapply::pblapply(seq_len(n_pair), function(i) {
      tryCatch({
        fit <- getAB_Ctrue(data, id, As[i], Bs[i],
          C, adjusted_D,
          start_date, matchvar_exact, matchvar_fuzzy, caliper_fuzzy,
          covariates_ps, caliper_ps_sd, ratio, window_time,
          output_match = FALSE, trace = FALSE, cv_for_lambda,
          seed = seed + i)
        
        if (i %% max(1, floor(n_pair / 500)) == 0 || i == n_pair)
          cat(sprintf("[%s] %d/%d (%.0f%%) %s -> %s\n", 
                      Sys.time(), i, n_pair, round(i/n_pair*100), first_Ds[i], second_Ds[i]),
              file = "logs/progress.txt", append = TRUE)
        
        tibble(first_D = first_Ds[i], second_D = second_Ds[i],
               n_exposure = fit$n_exposure, n_unmatched = fit$n_unmatched,
               n_control = fit$n_control, n_event = fit$n_event) %>%
          cbind(fit$res_coef) %>%
          { cbind(coef_type = rownames(.), ., error = NA_character_) }
      }, error = function(e) {
        if (stop_on_error) stop(e)
        C_event_col <- paste0(C, "_event")
        C_date_col  <- paste0(C, "_date")
        has_C_event  <- data[[C_event_col]] == 1
        has_C_event[is.na(has_C_event)] <- FALSE
        dA <- tryCatch(as.Date(data[[A_dates[i]]]), error = function(.) NULL)
        dB <- tryCatch(as.Date(data[[B_dates[i]]]), error = function(.) NULL)
        dC <- tryCatch(as.Date(data[[C_date_col]]),  error = function(.) NULL)
        valid <- if (!is.null(dA) && !is.null(dB) && !is.null(dC)) {
          has_C_event & as.integer(difftime(dC, dA, units = "days")) <= 0 &
            as.integer(difftime(dB, dA, units = "days")) > window_time
        } else rep(FALSE, nrow(data))
        valid[is.na(valid)] <- FALSE
        tibble(
          coef_type = NA_character_, first_D = first_Ds[i],
          second_D = second_Ds[i],
          n_exposure = sum(data[[A_events[i]]][valid], na.rm = TRUE),
          n_unmatched = NA_integer_, n_control = NA_integer_,
          n_event = sum(data[[B_events[i]]][valid], na.rm = TRUE),
          coef = NA_real_, `exp(coef)` = NA_real_, `se(coef)` = NA_real_,
          z = NA_real_, `Pr(>|z|)` = NA_real_, error = conditionMessage(e))
      })
    }, cl = ncore)
    results <- rbindlist(results[!vapply(results, is.null, logical(1L))])

  } else {
    # Windows: socket 模式，clusterExport 导出 data（仅一次）
    ncore <- min(ncore, parallel::detectCores() - 1L, n_pair)
    cl <- snow::makeCluster(ncore, type = "SOCK")
    doSNOW::registerDoSNOW(cl)

    # 导出函数和数据（data 通过 clusterExport 而非词法作用域，避免重复序列化）
    snow::clusterExport(cl, c("getAB_Ctrue", "run_eds_sampling",
                              "calculate_propensity_score",
                              "build_matchvar_categories",
                              "resolve_disease_pairs",
                              "cut_cat", "extractAICBIC.glmnet",
                              "select_lasso_parameters", "run_eds_sampling_rcpp",
                              "extract_coxph_coef",
                              "SCRIPT_DIR"), envir = environment())
    snow::clusterEvalQ(cl, {
      library(tidyverse)
      library(data.table)
      library(survival)
      library(glmnet)
      library(Rcpp)
      Rcpp::sourceCpp(file.path(SCRIPT_DIR, "eds_rcpp.cpp"))
    })

    pb <- txtProgressBar(max = n_pair, style = 3)
    progress <- function(nn) setTxtProgressBar(pb, nn)
    opts <- list(progress = progress)

    results <- foreach(i = seq_len(n_pair),
                       .combine = function(x, y) rbindlist(list(x, y)),
                       .options.snow = opts) %dopar% {
      tryCatch({
        fit <- getAB_Ctrue(data, id, As[i], Bs[i],
          C, adjusted_D,
          start_date, matchvar_exact, matchvar_fuzzy, caliper_fuzzy,
          covariates_ps, caliper_ps_sd, ratio, window_time,
          output_match = FALSE, trace = FALSE, cv_for_lambda,
          seed = seed + i)
        tibble(first_D = first_Ds[i], second_D = second_Ds[i],
               n_exposure = fit$n_exposure, n_unmatched = fit$n_unmatched,
               n_control = fit$n_control, n_event = fit$n_event) %>%
          cbind(fit$res_coef) %>%
          { cbind(coef_type = rownames(.), ., error = NA_character_) }
      }, error = function(e) {
        if (stop_on_error) stop(e)
        C_event_col <- paste0(C, "_event")
        C_date_col  <- paste0(C, "_date")
        has_C_event  <- data[[C_event_col]] == 1
        has_C_event[is.na(has_C_event)] <- FALSE
        dA <- tryCatch(as.Date(data[[A_dates[i]]]), error = function(.) NULL)
        dB <- tryCatch(as.Date(data[[B_dates[i]]]), error = function(.) NULL)
        dC <- tryCatch(as.Date(data[[C_date_col]]),  error = function(.) NULL)
        valid <- if (!is.null(dA) && !is.null(dB) && !is.null(dC)) {
          has_C_event & as.integer(difftime(dC, dA, units = "days")) <= 0 &
            as.integer(difftime(dB, dA, units = "days")) > window_time
        } else rep(FALSE, nrow(data))
        valid[is.na(valid)] <- FALSE
        tibble(
          coef_type = NA_character_, first_D = first_Ds[i],
          second_D = second_Ds[i],
          n_exposure = sum(data[[A_events[i]]][valid], na.rm = TRUE),
          n_unmatched = NA_integer_, n_control = NA_integer_,
          n_event = sum(data[[B_events[i]]][valid], na.rm = TRUE),
          coef = NA_real_, `exp(coef)` = NA_real_, `se(coef)` = NA_real_,
          z = NA_real_, `Pr(>|z|)` = NA_real_, error = conditionMessage(e))
      })
    }
    close(pb)
    snow::stopCluster(cl)
  }

  rownames(results) <- NULL
  as.data.frame(results)
}





#' Estimate the effect of A on B given that C has not yet occurred
#'
#' Uses exposure density sampling (EDS) to estimate the causal effect of
#' disease A on disease B, conditional on disease C not having occurred
#' before A. Patients are right-censored at \code{C_time} if C occurs
#' before B. C is excluded from \code{adjusted_D}.
#'
#' @param data A data frame containing patient-level data.
#' @param id Character. The column name of the unique patient identifier.
#' @param A Character. The name of the first disease (exposure).
#' @param B Character. The name of the second disease (outcome).
#' @param C Character. The name of the conditional disease. Only patients
#'   with \code{C_time > A_time} (C not yet occurred) are included.
#' @param start_date Character. Column name of the study start date (Date).
#' @param matchvar_exact Character vector. Exact matching variables.
#'   Use \code{NULL} to skip.
#' @param matchvar_fuzzy Character vector. Caliper matching variables
#'   (numeric only). Use \code{NULL} to skip.
#' @param caliper_fuzzy Numeric vector. Caliper bandwidths for
#'   \code{matchvar_fuzzy}. Must have the same length.
#' @param covariates_ps Character vector. Covariates for the PS model.
#'   Use \code{NULL} to skip.
#' @param caliper_ps_sd Numeric scalar. PS caliper in SD units.
#'   Default: 0.2.
#' @param ratio Integer. Controls per exposed case. Default: 1.
#' @param window_time Numeric. Minimum interval (days) between A and B.
#'   Default: 0.
#' @param adjusted_D Character vector. Disease names as candidate
#'   confounders (C is automatically excluded).
#' @param output_match Logical. If \code{TRUE}, include matched dataset
#'   in output. Default: \code{FALSE}.
#' @param trace Logical. Print EDS progress. Default: \code{FALSE}.
#' @param cv_for_lambda Logical. Compute CV-based lambda variants.
#'   Default: \code{FALSE}.
#' @param seed Integer. Random seed. Default: 1.
#' @return A list with the same structure as \code{getAB()}.
getAB_Cfalse <- function(data,id,A,B,C,adjusted_D,
                        start_date,matchvar_exact = NULL,matchvar_fuzzy = NULL,
                        caliper_fuzzy = NULL,covariates_ps = NULL,caliper_ps_sd = 0.2,ratio = 1,
                        window_time = 0,output_match = FALSE,
                        trace = FALSE, cv_for_lambda = F,seed = 1){
  # 参数标准化：NULL → 零长度向量
  matchvar_exact  <- if (is.null(matchvar_exact))  character(0) else matchvar_exact
  matchvar_fuzzy  <- if (is.null(matchvar_fuzzy))  character(0) else matchvar_fuzzy
  caliper_fuzzy   <- if (is.null(caliper_fuzzy))   numeric(0)   else caliper_fuzzy
  covariates_ps   <- if (is.null(covariates_ps))   character(0) else covariates_ps
  data <- as.data.frame(data)
  if(any(sapply(c(id,A,B,C,start_date,ratio,caliper_ps_sd),length) != 1))
    stop("The length of id,A,B,C,start_date,ratio,caliper_ps_sd must be 1")
  if(length(matchvar_fuzzy) != length(caliper_fuzzy)) stop("matchvar_fuzzy and caliper_fuzzy have different length")
  
  A_date <- paste0(A,"_date")
  B_date <- paste0(B,"_date")
  C_date <- paste0(C,"_date")
  A_event <- paste0(A,"_event")
  B_event <- paste0(B,"_event")
  C_event <- paste0(C,"_event")
  adjusted_D <- setdiff(adjusted_D, c(A, B, C))
  adjusted_D_date <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_date") else character(0)
  adjusted_D_event <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_event") else character(0)
  
  if (!start_date %in% names(data)) stop("start_date column not found in data")
  if (class(data[[start_date]]) != "Date") stop("start_date must be Date type")
  if (any(data[[start_date]] > data[[A_date]] | data[[start_date]] > data[[B_date]], na.rm = TRUE))
    stop("start_date must be <= A_date and B_date for all rows")
  
  keepvars <- c(id,A_date,A_event,B_date,B_event,C_date,C_event,start_date,matchvar_exact,
                matchvar_fuzzy,covariates_ps,adjusted_D_date,adjusted_D_event)
  if(!all(keepvars %in% names(data))) 
    stop("There were variables not being included in the data")
  if(any(sapply(data[,c(A_date,B_date,C_date,start_date,adjusted_D_date)], class) != "Date")) stop("There are date variables not being Date type")
  if(class(data[,A_event]) == "logical") data[,A_event] <- as.integer(data[,A_event])
  if(class(data[,B_event]) == "logical") data[,B_event] <- as.integer(data[,B_event])
  if(class(data[,C_event]) == "logical") data[,C_event] <- as.integer(data[,C_event])
  if(any(duplicated(c(A,B,C)))) stop("There is duplicated variable for A, B and C")
  if(!all(sapply(data[,matchvar_fuzzy,drop = FALSE], class) %in% c("integer","numeric"))) stop("matchvar_fuzzy must be a numeric or integer variable")
  
  origin_data <- na.omit(data[,c(id,A_date,A_event,B_date,B_event,C_date, C_event,start_date,matchvar_exact,
                                 matchvar_fuzzy,covariates_ps,adjusted_D_date,adjusted_D_event)])  
  for (i in adjusted_D) {
    origin_data[,paste0(i,"_time")] <- difftime(origin_data[,paste0(i,"_date")],origin_data[,start_date],units = "days") %>% as.integer
  }
  
  data <- origin_data[,c(id,A_date,A_event,B_date,B_event,C_date,C_event,start_date,matchvar_exact,matchvar_fuzzy,covariates_ps)]  
  data[,"A_time"] <- difftime(data[,A_date],data[,start_date],units = "days") %>% as.integer
  data[,"B_time"] <- difftime(data[,B_date],data[,start_date],units = "days") %>% as.integer
  data[,"C_time"] <- difftime(data[,C_date],data[,start_date],units = "days") %>% as.integer
  data[,"A_event"] <- data[,A_event]
  data[,"B_event"] <- data[,B_event]
  data[,"C_event"] <- data[,C_event]

  # 过滤负时间值（A_date或B_date在start_date之前的情况）
  data <- subset(data, A_time > 0 & B_time > 0)

  # 根据matchvar_fuzzy划分类别并计算倾向性评分（支持空参数）
  if (length(matchvar_fuzzy) > 0) {
    data <- build_matchvar_categories(data, matchvar_fuzzy, caliper_fuzzy)
  }
  # 占位符：确保 matchvar_names 永不为空
  if (length(matchvar_exact) == 0) {
    data[["_placeholder_exact"]] <- 1L
  }
  if (length(matchvar_fuzzy) == 0) {
    data[["_placeholder_fuzzy_cat"]] <- 1L
  }
  matchvar_fuzzy_cat <- if (length(matchvar_fuzzy) > 0) paste0(matchvar_fuzzy, "_cat") else character(0)
  if (length(covariates_ps) > 0) {
    data <- calculate_propensity_score(data, covariates_ps, matchvar_exact,
                                        matchvar_fuzzy_cat, caliper_ps_sd,
                                        ps_time_var = "A_time",
                                        ps_filter_expr = "B_time >= A_time & C_time > A_time")
  } else {
    data[["ps"]] <- 1
    data[["ps_cat"]] <- 1L
  }


  ######## EDS抽样
  matchvar_names <- c(
    if (length(matchvar_exact) > 0) matchvar_exact else "_placeholder_exact",
    if (length(matchvar_fuzzy) > 0) paste0(matchvar_fuzzy, "_cat") else "_placeholder_fuzzy_cat",
    "ps_cat"
  )
  eds_result <- run_eds_sampling_rcpp(
    data = data, matchvar_names = matchvar_names,
    ratio = ratio, window_time = window_time, seed = seed,
    exposure_filter_expr = "C_time > A_time & A_time + window_time < B_time & A_event == 1",
    control_filter_expr = "C_time > start_time_i & A_time > start_time_i + window_time & B_time > start_time_i + window_time",
    id = id, include_C_time = TRUE, trace = trace
  )
  data_matched_list <- eds_result$data_matched_list
  exposure_dt <- eds_result$exposure_dt
  
  # 剔除没有匹配到对照的list
  NA_index <- sapply(data_matched_list, function(x) anyNA(x[[id]])) %>% which
  if(length(NA_index) >= 1) data_matched_list <- data_matched_list[-NA_index]
  
  
  adjusted_D_time <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_time") else character(0)
  data_matched <- rbindlist(data_matched_list) %>%
    left_join(data[,c(id,matchvar_names,matchvar_fuzzy,"ps","A_time","B_time","C_time","A_event","B_event","C_event")],by = id) %>%
    left_join(origin_data[,c(id,adjusted_D_time,adjusted_D_event),drop = FALSE],by = id) %>%
    group_by(!!!syms(matchvar_names)) %>%
    mutate(strata = cur_group_id()) %>%
    ungroup()
  
  # 处理对照与暴露组
  # 暴露组：if C_time < B_time, start_time = A_time_e, end_time = C_time，status = 0
  # If C_time ≥ B_time, start_time = A_time_e, end_time = B_time，status = B_event
  # 对照组:  if C_time ≤ A_time & C_time < B_time, start_time = A_time_e，end_time = C_time，status = 0
  # if A_time < C_time & A_time < B_time, start_time = A_time_e，end_time = A_time，status = 0
  # if B_time ≤ C_time & B_time ≤  A_time, start_time = A_time_e，end_time = B_time，status = B_event
  
  data_matched <- data_matched %>% mutate(
    end_time = fifelse(
      (exposure == 1 & C_time < B_time) | (exposure == 0 & C_time <= A_time & C_time < B_time), C_time,
      fifelse(exposure == 0 & A_time < C_time & A_time < B_time, A_time, B_time)),
    status = fifelse(
      (exposure == 1 & C_time < B_time) | (exposure == 0 & C_time <= A_time & C_time < B_time) | (exposure == 0 & A_time < C_time & A_time < B_time),
      0L, B_event)
  )
  
  # 去掉strata内部的重复个体，保留start_time最早的
  data_matched <- data_matched %>% 
    group_by(strata,!!sym(id),exposure) %>%
    slice_min(start_time, n = 1, with_ties = FALSE) %>%  # 保留 start_time 最小的行
    ungroup()
  
  # 去掉没有事件发生的strata
  data_matched_model <- data_matched %>%
    group_by(strata) %>%
    filter(any(status == 1)) %>%  # 保留至少有一个事件的 strata
    ungroup()
  
  # 添加其它condition的时依协变量
  data_matched_model$unique_row <- 1:nrow(data_matched_model)
  tdc_part <- if (length(adjusted_D) > 0) {
    paste0(",", paste0(paste0("`", adjusted_D, "` = tdc(`", adjusted_D, "_time`)"), collapse = ","))
  } else {
    ""
  }
  tdc_expr <- paste0(
    "tmerge(data_matched_model, data_matched_model, id=unique_row,tstart = start_time,tstop = end_time,status = event(end_time,status)",
    tdc_part, ")"
  )
  td_matched_data <- eval(parse(text = tdc_expr)) 
  
  n_unmatched <- nrow(exposure_dt) - sum(data_matched$exposure == 1)
  if(n_unmatched > 0) warning(sprintf("There are %.0f participants having not been matched for %.0f participants",n_unmatched,nrow(exposure_dt)))
  n_exposure <- nrow(exposure_dt)
  n_control <- sum(data_matched$exposure==0)
  n_event <- sum(data_matched$status==1)
  #### fit model: cox + lasso

  # not adjusting for other conditions
  fit_unadjusted_D <- coxph(Surv(tstart,tstop,status) ~ exposure + strata(strata),data = td_matched_data)
  coef_unadjusted_D <- extract_coxph_coef(fit_unadjusted_D, "exposure")[1L, ]
  if (length(adjusted_D) > 0) {
  
  # adjusting for other conditions
  n_total <- nrow(td_matched_data)
  params <- select_lasso_parameters(n_total)
  nlambda_val <- params$nlambda_val
  lambda_min_ratio_val <- params$lambda_min_ratio_val
  nfolds_val <- params$nfolds_val
  control_list <- params$control_list
  x0 <- td_matched_data[,c("exposure",adjusted_D)] %>% as.matrix()
  y0 <- stratifySurv(Surv(td_matched_data$tstart,td_matched_data$tstop,td_matched_data$status),strata = td_matched_data$strata)
  args_list <- list(
    x = x0,y = y0,family = "cox",penalty.factor = c(0,rep(1, length(adjusted_D))),
    relax = T,nlambda = nlambda_val,cox.ties = "breslow",
    lambda.min.ratio = lambda_min_ratio_val,control = control_list,
    standardize = FALSE
  )
  fit_lasso <- do.call(glmnet,args_list)
  
  
  lambda_aic_bic <- extractAICBIC.glmnet(fit_lasso,n_event)
  
  coef_lasso <- coef(fit_lasso,s = lambda_aic_bic$lambda.AIC)
  vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
  fit_cox <- coxph(cox_formula,data = td_matched_data)
  coef_AIC_lambda <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  
  coef_lasso <- coef(fit_lasso,s = lambda_aic_bic$lambda.BIC)
  vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
  fit_cox <- coxph(cox_formula,data = td_matched_data)
  coef_BIC_lambda <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  
  if(cv_for_lambda){
    # fit_lasso_cv <- cv.glmnet(x = x0,y = y0,family = "cox",penalty.factor = c(0,rep(1, length(adjusted_D))),nfolds = 5,relax = T,gamma = 0,nlambda = 30)
    fit_lasso_cv <- do.call(cv.glmnet,modifyList(args_list,val = list(gamma = 0,nfolds = nfolds_val)))
    coef_lasso <- coef(fit_lasso_cv,s = fit_lasso_cv$lambda.min)
    vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
    fit_cox <- coxph(cox_formula,data = td_matched_data)
    coef_cv_lambda.min <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
    
    coef_lasso <- coef(fit_lasso_cv,s = fit_lasso_cv$lambda.1se)
    vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
    fit_cox <- coxph(cox_formula,data = td_matched_data)
    coef_cv_lambda.1se <- extract_coxph_coef(fit_cox, "exposure")[1L, ]
  } else{
    coef_cv_lambda.1se <- coef_cv_lambda.min <- NA
  }
  
  } else {
    coef_AIC_lambda <- coef_BIC_lambda <- coef_cv_lambda.min <- coef_cv_lambda.1se <- NA
  }
  res_coef <- rbind(coef_unadjusted_D,coef_AIC_lambda,coef_BIC_lambda,coef_cv_lambda.min,coef_cv_lambda.1se)
  
  if(output_match) {
    return(list(res_coef = res_coef,n_exposure = n_exposure,n_control = n_control,
                 n_unmatched = n_unmatched,n_event = n_event,data_matched = data_matched))
  } else {
    return(list(res_coef = res_coef,n_exposure = n_exposure,n_control = n_control,
                 n_unmatched = n_unmatched,n_event = n_event))
  }
}



#' Batch estimation for C-false conditional effects
#'
#' Runs \code{getAB_Cfalse()} for multiple disease pairs in parallel,
#' estimating the effect of A on B given that C has not yet occurred.
#'
#' @param data A data frame containing patient-level data.
#' @param id Character. The column name of the unique patient identifier.
#' @param D_all Character vector. All disease names for permutation.
#' @param As Character vector. First disease of each pair.
#' @param Bs Character vector. Second disease of each pair.
#' @param C Character. The conditional disease name.
#' @param start_date Character. Column name of the study start date (Date).
#' @param matchvar_exact Character vector. Exact matching variables.
#'   Use \code{NULL} to skip.
#' @param matchvar_fuzzy Character vector. Caliper matching variables
#'   (numeric only). Use \code{NULL} to skip.
#' @param caliper_fuzzy Numeric vector. Caliper bandwidths for
#'   \code{matchvar_fuzzy}. Must have the same length.
#' @param covariates_ps Character vector. Covariates for the PS model.
#'   Use \code{NULL} to skip.
#' @param caliper_ps_sd Numeric scalar. PS caliper in SD units.
#'   Default: 0.2.
#' @param ratio Integer. Controls per exposed case. Default: 1.
#' @param window_time Numeric. Minimum interval (days) between A and B.
#'   Default: 0.
#' @param adjusted_D Character vector. Disease names as candidate
#'   confounders (C is automatically excluded).
#' @param ncore Integer. Number of CPU cores for parallel execution.
#'   Default: 1.
#' @param cv_for_lambda Logical. Compute CV-based lambda variants.
#'   Default: \code{FALSE}.
#' @param seed Integer. Base random seed (\code{seed + i} per pair).
#'   Default: 1.
#' @param parallel_backend Character. Parallel backend:
#'   \code{"auto"} (default), \code{"fork"}, or \code{"socket"}.
#' @param stop_on_error Logical. If \code{FALSE} (default), return
#'   \code{NA} rows for failing pairs.
#' @param ... Additional arguments (currently unused).
#' @return A data frame with the same structure as \code{getAB_batch()}.
getAB_Cfalse_batch <- function(data,id,D_all,As,Bs,C,start_date,matchvar_exact = NULL,
                               matchvar_fuzzy = NULL,caliper_fuzzy = NULL,covariates_ps = NULL,caliper_ps_sd = 0.2,
                               ratio = 1,window_time = 0,adjusted_D,ncore = 1,
                               cv_for_lambda = F,seed = 1,
                               parallel_backend = c("auto", "fork", "socket"),
                               stop_on_error = FALSE,...){
  # 参数标准化
  matchvar_exact  <- if (is.null(matchvar_exact))  character(0) else matchvar_exact
  matchvar_fuzzy  <- if (is.null(matchvar_fuzzy))  character(0) else matchvar_fuzzy
  caliper_fuzzy   <- if (is.null(caliper_fuzzy))   numeric(0)   else caliper_fuzzy
  covariates_ps   <- if (is.null(covariates_ps))   character(0) else covariates_ps
  if(length(matchvar_fuzzy) != length(caliper_fuzzy)) stop("matchvar_fuzzy and caliper_fuzzy have different length")

  # 解析疾病对
  pairs <- resolve_disease_pairs(data, D_all, As, Bs, adjusted_D, C = C)
  As <- pairs$As; Bs <- pairs$Bs; adjusted_D <- pairs$adjusted_D
  n_pair <- pairs$n_pair
  A_dates <- pairs$A_dates; B_dates <- pairs$B_dates
  A_events <- pairs$A_events; B_events <- pairs$B_events
  first_Ds <- pairs$first_Ds; second_Ds <- pairs$second_Ds

  # 选择并行后端
  parallel_backend <- match.arg(parallel_backend)
  if (parallel_backend == "auto") {
    parallel_backend <- if (.Platform$OS.type == "windows") "socket" else "fork"
  }

  if (parallel_backend == "fork") {
    # Linux/macOS: fork 模式，共享内存，无需序列化 data
    ncore <- min(ncore, parallel::detectCores() - 1L, n_pair)
    results <- pbapply::pblapply(seq_len(n_pair), function(i) {
      tryCatch({
        fit <- getAB_Cfalse(data, id, As[i], Bs[i],
          C, adjusted_D,
          start_date, matchvar_exact, matchvar_fuzzy, caliper_fuzzy,
          covariates_ps, caliper_ps_sd, ratio, window_time,
          output_match = FALSE, trace = FALSE, cv_for_lambda,
          seed = seed + i)
        
        if (i %% max(1, floor(n_pair / 500)) == 0 || i == n_pair)
          cat(sprintf("[%s] %d/%d (%.0f%%) %s -> %s\n", 
                      Sys.time(), i, n_pair, round(i/n_pair*100), first_Ds[i], second_Ds[i]),
              file = "logs/progress.txt", append = TRUE)
        
        tibble(first_D = first_Ds[i], second_D = second_Ds[i],
               n_exposure = fit$n_exposure, n_unmatched = fit$n_unmatched,
               n_control = fit$n_control, n_event = fit$n_event) %>%
          cbind(fit$res_coef) %>%
          { cbind(coef_type = rownames(.), ., error = NA_character_) }
      }, error = function(e) {
        if (stop_on_error) stop(e)
        C_event_col <- paste0(C, "_event")
        C_date_col  <- paste0(C, "_date")
        dA <- tryCatch(as.Date(data[[A_dates[i]]]), error = function(.) NULL)
        dB <- tryCatch(as.Date(data[[B_dates[i]]]), error = function(.) NULL)
        dC <- tryCatch(as.Date(data[[C_date_col]]),  error = function(.) NULL)
        valid <- if (!is.null(dA) && !is.null(dB) && !is.null(dC)) {
          as.integer(difftime(dC, dA, units = "days")) > 0 &
            as.integer(difftime(dB, dA, units = "days")) > window_time
        } else rep(FALSE, nrow(data))
        valid[is.na(valid)] <- FALSE
        tibble(
          coef_type = NA_character_, first_D = first_Ds[i],
          second_D = second_Ds[i],
          n_exposure = sum(data[[A_events[i]]][valid], na.rm = TRUE),
          n_unmatched = NA_integer_, n_control = NA_integer_,
          n_event = sum(data[[B_events[i]]][valid], na.rm = TRUE),
          coef = NA_real_, `exp(coef)` = NA_real_, `se(coef)` = NA_real_,
          z = NA_real_, `Pr(>|z|)` = NA_real_, error = conditionMessage(e))
      })
    }, cl = ncore)
    results <- rbindlist(results[!vapply(results, is.null, logical(1L))])

  } else {
    # Windows: socket 模式，函数和数据通过词法作用域传递
    ncore <- min(ncore, parallel::detectCores() - 1L, n_pair)
    cl <- snow::makeCluster(ncore, type = "SOCK")
    doSNOW::registerDoSNOW(cl)

    snow::clusterExport(cl, c("getAB_Cfalse", "run_eds_sampling",
                              "calculate_propensity_score",
                              "build_matchvar_categories",
                              "resolve_disease_pairs",
                              "cut_cat", "extractAICBIC.glmnet",
                              "select_lasso_parameters", "run_eds_sampling_rcpp",
                              "extract_coxph_coef",
                              "SCRIPT_DIR"), envir = environment())
    snow::clusterEvalQ(cl, {
      library(data.table)
      library(survival)
      library(glmnet)
      library(dplyr)
      library(tidyr)
      library(Rcpp)
      suppressMessages(source(file.path(SCRIPT_DIR, "utils.R")))
    })

    pb <- txtProgressBar(max = n_pair, style = 3)
    progress <- function(nn) setTxtProgressBar(pb, nn)
    opts <- list(progress = progress)

    results <- foreach(i = seq_len(n_pair),
                       .combine = function(x, y) rbindlist(list(x, y)),
                       .options.snow = opts) %dopar% {
      tryCatch({
        fit <- getAB_Cfalse(data, id, As[i], Bs[i],
          C, adjusted_D,
          start_date, matchvar_exact, matchvar_fuzzy, caliper_fuzzy,
          covariates_ps, caliper_ps_sd, ratio, window_time,
          output_match = FALSE, trace = FALSE, cv_for_lambda,
          seed = seed + i)
        tibble(first_D = first_Ds[i], second_D = second_Ds[i],
               n_exposure = fit$n_exposure, n_unmatched = fit$n_unmatched,
               n_control = fit$n_control, n_event = fit$n_event) %>%
          cbind(fit$res_coef) %>%
          { cbind(coef_type = rownames(.), ., error = NA_character_) }
      }, error = function(e) {
        if (stop_on_error) stop(e)
        C_event_col <- paste0(C, "_event")
        C_date_col  <- paste0(C, "_date")
        dA <- tryCatch(as.Date(data[[A_dates[i]]]), error = function(.) NULL)
        dB <- tryCatch(as.Date(data[[B_dates[i]]]), error = function(.) NULL)
        dC <- tryCatch(as.Date(data[[C_date_col]]),  error = function(.) NULL)
        valid <- if (!is.null(dA) && !is.null(dB) && !is.null(dC)) {
          as.integer(difftime(dC, dA, units = "days")) > 0 &
            as.integer(difftime(dB, dA, units = "days")) > window_time
        } else rep(FALSE, nrow(data))
        valid[is.na(valid)] <- FALSE
        tibble(
          coef_type = NA_character_, first_D = first_Ds[i],
          second_D = second_Ds[i],
          n_exposure = sum(data[[A_events[i]]][valid], na.rm = TRUE),
          n_unmatched = NA_integer_, n_control = NA_integer_,
          n_event = sum(data[[B_events[i]]][valid], na.rm = TRUE),
          coef = NA_real_, `exp(coef)` = NA_real_, `se(coef)` = NA_real_,
          z = NA_real_, `Pr(>|z|)` = NA_real_, error = conditionMessage(e))
      })
    }
    close(pb)
    snow::stopCluster(cl)
  }

  rownames(results) <- NULL
  as.data.frame(results)
}



#' Estimate the A x C interaction effect on B
#'
#' Uses exposure density sampling (EDS) to estimate the interaction effect
#' of disease A and disease C on the risk of disease B. In the LASSO-Cox
#' model, \code{exposure}, \code{C}, and \code{exposure:C} are all
#' unpenalized. Uses standard EDS (no C-conditioning on risk sets).
#'
#' @param data A data frame containing patient-level data.
#' @param id Character. The column name of the unique patient identifier.
#' @param A Character. The name of the first disease (exposure).
#' @param B Character. The name of the second disease (outcome).
#' @param C Character. The name of the interacting disease.
#' @param start_date Character. Column name of the study start date (Date).
#' @param matchvar_exact Character vector. Exact matching variables.
#'   Use \code{NULL} to skip.
#' @param matchvar_fuzzy Character vector. Caliper matching variables
#'   (numeric only). Use \code{NULL} to skip.
#' @param caliper_fuzzy Numeric vector. Caliper bandwidths for
#'   \code{matchvar_fuzzy}. Must have the same length.
#' @param covariates_ps Character vector. Covariates for the PS model.
#'   Use \code{NULL} to skip.
#' @param caliper_ps_sd Numeric scalar. PS caliper in SD units.
#'   Default: 0.2.
#' @param ratio Integer. Controls per exposed case. Default: 1.
#' @param window_time Numeric. Minimum interval (days) between A and B.
#'   Default: 0.
#' @param adjusted_D Character vector. Disease names as candidate
#'   confounders (A, C, and A:C are unpenalized in LASSO).
#' @param output_match Logical. If \code{TRUE}, include matched dataset
#'   in output. Default: \code{FALSE}.
#' @param trace Logical. Print EDS progress. Default: \code{FALSE}.
#' @param cv_for_lambda Logical. Compute CV-based lambda variants.
#'   Default: \code{FALSE}.
#' @param seed Integer. Random seed. Default: 1.
#' @return A list with \code{res_coef} containing 3 rows per model variant
#'   (exposure, C, exposure:C) instead of 1.
getAB_C_interact <- function(data,id,A,B,C,adjusted_D,start_date,matchvar_exact = NULL,
                         matchvar_fuzzy = NULL,caliper_fuzzy = NULL,covariates_ps = NULL,caliper_ps_sd = 0.2,
                         ratio = 1,window_time = 0,output_match = FALSE,trace = FALSE, cv_for_lambda = F,seed = 1){
  # 参数标准化：NULL → 零长度向量
  matchvar_exact  <- if (is.null(matchvar_exact))  character(0) else matchvar_exact
  matchvar_fuzzy  <- if (is.null(matchvar_fuzzy))  character(0) else matchvar_fuzzy
  caliper_fuzzy   <- if (is.null(caliper_fuzzy))   numeric(0)   else caliper_fuzzy
  covariates_ps   <- if (is.null(covariates_ps))   character(0) else covariates_ps
  data <- as.data.frame(data)
  if(any(sapply(c(id,A,B,C,start_date,ratio,caliper_ps_sd),length) != 1))
    stop("The length of id,A,B,C,start_date,ratio,caliper_ps_sd must be 1")
  if(length(matchvar_fuzzy) != length(caliper_fuzzy)) stop("matchvar_fuzzy and caliper_fuzzy have different length")
  
  A_date <- paste0(A,"_date")
  B_date <- paste0(B,"_date")
  C_date <- paste0(C,"_date")
  A_event <- paste0(A,"_event")
  B_event <- paste0(B,"_event")
  C_event <- paste0(C,"_event")
  adjusted_D <- setdiff(adjusted_D,c(A,B,C))
  adjusted_D_date <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_date") else character(0)
  adjusted_D_event <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_event") else character(0)
  
  if (!start_date %in% names(data)) stop("start_date column not found in data")
  if (class(data[[start_date]]) != "Date") stop("start_date must be Date type")
  if (any(data[[start_date]] > data[[A_date]] | data[[start_date]] > data[[B_date]], na.rm = TRUE))
    stop("start_date must be <= A_date and B_date for all rows")
  
  keepvars <- c(id,A_date,A_event,B_date,B_event,C_date,C_event,start_date,matchvar_exact,
                matchvar_fuzzy,covariates_ps,adjusted_D_date,adjusted_D_event)
  if(!all(keepvars %in% names(data))) 
    stop("There were variables not being included in the data")
  if(any(sapply(data[,c(A_date,B_date,C_date,start_date,adjusted_D_date)], class) != "Date")) stop("There are date variables not being Date type")
  if(class(data[,A_event]) == "logical") data[,A_event] <- as.integer(data[,A_event])
  if(class(data[,B_event]) == "logical") data[,B_event] <- as.integer(data[,B_event])
  if(class(data[,C_event]) == "logical") data[,C_event] <- as.integer(data[,C_event])
  if(any(duplicated(c(A,B,C)))) stop("There is duplicated variable for A, B and C")
  if(!all(sapply(data[,matchvar_fuzzy,drop = FALSE], class) %in% c("integer","numeric"))) stop("matchvar_fuzzy must be a numeric or integer variable")
  
  origin_data <- na.omit(data[,c(id,A_date,A_event,B_date,B_event,C_date, C_event,start_date,matchvar_exact,
                                 matchvar_fuzzy,covariates_ps,adjusted_D_date,adjusted_D_event)])  
  for (i in adjusted_D) {
    origin_data[,paste0(i,"_time")] <- difftime(origin_data[,paste0(i,"_date")],origin_data[,start_date],units = "days") %>% as.integer
  }
  
  data <- origin_data[,c(id,A_date,A_event,B_date,B_event,C_date,C_event,start_date,matchvar_exact,matchvar_fuzzy,covariates_ps)]  
  data[,"A_time"] <- difftime(data[,A_date],data[,start_date],units = "days") %>% as.integer
  data[,"B_time"] <- difftime(data[,B_date],data[,start_date],units = "days") %>% as.integer
  data[,"C_time"] <- difftime(data[,C_date],data[,start_date],units = "days") %>% as.integer
  data[,"A_event"] <- data[,A_event]
  data[,"B_event"] <- data[,B_event]
  data[,"C_event"] <- data[,C_event]
  

  # 根据matchvar_fuzzy划分类别并计算倾向性评分（支持空参数）
  if (length(matchvar_fuzzy) > 0) {
    data <- build_matchvar_categories(data, matchvar_fuzzy, caliper_fuzzy)
  }
  # 占位符：确保 matchvar_names 永不为空
  if (length(matchvar_exact) == 0) {
    data[["_placeholder_exact"]] <- 1L
  }
  if (length(matchvar_fuzzy) == 0) {
    data[["_placeholder_fuzzy_cat"]] <- 1L
  }
  matchvar_fuzzy_cat <- if (length(matchvar_fuzzy) > 0) paste0(matchvar_fuzzy, "_cat") else character(0)
  if (length(covariates_ps) > 0) {
    data <- calculate_propensity_score(data, covariates_ps, matchvar_exact,
                                        matchvar_fuzzy_cat, caliper_ps_sd,
                                        ps_time_var = "A_time",
                                        ps_filter_expr = "A_time > 0 & B_time > 0 & B_time >= A_time")
  } else {
    data[["ps"]] <- 1
    data[["ps_cat"]] <- 1L
  }

  
  
  ######## EDS抽样
  matchvar_names <- c(
    if (length(matchvar_exact) > 0) matchvar_exact else "_placeholder_exact",
    if (length(matchvar_fuzzy) > 0) paste0(matchvar_fuzzy, "_cat") else "_placeholder_fuzzy_cat",
    "ps_cat"
  )
  eds_result <- run_eds_sampling_rcpp(
    data = data, matchvar_names = matchvar_names,
    ratio = ratio, window_time = window_time, seed = seed,
    exposure_filter_expr = "A_time + window_time < B_time & A_event == 1",
    control_filter_expr = "A_time > start_time_i + window_time & B_time > start_time_i + window_time",
    id = id, include_C_time = FALSE, trace = trace
  )
  data_matched_list <- eds_result$data_matched_list
  exposure_dt <- eds_result$exposure_dt
  
  # 剔除没有匹配到对照的list
  NA_index <- sapply(data_matched_list, function(x) anyNA(x[[id]])) %>% which
  if(length(NA_index) >= 1) data_matched_list <- data_matched_list[-NA_index]
  
  
  adjusted_D_time <- if (length(adjusted_D) > 0) paste0(adjusted_D, "_time") else character(0)
  data_matched <- rbindlist(data_matched_list) %>%
    left_join(data[,c(id,matchvar_names,matchvar_fuzzy,"ps","A_time","B_time","A_event","B_event","C_time","C_event")],by = id) %>%
    left_join(origin_data[,c(id,adjusted_D_time,adjusted_D_event),drop = FALSE],by = id) %>%
    group_by(!!!syms(matchvar_names)) %>%
    mutate(strata = cur_group_id()) %>%
    ungroup()
  
  # 处理对照组中B发生在A之前的：只保留区间1，因为区间2已经在exposure组了 
  data_matched <- data_matched %>% mutate(
    end_time = fifelse(
      exposure == 0 & A_time < B_time, A_time,
      B_time),
    status = fifelse(
      exposure == 0 & A_time < B_time, 0L,
      B_event),
  )
  
  # 去掉strata内部的重复个体，保留start_time最早的
  data_matched <- data_matched %>% 
    group_by(strata,!!sym(id),exposure) %>%
    slice_min(start_time, n = 1, with_ties = FALSE) %>%  # 保留 start_time 最小的行
    ungroup()
  
  # 去掉没有事件发生的strata
  data_matched_model <- data_matched %>%
    group_by(strata) %>%
    filter(any(status == 1)) %>%  # 保留至少有一个事件的 strata
    ungroup()
  
  # 添加其它condition的时依协变量
  data_matched_model$unique_row <- 1:nrow(data_matched_model)
  tdc_expr <- paste0(
    "tmerge(data_matched_model, data_matched_model, id=unique_row,tstart = start_time,tstop = end_time,status = event(end_time,status),",
    paste0(paste0("`", c("C",adjusted_D), "` = tdc(`", c("C",adjusted_D), "_time`)"),collapse = ","),")"
  )
  td_matched_data <- eval(parse(text = tdc_expr)) 
  
  n_unmatched <- nrow(exposure_dt) - sum(data_matched$exposure == 1)
  if(n_unmatched > 0) warning(sprintf("There are %.0f participants having not been matched for %.0f participants",n_unmatched,nrow(exposure_dt)))
  n_exposure <- nrow(exposure_dt)
  n_control <- sum(data_matched$exposure==0)
  n_event <- sum(data_matched$status==1)
  #### fit model: cox + lasso

  # not adjusting for other conditions
  fit_unadjusted_D <- coxph(Surv(tstart,tstop,status) ~ exposure*C + strata(strata),data = td_matched_data)
  coef_unadjusted_D <- extract_coxph_coef(fit_unadjusted_D, c("exposure","C","exposure:C"))
  
  # adjusting for other conditions
  n_total <- nrow(td_matched_data)
  params <- select_lasso_parameters(n_total)
  nlambda_val <- params$nlambda_val
  lambda_min_ratio_val <- params$lambda_min_ratio_val
  nfolds_val <- params$nfolds_val
  control_list <- params$control_list
  
  
  td_matched_data$exposure_C <- td_matched_data$exposure * td_matched_data$C
  x0 <- td_matched_data[,c("exposure","C","exposure_C",adjusted_D)] %>% as.matrix()
  y0 <- stratifySurv(Surv(td_matched_data$tstart,td_matched_data$tstop,td_matched_data$status),strata = td_matched_data$strata)
  # fit_lasso <- glmnet(x = x0,y = y0,family = "cox",penalty.factor = c(0,0,0,rep(1, length(adjusted_D))),relax = T,nlambda = 50)
  args_list <- list(
    x = x0,y = y0,family = "cox",penalty.factor = c(0,0,0,rep(1, length(adjusted_D))),
    relax = T,nlambda = nlambda_val,cox.ties = "breslow",
    lambda.min.ratio = lambda_min_ratio_val,control = control_list,
    standardize = FALSE
  )
  fit_lasso <- do.call(glmnet,args_list)
  
  
  lambda_aic_bic <- extractAICBIC.glmnet(fit_lasso,n_event)
  
  coef_lasso <- coef(fit_lasso,s = lambda_aic_bic$lambda.AIC)
  vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
  fit_cox <- coxph(cox_formula,data = td_matched_data)
  coef_AIC_lambda <- extract_coxph_coef(fit_cox, c("exposure","C","exposure_C")) %>% {rownames(.)[3L] <- "exposure:C"; .}

  
  coef_lasso <- coef(fit_lasso,s = lambda_aic_bic$lambda.BIC)
  vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
  fit_cox <- coxph(cox_formula,data = td_matched_data)
  coef_BIC_lambda <- extract_coxph_coef(fit_cox, c("exposure","C","exposure_C")) %>% {rownames(.)[3L] <- "exposure:C"; .}
  
  if(cv_for_lambda){
    # fit_lasso_cv <- cv.glmnet(x = x0,y = y0,family = "cox",penalty.factor = c(0,0,0,rep(1, length(adjusted_D))),nfolds = 5,relax = T,gamma = 0,nlambda = 30)
    fit_lasso_cv <- do.call(cv.glmnet,modifyList(args_list,val = list(gamma = 0,nfolds = nfolds_val)))
    coef_lasso <- coef(fit_lasso_cv,s = fit_lasso_cv$lambda.min)
    vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
    fit_cox <- coxph(cox_formula,data = td_matched_data)
    coef_cv_lambda.min <- extract_coxph_coef(fit_cox, c("exposure","C","exposure_C")) %>% {rownames(.)[3L] <- "exposure:C"; .}
    
    coef_lasso <- coef(fit_lasso_cv,s = fit_lasso_cv$lambda.1se)
    vars_remain <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    cox_formula <- paste0("Surv(tstart,tstop,status) ~",paste0("`", vars_remain, "`", collapse = "+"),"+strata(strata)") %>% as.formula()
    fit_cox <- coxph(cox_formula,data = td_matched_data)
    coef_cv_lambda.1se <- extract_coxph_coef(fit_cox, c("exposure","C","exposure_C")) %>% {rownames(.)[3L] <- "exposure:C"; .}
  } else{
    coef_cv_lambda.1se <- coef_cv_lambda.min <- matrix(NA,ncol = 5,nrow = 3) %>% {rownames(.) <- c("exposure","C","exposure:C");.}
  }
  
  res_coef <- rbind(coef_unadjusted_D,coef_AIC_lambda,coef_BIC_lambda,coef_cv_lambda.min,coef_cv_lambda.1se) 
  res_coef <- cbind(coef_type = rep(c("coef_unadjusted_D","coef_AIC_lambda","coef_BIC_lambda","coef_cv_lambda.min","coef_cv_lambda.1se"),each = 3),
                    variable = rownames(res_coef),as.data.frame(res_coef)) %>% {rownames(.) <- NULL;.}
  
  if(output_match) {
    return(list(res_coef = res_coef,n_exposure = n_exposure,n_control = n_control,
                 n_unmatched = n_unmatched,n_event = n_event,data_matched = data_matched))
  } else {
    return(list(res_coef = res_coef,n_exposure = n_exposure,n_control = n_control,
                 n_unmatched = n_unmatched,n_event = n_event))
  }
}



#' Batch estimation for A x C interaction effects
#'
#' Runs \code{getAB_C_interact()} for multiple disease pairs in parallel,
#' estimating the interaction effect of A and C on the risk of B.
#'
#' @param data A data frame containing patient-level data.
#' @param id Character. The column name of the unique patient identifier.
#' @param D_all Character vector. All disease names for permutation.
#' @param As Character vector. First disease of each pair.
#' @param Bs Character vector. Second disease of each pair.
#' @param C Character. The interacting disease name.
#' @param start_date Character. Column name of the study start date (Date).
#' @param matchvar_exact Character vector. Exact matching variables.
#'   Use \code{NULL} to skip.
#' @param matchvar_fuzzy Character vector. Caliper matching variables
#'   (numeric only). Use \code{NULL} to skip.
#' @param caliper_fuzzy Numeric vector. Caliper bandwidths for
#'   \code{matchvar_fuzzy}. Must have the same length.
#' @param covariates_ps Character vector. Covariates for the PS model.
#'   Use \code{NULL} to skip.
#' @param caliper_ps_sd Numeric scalar. PS caliper in SD units.
#'   Default: 0.2.
#' @param ratio Integer. Controls per exposed case. Default: 1.
#' @param window_time Numeric. Minimum interval (days) between A and B.
#'   Default: 0.
#' @param adjusted_D Character vector. Disease names as candidate
#'   confounders (A, C, and A:C are unpenalized in LASSO).
#' @param ncore Integer. Number of CPU cores for parallel execution.
#'   Default: 1.
#' @param cv_for_lambda Logical. Compute CV-based lambda variants.
#'   Default: \code{FALSE}.
#' @param seed Integer. Base random seed (\code{seed + i} per pair).
#'   Default: 1.
#' @param parallel_backend Character. Parallel backend:
#'   \code{"auto"} (default), \code{"fork"}, or \code{"socket"}.
#' @param stop_on_error Logical. If \code{FALSE} (default), return
#'   \code{NA} rows for failing pairs.
#' @param ... Additional arguments (currently unused).
#' @return A data frame with 3 rows per disease pair (exposure, C,
#'   exposure:C) and the same columns as \code{getAB_batch()}.
getAB_C_interact_batch <- function(data,id,D_all,As,Bs,C,start_date,matchvar_exact = NULL,
                               matchvar_fuzzy = NULL,caliper_fuzzy = NULL,covariates_ps = NULL,caliper_ps_sd = 0.2,
                               ratio = 1,window_time = 0,adjusted_D,ncore = 1,
                               cv_for_lambda = F,seed = 1,
                               parallel_backend = c("auto", "fork", "socket"),
                               stop_on_error = FALSE,...){
  # 参数标准化
  matchvar_exact  <- if (is.null(matchvar_exact))  character(0) else matchvar_exact
  matchvar_fuzzy  <- if (is.null(matchvar_fuzzy))  character(0) else matchvar_fuzzy
  caliper_fuzzy   <- if (is.null(caliper_fuzzy))   numeric(0)   else caliper_fuzzy
  covariates_ps   <- if (is.null(covariates_ps))   character(0) else covariates_ps
  if(length(matchvar_fuzzy) != length(caliper_fuzzy)) stop("matchvar_fuzzy and caliper_fuzzy have different length")

  # 解析疾病对
  pairs <- resolve_disease_pairs(data, D_all, As, Bs, adjusted_D, C = C)
  As <- pairs$As; Bs <- pairs$Bs; adjusted_D <- pairs$adjusted_D
  n_pair <- pairs$n_pair
  A_dates <- pairs$A_dates; B_dates <- pairs$B_dates
  A_events <- pairs$A_events; B_events <- pairs$B_events
  first_Ds <- pairs$first_Ds; second_Ds <- pairs$second_Ds

  # 选择并行后端
  parallel_backend <- match.arg(parallel_backend)
  if (parallel_backend == "auto") {
    parallel_backend <- if (.Platform$OS.type == "windows") "socket" else "fork"
  }

  if (parallel_backend == "fork") {
    # Linux/macOS: fork 模式，共享内存，无需序列化 data
    ncore <- min(ncore, parallel::detectCores() - 1L, n_pair)
    results <- pbapply::pblapply(seq_len(n_pair), function(i) {
      tryCatch({
        fit <- getAB_C_interact(data, id, As[i], Bs[i],
          C, adjusted_D,
          start_date, matchvar_exact, matchvar_fuzzy, caliper_fuzzy,
          covariates_ps, caliper_ps_sd, ratio, window_time,
          output_match = FALSE, trace = FALSE, cv_for_lambda,
          seed = seed + i)
        
        if (i %% max(1, floor(n_pair / 500)) == 0 || i == n_pair)
          cat(sprintf("[%s] %d/%d (%.0f%%) %s -> %s\n", 
                      Sys.time(), i, n_pair, round(i/n_pair*100), first_Ds[i], second_Ds[i]),
              file = "logs/progress.txt", append = TRUE)
        
        tibble(first_D = first_Ds[i], second_D = second_Ds[i],
               n_exposure = fit$n_exposure, n_unmatched = fit$n_unmatched,
               n_control = fit$n_control, n_event = fit$n_event) %>%
          cbind(fit$res_coef) %>%
          { cbind(coef_type = rownames(.), ., error = NA_character_) }
      }, error = function(e) {
        if (stop_on_error) stop(e)
        dA <- tryCatch(as.Date(data[[A_dates[i]]]), error = function(.) NULL)
        dB <- tryCatch(as.Date(data[[B_dates[i]]]), error = function(.) NULL)
        valid <- if (!is.null(dA) && !is.null(dB)) {
          idx <- as.integer(difftime(dB, dA, units = "days")) > window_time
          idx[is.na(idx)] <- FALSE; idx
        } else rep(FALSE, nrow(data))
        tibble(
          coef_type = NA_character_, first_D = first_Ds[i],
          second_D = second_Ds[i],
          n_exposure = sum(data[[A_events[i]]][valid], na.rm = TRUE),
          n_unmatched = NA_integer_, n_control = NA_integer_,
          n_event = sum(data[[B_events[i]]][valid], na.rm = TRUE) * (ratio + 1),
          coef = NA_real_, `exp(coef)` = NA_real_, `se(coef)` = NA_real_,
          z = NA_real_, `Pr(>|z|)` = NA_real_, error = conditionMessage(e))
      })
    }, cl = ncore)
    results <- rbindlist(results[!vapply(results, is.null, logical(1L))])

  } else {
    # Windows: socket 模式，clusterExport 导出 data（仅一次）
    ncore <- min(ncore, parallel::detectCores() - 1L, n_pair)
    cl <- snow::makeCluster(ncore, type = "SOCK")
    doSNOW::registerDoSNOW(cl)

    # 导出函数和数据（data 通过 clusterExport 而非词法作用域，避免重复序列化）
    snow::clusterExport(cl, c("getAB_C_interact", "run_eds_sampling",
                              "calculate_propensity_score",
                              "build_matchvar_categories",
                              "resolve_disease_pairs",
                              "cut_cat", "extractAICBIC.glmnet",
                              "select_lasso_parameters", "run_eds_sampling_rcpp",
                              "extract_coxph_coef",
                              "SCRIPT_DIR"), envir = environment())
    snow::clusterEvalQ(cl, {
      library(tidyverse)
      library(data.table)
      library(survival)
      library(glmnet)
      library(Rcpp)
      Rcpp::sourceCpp(file.path(SCRIPT_DIR, "eds_rcpp.cpp"))
    })

    pb <- txtProgressBar(max = n_pair, style = 3)
    progress <- function(nn) setTxtProgressBar(pb, nn)
    opts <- list(progress = progress)

    results <- foreach(i = seq_len(n_pair),
                       .combine = function(x, y) rbindlist(list(x, y)),
                       .options.snow = opts) %dopar% {
      tryCatch({
        fit <- getAB_C_interact(data, id, As[i], Bs[i],
          C, adjusted_D,
          start_date, matchvar_exact, matchvar_fuzzy, caliper_fuzzy,
          covariates_ps, caliper_ps_sd, ratio, window_time,
          output_match = FALSE, trace = FALSE, cv_for_lambda,
          seed = seed + i)
        tibble(first_D = first_Ds[i], second_D = second_Ds[i],
               n_exposure = fit$n_exposure, n_unmatched = fit$n_unmatched,
               n_control = fit$n_control, n_event = fit$n_event) %>%
          cbind(fit$res_coef) %>%
          { cbind(coef_type = rownames(.), ., error = NA_character_) }
      }, error = function(e) {
        if (stop_on_error) stop(e)
        dA <- tryCatch(as.Date(data[[A_dates[i]]]), error = function(.) NULL)
        dB <- tryCatch(as.Date(data[[B_dates[i]]]), error = function(.) NULL)
        valid <- if (!is.null(dA) && !is.null(dB)) {
          idx <- as.integer(difftime(dB, dA, units = "days")) > window_time
          idx[is.na(idx)] <- FALSE; idx
        } else rep(FALSE, nrow(data))
        tibble(
          coef_type = NA_character_, first_D = first_Ds[i],
          second_D = second_Ds[i],
          n_exposure = sum(data[[A_events[i]]][valid], na.rm = TRUE),
          n_unmatched = NA_integer_, n_control = NA_integer_,
          n_event = sum(data[[B_events[i]]][valid], na.rm = TRUE) * (ratio + 1),
          coef = NA_real_, `exp(coef)` = NA_real_, `se(coef)` = NA_real_,
          z = NA_real_, `Pr(>|z|)` = NA_real_, error = conditionMessage(e))
      })
    }
    close(pb)
    snow::stopCluster(cl)
  }

  rownames(results) <- NULL
  as.data.frame(results)
}












