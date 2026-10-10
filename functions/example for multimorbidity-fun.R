###############################################################################-#
# 暴露密度抽样（EDS）示例 — 共病因果效应估计                                -#
#                                                                              -#
# 演示内容：                                                                   -#
#   1. 利用UKB的示列数据（9 种疾病 + 协变量）                                  -#
#   2. 单对疾病分析：getAB(D50 → D64)                                          -#
#   3. getAB_batch：全部排列（D_all）                                         -#
#   4. getAB_batch：指定疾病对（As/Bs）                                       -#
#   5. 单对条件分析：Ctrue / Cfalse / Interact                               -#
#   6. 轨迹驱动批处理：全部排列（D_all）                                      -#
#   7. 轨迹驱动批处理：指定疾病对（As/Bs）                                    -#
###############################################################################-#

# ===========================================================================
# 0. 环境设置
# ===========================================================================
source("version5/comorbidity-fun.R")

load("data_sample.Rdata")
DT <- data_sample %>% rename(date_recruitment = date_recuitment)
ps_vars <- c("alcohol_status","is_urban","Townsend_deprivation","household_income",
             "education","ethnic","smoking_status","smoking_positive",
             "met_minutes","BMI","waist","diet_total_score")
D_all <- c("D50","D64","I10","I20","I21","I25","I26","I44","I48")
# ===========================================================================
# 1. 单对分析：D50 → D64
# ===========================================================================
cat("\n========== 1. 单对分析: D50 → D64 ==========\n")

system.time({
  fit_D50_D64 <- getAB(
    data           = DT,
    id             = "Participant_ID",
    A              = "D50",
    B              = "D64",
    start_date     = "date_recruitment",
    matchvar_exact = "region",
    matchvar_fuzzy = "age_at_recruitment",
    caliper_fuzzy  = 4,
    covariates_ps  = ps_vars,
    caliper_ps_sd  = 0.2,
    ratio          = 4,
    window_time    = 30,
    adjusted_D     = D_all,
    output_match   = FALSE,
    trace          = FALSE,
    cv_for_lambda  = TRUE,
    seed           = 1
  )
})

cat("\n--- 结果: D50 → D64 ---\n")
print(fit_D50_D64)


# ===========================================================================
# 2. getAB_batch：全部疾病排列（并行）
# ===========================================================================
cat("\n========== 2. getAB_batch: 全部疾病排列 ==========\n")

n_cores <- 5L  # 并行核心数
cat(sprintf("疾病: %s\n", paste(D_all, collapse = ", ")))
cat(sprintf("预计 %d 对, 并行核心: %d\n",
            length(D_all) * (length(D_all) - 1L), n_cores))

system.time({
  fit_batch_all <- getAB_batch(
    data           = DT,
    id             = "Participant_ID",
    D_all          = D_all,
    start_date     = "date_recruitment",
    matchvar_exact = "region",
    matchvar_fuzzy = "age_at_recruitment",
    caliper_fuzzy  = 4,
    covariates_ps  = ps_vars,
    caliper_ps_sd  = 0.2,
    ratio          = 4,
    window_time    = 30,
    ncore          = n_cores,
    cv_for_lambda  = F,
    seed           = 1
  )
})

cat(sprintf("\n--- 批量结果: %d 对疾病 ---\n", nrow(fit_batch_all)))
print(fit_batch_all, nrows = 30)


# ===========================================================================
# 3. getAB_batch：指定疾病对（并行）
# ===========================================================================
cat("\n========== 3. getAB_batch: 指定疾病对 ==========\n")

cat("疾病对: D50→D64, D50→I25, I25→I26\n")

system.time({
  fit_batch_spec <- getAB_batch(
    data           = DT,
    id             = "Participant_ID",
    As             = c("D50", "D50", "I25"),
    Bs             = c("D64", "I25", "I26"),
    start_date     = "date_recruitment",
    adjusted_D     = D_all,
    matchvar_exact = "region",
    matchvar_fuzzy = "age_at_recruitment",
    caliper_fuzzy  = 4,
    covariates_ps  = ps_vars,
    caliper_ps_sd  = 0.2,
    ratio          = 4,
    window_time    = 30,
    ncore          = n_cores,
    cv_for_lambda  = F,
    seed           = 1
  )
})

cat("\n--- 指定疾病对结果 ---\n")
print(fit_batch_spec)

cat("\n========== 完成 Part 1 ==========\n")

# ===========================================================================
# Part 2: 轨迹驱动条件分析 (trajectory-driven)
# ===========================================================================
cat("\n========== Part 2: 轨迹驱动条件分析 ==========\n")
# ---- 4. 单对条件分析: Ctrue / Cfalse / Interact ----
cat("\n--- 4. 单对条件分析 ---\n")

system.time({
  # Ctrue: 需C已发生后A对B的效应
  fit_CT <- getAB_Ctrue(
    data           = DT,
    id             = "Participant_ID",
    A              = "D50",
    B              = "D64",
    C              = "I10",
    start_date     = "date_recruitment",
    matchvar_exact = "region",
    matchvar_fuzzy = "age_at_recruitment",
    caliper_fuzzy  = 4,
    covariates_ps  = ps_vars,
    caliper_ps_sd  = 0.2,
    ratio          = 4,
    window_time    = 30,
    adjusted_D     = D_all,
    output_match   = FALSE,
    trace          = FALSE,
    cv_for_lambda  = F,
    seed           = 1
  )

  # Cfalse: 需C未发生时A对B的效应
  fit_CF <- getAB_Cfalse(
    data           = DT,
    id             = "Participant_ID",
    A              = "D50",
    B              = "D64",
    C              = "I10",
    start_date     = "date_recruitment",
    matchvar_exact = "region",
    matchvar_fuzzy = "age_at_recruitment",
    caliper_fuzzy  = 4,
    covariates_ps  = ps_vars,
    caliper_ps_sd  = 0.2,
    ratio          = 4,
    window_time    = 30,
    adjusted_D     = D_all,
    output_match   = FALSE,
    trace          = FALSE,
    cv_for_lambda  = F,
    seed           = 1
  )

  # Interact: A×C交互效应，D1⩒D2有因果结构，D3作为交互项
  fit_interact <- getAB_C_interact(
    data           = DT,
    id             = "Participant_ID",
    A              = "D50",
    B              = "D64",
    C              = "I10",
    start_date     = "date_recruitment",
    matchvar_exact = "region",
    matchvar_fuzzy = "age_at_recruitment",
    caliper_fuzzy  = 4,
    covariates_ps  = ps_vars,
    caliper_ps_sd  = 0.2,
    ratio          = 4,
    window_time    = 30,
    adjusted_D     = D_all,
    output_match   = FALSE,
    trace          = FALSE,
    cv_for_lambda  = F,
    seed           = 1
  )
})

cat("\n=== getAB_Ctrue (D50->D64 | C=I10 occurred) ===\n")
print(fit_CT)
cat("\n=== getAB_Cfalse (D50->D64 | C=I10 not occurred) ===\n")
print(fit_CF)
cat("\n=== getAB_C_interact (D50->D64 x C=I10) ===\n")
print(fit_interact)


# ---- 5. 轨迹驱动批处理：全部排列（D_all） ----
cat("\n--- 5. 轨迹驱动批处理: 全部排列 ---\n")

# Ctrue batch
fit_CT_batch <- getAB_Ctrue_batch(
  data           = DT,
  id             = "Participant_ID",
  D_all          = D_all,
  C = "I10",
  start_date     = "date_recruitment",
  matchvar_exact = "region",
  matchvar_fuzzy = "age_at_recruitment",
  caliper_fuzzy  = 4,
  covariates_ps  = ps_vars,
  caliper_ps_sd  = 0.2,
  ratio          = 4,
  window_time    = 30,
  ncore          = n_cores,
  cv_for_lambda  = F,
  seed           = 1
)
print(fit_CT_batch)

# Cfalse batch
fit_CF_batch <- getAB_Cfalse_batch(
  data           = DT,
  id             = "Participant_ID",
  D_all          = D_all,
  C = "I10",
  start_date     = "date_recruitment",
  matchvar_exact = "region",
  matchvar_fuzzy = "age_at_recruitment",
  caliper_fuzzy  = 4,
  covariates_ps  = ps_vars,
  caliper_ps_sd  = 0.2,
  ratio          = 4,
  window_time    = 30,
  ncore          = n_cores,
  cv_for_lambda  = F,
  seed           = 1
)

print(fit_CF_batch)

# Interact batch
fit_interact_batch <- getAB_C_interact_batch(
  data           = DT,
  id             = "Participant_ID",
  D_all          = D_all,
  C = "I10",
  start_date     = "date_recruitment",
  matchvar_exact = "region",
  matchvar_fuzzy = "age_at_recruitment",
  caliper_fuzzy  = 4,
  covariates_ps  = ps_vars,
  caliper_ps_sd  = 0.2,
  ratio          = 4,
  window_time    = 30,
  ncore          = n_cores,
  cv_for_lambda  = F,
  seed           = 1
)
print(fit_interact_batch)


# ---- 6. 轨迹驱动批处理：指定疾病对（As/Bs） ----
cat("\n--- 6. 轨迹驱动批处理: 指定疾病对 ---\n")

cat("疾病对: D50→D64, D50→I25, I25→I26, 条件疾病: I10\n")
n_cores3 <- max(1L, n_cores %/% 2L)

# Ctrue batch (指定)
fit_CT_batch_spec <- getAB_Ctrue_batch(
  data           = DT,
  id             = "Participant_ID",
  As             = c("D50", "D50", "I25"),
  Bs             = c("D64", "I25", "I26"),
  C = "I10",
  start_date     = "date_recruitment",
  adjusted_D     = D_all,
  matchvar_exact = "region",
  matchvar_fuzzy = "age_at_recruitment",
  caliper_fuzzy  = 4,
  covariates_ps  = ps_vars,
  caliper_ps_sd  = 0.2,
  ratio          = 4,
  window_time    = 30,
  ncore          = n_cores,
  cv_for_lambda  = F,
  seed           = 1
)
cat("\n=== getAB_Ctrue_batch (指定: D50→D64, D50→I25, I25→I26, C=I10) ===\n")
print(fit_CT_batch_spec)

# Cfalse batch (指定)
fit_CF_batch_spec <- getAB_Cfalse_batch(
  data           = DT,
  id             = "Participant_ID",
  As             = c("D50", "D50", "I25"),
  Bs             = c("D64", "I25", "I26"),
  C = "I10",
  start_date     = "date_recruitment",
  adjusted_D     = D_all,
  matchvar_exact = "region",
  matchvar_fuzzy = "age_at_recruitment",
  caliper_fuzzy  = 4,
  covariates_ps  = ps_vars,
  caliper_ps_sd  = 0.2,
  ratio          = 4,
  window_time    = 30,
  ncore          = n_cores,
  cv_for_lambda  = F,
  seed           = 1
)
cat("\n=== getAB_Cfalse_batch (指定: D50→D64, D50→I25, I25→I26, C=I10) ===\n")
print(fit_CF_batch_spec)

# Interact batch (指定)
fit_interact_batch_spec <- getAB_C_interact_batch(
  data           = DT,
  id             = "Participant_ID",
  As             = c("D50", "D50", "I25"),
  Bs             = c("D64", "I25", "I26"),
  C = "I10",
  start_date     = "date_recruitment",
  adjusted_D     = D_all,
  matchvar_exact = "region",
  matchvar_fuzzy = "age_at_recruitment",
  caliper_fuzzy  = 4,
  covariates_ps  = ps_vars,
  caliper_ps_sd  = 0.2,
  ratio          = 4,
  window_time    = 30,
  ncore          = n_cores,
  cv_for_lambda  = F,
  seed           = 1
)
print(fit_interact_batch_spec)


# ---- 7. 保存全部结果 ----
cat("\n========== 全部完成 ==========\n")
