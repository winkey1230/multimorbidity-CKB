###############################################################################-#
###### get pair effects for any two diseases
###############################################################################-#
# args <- c("male","4","True")
args <- commandArgs(trailingOnly = TRUE) 
if(length(args) != 3) stop("The length of args is not 3")
sex <- args[1]  # male or female
ncores <- as.integer(args[2]) # integer
exclude_diagnose_at_baseline <- as.logical(args[3]) # True or False
window_time_day <- 180


cat("sex:",sex,"\n")
cat("ncores:",ncores,"\n")
cat("exclude_diagnose_at_baseline:",exclude_diagnose_at_baseline,"\n")
cat("window_time_day",window_time_day,"\n")

setwd("/public/home/hanyuting/2025-08-07-MM-trajectory/3.WW0823")
source("functions/comorbidity-fun.R")
data0_all <- fread("../1.Data/traj_260823.csv")

# all disease names
Dall <- names(data0_all)[grepl("_date$", names(data0_all))] %>%
  setdiff(c("study_date","censoring_date","age_at_study_date")) %>%
  gsub(pattern = "_date$",replacement = "")

data_sex <- data0_all %>%
  mutate(across((ends_with("date") & !matches("age_at_study_date")) | matches("dob_anon"), 
                ~ as.Date(., origin = "1970-01-01"))) %>%
  mutate(across(c(highest_education, smoking_4groups, alcohol_5groups, diet_5groups), as.factor)) %>%
  filter(is_female == ifelse(sex == "female",1,0)) %>%
  mutate(across(all_of(paste0(Dall, "_date")), 
                ~ if_else(. <= study_date, study_date, .))) %>%
  rename_with(~ paste0(., "_event"), all_of(Dall))


ps_vars <- c("highest_education", "smoking_4groups", "alcohol_5groups", "diet_5groups", "bmi","waist", "met")
match_exact <- "region_code"


if(exclude_diagnose_at_baseline){
  data_sex$start_date <- data_sex$study_date
} else {
  data_sex$start_date <- data_sex$study_date - 1
}


fit_batch_sex <- getAB_batch(
  data           = data_sex,
  id             = "csid",
  D_all          = Dall,
  start_date     = "start_date",
  adjusted_D     = setdiff(Dall,"death"),
  matchvar_exact = match_exact,
  matchvar_fuzzy = "age_at_study_date",
  caliper_fuzzy  = 4,
  covariates_ps  = ps_vars,
  caliper_ps_sd  = 0.2,
  ratio          = 4,
  window_time    = window_time_day,
  ncore          = ncores,
  cv_for_lambda  = F,
  seed           = 1
)
save(fit_batch_sex,file = paste0("results/pairs-effect-",sex,"-exclude",ifelse(exclude_diagnose_at_baseline,1,0),".Rdata"))
write.csv(fit_batch_sex,file = paste0("results/pairs-effect-",sex,"-exclude",ifelse(exclude_diagnose_at_baseline,1,0),".csv"))








