###############################################################################
############### extract trajectories
###############################################################################
library(tidyverse)
library(openxlsx)
library(data.table)
setwd("/public/home/hanyuting/2025-08-07-MM-trajectory/3.WW0823")
source("functions/extract-trajectories-fun.R")

# args <- c("male","True")
args <- commandArgs(trailingOnly = TRUE) 
if(length(args) != 2 & length(args) != 0) stop("The length of args is not 2")
sex <- ifelse(length(args) == 0,"male",args[1])    # male or female
exclude_diagnose_at_baseline <- ifelse(length(args) == 0,TRUE,as.logical(args[2])) # True or False


n_at_least <- 100
coef_type_target <- "coef_BIC_lambda"

data0_all <- fread("../1.Data/traj_260823.csv")
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
  rename_with(~ paste0(., "_event"), all_of(Dall)) %>% as.data.frame


data_target_pair <- read.csv(paste0("results/pairs-effect-",sex,"-exclude",ifelse(exclude_diagnose_at_baseline,1,0),".csv")) %>%
  filter(coef_type == coef_type_target | is.na(coef_type)) %>%
  filter(n_event >= n_at_least & !is.na(coef_type)) %>%
  rename(p_effect = Pr...z..,HR = exp.coef.) %>%
  filter(p_effect < 0.05/nrow(.)) %>%
  filter(HR > 1) %>%
  select(first_D, second_D)

# ----------------------------------------------------------------------------
# length = 2,3,4,5
# ----------------------------------------------------------------------------
start_disease <- NULL
n_minimum <- 10
trace_step <- 200
max_trajectories = Inf

traj_length <- 2
data_trajectory2 <- extract_trajectory(
  Ddata              = data_target_pair,
  idata              = data_sex,
  length             = traj_length,
  start_date         = "study_date",
  max_trajectories   = max_trajectories,
  n_minimum = n_minimum,
  trace_step = trace_step,
  start_disease = start_disease
)

traj_length <- 3
data_trajectory3 <- extract_trajectory(
  Ddata              = data_target_pair,
  idata              = data_sex,
  length             = traj_length,
  start_date         = "study_date",
  max_trajectories   = max_trajectories,
  trace_step = trace_step,
  n_minimum = n_minimum,
  data_traj_previous = data_trajectory2[,paste0("D",1:(traj_length-1))],
  start_disease = start_disease
)

traj_length <- 4
data_trajectory4 <- extract_trajectory(
  Ddata              = data_target_pair,
  idata              = data_sex,
  length             = traj_length,
  start_date         = "study_date",
  max_trajectories   = max_trajectories,
  trace_step = trace_step,
  n_minimum = n_minimum,
  data_traj_previous = data_trajectory3[,paste0("D",1:(traj_length-1))],
  start_disease = start_disease
)

traj_length <- 5
data_trajectory5 <- extract_trajectory(
  Ddata              = data_target_pair,
  idata              = data_sex,
  length             = traj_length,
  start_date         = "study_date",
  max_trajectories   = max_trajectories,
  trace_step = trace_step,
  n_minimum = n_minimum,
  data_traj_previous = data_trajectory4[,paste0("D",1:(traj_length-1))],
  start_disease = start_disease
)

data_traj <- data_trajectory4
data_interval <- data_traj[,grepl("median$",names(data_traj))] 
data_interval_p <- data_traj[,grepl("^p_",names(data_traj))] 
data_dif <- apply(data_interval, 1, diff) %>% t
mean(data_dif > 0 & data_interval_p < 0.05)
mean(data_dif < 0 & data_interval_p < 0.05)


write.xlsx(list(length2 = data_trajectory2,length3 = data_trajectory3,
                    length4 = data_trajectory4,length5 = data_trajectory5),
           paste0("results/data-trajectory-",sex,".xlsx"))

