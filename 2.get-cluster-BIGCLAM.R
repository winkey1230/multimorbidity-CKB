################################################################################
## BIGCLAM 
################################################################################
library(tidyverse)
library(readxl)
library(igraph)
library(openxlsx)
library(patchwork)
source("E:\\BaiduNetdiskWorkspace\\research\\multimorbidity\\bigclam-fun\\function\\bigclam_fun.R")
###### parameter setting-------
# args <- c("male","True")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 & length(args) != 0) stop("The length of args is not 2")
sex <- ifelse(length(args) == 0, "female", args[1])    # male or female
exclude_diagnose_at_baseline <- ifelse(length(args) == 0, TRUE, as.logical(args[2])) # True or False

n_at_least <- 100
coef_type_target <- "coef_BIC_lambda"
c_scan <- 2:8
scale_w <- 100
n_boot <- 50
n_starts <- 20
n_starts_boot <- 3
eps <- 1e-6
final_c <- 6
###### data ------
setwd("E:\\BaiduNetdiskWorkspace\\research\\multimorbidity\\CKB\\overall")
data_disease <- read_xlsx("研究设计-20260823.xlsx", sheet = "Conditions", skip = 1) %>%
  setNames(c("abbreviation_coding","full_names","abbrevation_writing","N_total",
             "N_male","N_female","excluded_male","excluded_female","label","category"))
data_effect0 <- read.csv(paste0("results/pairs-effect-", sex, "-exclude",
                                ifelse(exclude_diagnose_at_baseline, 1, 0), ".csv"))

###### selecting significant edge with HR > 1----------
data_effect <- data_effect0 %>%
  filter(coef_type == coef_type_target | is.na(coef_type)) %>%
  filter(n_event >= n_at_least & !is.na(coef_type)) %>%
  rename(p_effect = Pr...z.., HR = exp.coef.) %>%
  filter(p_effect < 0.05/nrow(.)) %>%
  filter(HR > 1) %>%
  filter(first_D != "death" & second_D != "death")


## Sum of bidirectional edges
edges <- data_effect %>%
  select(first_D, second_D, coef) %>%
  rowwise() %>%
  mutate(node1 = min(first_D, second_D), node2 = max(first_D, second_D)) %>%
  ungroup() %>%
  group_by(node1, node2) %>%
  summarise(weight = sum(coef), .groups = "drop")


## fitting bigclam
res_fix <- bigclam(edges, c = final_c, n_starts = n_starts, scale_w = scale_w)
print(res_fix)
cl_tab <- extract_clusters(res_fix,row_rel = 0,column_rel = 0,maxmean_rel = 0.5,drop_max = FALSE)
write.csv(cl_tab,file = paste0("results/cluster-bigclam-", sex,"-k",final_c, ".csv"))
p_fv <- plot_fvalues(x = res_fix$F,clusters = cl_tab,save = paste0("results/cluster-bigclam-", sex,"-k",final_c, ".tiff"))












