library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- here("final-figures", "figure-5")
save_dir <- file.path(analysis_dir, "rdata_files")
input_dir <- file.path(analysis_dir,"input_files")
plot_dir <- file.path(analysis_dir, "plots")
output_dir <- file.path(analysis_dir, "output_files")
gtac_sfpq_rci_expr_df_list <- readRDS(file.path(save_dir, "06-03-save-gtac-blood-pbmc-sfpq-affected-rci-expression-matrix.rds"))
model_validation_list <- readRDS(file.path(save_dir, "06-03-save-all-ridge-coefficient-models-from-grima-et-al.rds"))
selected_drosha_repeat_df <- readRDS("30-sequence-analysis-of-introns/rdata_files/10-01-hexanucleotide-repeats-with-increased-enrichment-for-drosha.rds")
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
# source("../../useful-functions/plot-utils.R")
# source("../../useful-functions/basepair_coverage_functions.R")
# source("../src/colours.R")
require("HiCBricks")
require("GenomicRanges")
require("stringr")
require("ggplot2")
library(rlang)
library(riskRegression)
library(survival)
library(survminer)
library(transformr)
library(plotly)
library(reshape2)
library(knitr)
library(rmarkdown)
require("RColorBrewer")
require("randomForestSRC")
require("parallel")
require("matrixStats")
require("survival")
require("patchwork")
library(broom)
# ==========================================================================
# Function
# ==========================================================================

prepare_progression_rsf_df <- function(df){
    df$time <- df$disease_duration
    df$status <- 1
    df$status[is.na(df$disease_duration)] <- 0
    return(df)
}

# ==========================================================================
# Load data
# ==========================================================================

patient_df <- gtac_sfpq_rci_expr_df_list$patient_df
sfpq_rci_expr_matrix <- gtac_sfpq_rci_expr_df_list$sfpq_affected_rcis_matrix
sfpq_affected_rcis <- gtac_sfpq_rci_expr_df_list$sfpq_affected_rcis
sfpq_affected_rloop_rcis <- gtac_sfpq_rci_expr_df_list$sfpq_affected_rloop_containing_rcis
drosha_responsive_rcis <- unique(selected_drosha_repeat_df$hash_id)

# ==========================================================================
# Create surv data
# ==========================================================================

patient_surv_df <- patient_df
patient_surv_df$disease_duration <- as.numeric(patient_surv_df$DZDURRNADRAW)
patient_surv_df$age_at_onset <- as.numeric(patient_surv_df$AGE_AT_ONSET)

# ==========================================================================
# Prepare RSF df
# ==========================================================================

non_na_patient_surv_df <- patient_surv_df[!is.na(patient_surv_df$age_at_onset),]
als_rsf_surv_df <- prepare_progression_rsf_df(non_na_patient_surv_df)

# ==========================================================================
# Get the models
# ==========================================================================

select_proteins <- c("DROSHA hexanucleotide", "CPEB4 trinucleotide", "CPSF6 trinucleotide")
selected_fits <- model_validation_list$rbp_fits[select_proteins]
selected_rbp_enriched_rcis <- model_validation_list$selected_rbp_enriched_rcis[select_proteins]

# ==========================================================================
# Create disease duration data
# ==========================================================================
# Test DROSHA 

base_fit <- coxph(Surv(time, status) ~ age_at_onset + sex, data = als_rsf_surv_df, x = TRUE)
summary(base_fit)
# Call:
# coxph(formula = Surv(time, status) ~ age_at_onset + sex, data = als_rsf_surv_df,
#     x = TRUE)

#   n= 226, number of events= 226

#                   coef exp(coef)  se(coef)      z Pr(>|z|)
# age_at_onset  0.031128  1.031618  0.006124  5.083 3.71e-07 ***
# sexmale      -0.184565  0.831466  0.140034 -1.318    0.188
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#              exp(coef) exp(-coef) lower .95 upper .95
# age_at_onset    1.0316     0.9694    1.0193     1.044
# sexmale         0.8315     1.2027    0.6319     1.094

# Concordance= 0.603  (se = 0.022 )
# Likelihood ratio test= 28.96  on 2 df,   p=5e-07
# Wald test            = 27.04  on 2 df,   p=1e-06
# Score (logrank) test = 26.81  on 2 df,   p=2e-06


require(glmnet)

use_rci_matrix_t <- t(sfpq_rci_expr_matrix[,als_rsf_surv_df$column_id])
X <- use_rci_matrix_t
y <- Surv(als_rsf_surv_df$time, als_rsf_surv_df$status)


use_hash_ids <- rownames(coef(selected_fits[["DROSHA hexanucleotide"]], s = "lambda.min"))
cur_model <- selected_fits[["DROSHA hexanucleotide"]]
cur_data <- X[,use_hash_ids]

als_rsf_surv_df$drosha_score <- as.numeric(predict(cur_model, newx = cur_data, s="lambda.min", type="link"))
als_rsf_surv_df$drosha_score_scaled <- scale(als_rsf_surv_df$drosha_score)

drosha_cox_fit <- coxph(Surv(time, status) ~ drosha_score_scaled, data = als_rsf_surv_df, x = TRUE)
summary(drosha_cox_fit)
# Call:
# coxph(formula = Surv(time, status) ~ drosha_score_scaled, data = als_rsf_surv_df,
#     x = TRUE)

#   n= 226, number of events= 226

#                        coef exp(coef) se(coef)     z Pr(>|z|)
# drosha_score_scaled 0.18551   1.20383  0.06369 2.913  0.00358 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#                     exp(coef) exp(-coef) lower .95 upper .95
# drosha_score_scaled     1.204     0.8307     1.063     1.364

# Concordance= 0.559  (se = 0.023 )
# Likelihood ratio test= 8.41  on 1 df,   p=0.004
# Wald test            = 8.48  on 1 df,   p=0.004
# Score (logrank) test = 8.38  on 1 df,   p=0.004


full_drosha_fit <- coxph(Surv(time, status) ~ age_at_onset + sex + drosha_score_scaled, data = als_rsf_surv_df, x = TRUE)
summary(full_drosha_fit)
# Call:
# coxph(formula = Surv(time, status) ~ age_at_onset + sex + drosha_score_scaled,
#     data = als_rsf_surv_df, x = TRUE)

#   n= 226, number of events= 226

#                          coef exp(coef)  se(coef)      z Pr(>|z|)
# age_at_onset         0.027797  1.028187  0.006331  4.391 1.13e-05 ***
# sexmale             -0.242960  0.784303  0.143647 -1.691   0.0908 .
# drosha_score_scaled  0.123681  1.131655  0.067175  1.841   0.0656 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#                     exp(coef) exp(-coef) lower .95 upper .95
# age_at_onset           1.0282     0.9726    1.0155     1.041
# sexmale                0.7843     1.2750    0.5918     1.039
# drosha_score_scaled    1.1317     0.8837    0.9921     1.291

# Concordance= 0.608  (se = 0.022 )
# Likelihood ratio test= 32.32  on 3 df,   p=4e-07
# Wald test            = 31.08  on 3 df,   p=8e-07
# Score (logrank) test = 31.12  on 3 df,   p=8e-07


# Test CPSF6 
use_hash_ids <- rownames(coef(selected_fits[["CPSF6 trinucleotide"]], s = "lambda.min"))
cur_model <- selected_fits[["CPSF6 trinucleotide"]]
cur_data <- X[,use_hash_ids]

als_rsf_surv_df$cpsf6_score <- as.numeric(predict(cur_model, newx = cur_data, s="lambda.min", type="link"))
als_rsf_surv_df$cpsf6_score_scaled <- scale(als_rsf_surv_df$cpsf6_score)

full_cpsf6_fit <- coxph(Surv(time, status) ~ age_at_onset + sex + cpsf6_score_scaled, data = als_rsf_surv_df, x = TRUE)
summary(full_cpsf6_fit)
# Call:
# coxph(formula = Surv(time, status) ~ age_at_onset + sex + cpsf6_score_scaled,
#     data = als_rsf_surv_df, x = TRUE)

#   n= 226, number of events= 226

#                         coef exp(coef)  se(coef)      z Pr(>|z|)
# age_at_onset        0.032810  1.033354  0.006148  5.337 9.46e-08 ***
# sexmale            -0.175294  0.839210  0.140081 -1.251   0.2108
# cpsf6_score_scaled  0.123711  1.131689  0.064782  1.910   0.0562 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#                    exp(coef) exp(-coef) lower .95 upper .95
# age_at_onset          1.0334     0.9677    1.0210     1.046
# sexmale               0.8392     1.1916    0.6377     1.104
# cpsf6_score_scaled    1.1317     0.8836    0.9967     1.285

# Concordance= 0.607  (se = 0.022 )
# Likelihood ratio test= 32.86  on 3 df,   p=3e-07
# Wald test            = 31.37  on 3 df,   p=7e-07
# Score (logrank) test = 31.24  on 3 df,   p=8e-07


# Test CPEB4
use_hash_ids <- rownames(coef(selected_fits[["CPEB4 trinucleotide"]], s = "lambda.min"))
cur_model <- selected_fits[["CPEB4 trinucleotide"]]
cur_data <- X[,use_hash_ids]

als_rsf_surv_df$cpeb4_score <- as.numeric(predict(cur_model, newx = cur_data, s="lambda.min", type="link"))
als_rsf_surv_df$cpeb4_score_scaled <- scale(als_rsf_surv_df$cpeb4_score)

cpeb4_cox_fit <- coxph(Surv(time, status) ~ cpeb4_score_scaled, data = als_rsf_surv_df, x = TRUE)
summary(cpeb4_cox_fit)
# Call:
# coxph(formula = Surv(time, status) ~ cpeb4_score_scaled, data = als_rsf_surv_df,
#     x = TRUE)

#   n= 226, number of events= 226

#                       coef exp(coef) se(coef)     z Pr(>|z|)
# cpeb4_score_scaled 0.12269   1.13053  0.06004 2.043    0.041 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#                    exp(coef) exp(-coef) lower .95 upper .95
# cpeb4_score_scaled     1.131     0.8845     1.005     1.272

# Concordance= 0.532  (se = 0.021 )
# Likelihood ratio test= 4.37  on 1 df,   p=0.04
# Wald test            = 4.18  on 1 df,   p=0.04
# Score (logrank) test = 4.15  on 1 df,   p=0.04


full_cpeb4_fit <- coxph(Surv(time, status) ~ age_at_onset + sex + cpeb4_score_scaled, data = als_rsf_surv_df, x = TRUE)
summary(full_cpeb4_fit)
# Call:
# coxph(formula = Surv(time, status) ~ age_at_onset + sex + cpeb4_score_scaled,
#     data = als_rsf_surv_df, x = TRUE)

#   n= 226, number of events= 226

#                         coef exp(coef)  se(coef)      z Pr(>|z|)
# age_at_onset        0.030817  1.031297  0.006093  5.058 4.24e-07 ***
# sexmale            -0.229184  0.795182  0.140862 -1.627   0.1037
# cpeb4_score_scaled  0.142096  1.152687  0.063021  2.255   0.0241 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#                    exp(coef) exp(-coef) lower .95 upper .95
# age_at_onset          1.0313     0.9697    1.0191     1.044
# sexmale               0.7952     1.2576    0.6033     1.048
# cpeb4_score_scaled    1.1527     0.8675    1.0188     1.304

# Concordance= 0.612  (se = 0.022 )
# Likelihood ratio test= 34.23  on 3 df,   p=2e-07
# Wald test            = 32.28  on 3 df,   p=5e-07
# Score (logrank) test = 32.37  on 3 df,   p=4e-07


# Test base fit vs CPEB4
anova(base_fit, full_cpeb4_fit)
# Analysis of Deviance Table
#  Cox model: response is  Surv(time, status)
#  Model 1: ~ age_at_onset + sex
#  Model 2: ~ age_at_onset + sex + cpeb4_score_scaled
#    loglik Chisq Df Pr(>|Chi|)
# 1 -988.19
# 2 -985.51 5.366  1    0.02053 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# Test base fit vs DROSHA
anova(base_fit, full_drosha_fit)
# Analysis of Deviance Table
#  Cox model: response is  Surv(time, status)
#  Model 1: ~ age_at_onset + sex
#  Model 2: ~ age_at_onset + sex + drosha_score_scaled
#    loglik  Chisq Df Pr(>|Chi|)
# 1 -988.19
# 2 -986.28 3.8096  1    0.05096 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

saveRDS(object = als_rsf_surv_df, file = file.path(save_dir, "06-04-save-gtac-predicted-score.rds"))


