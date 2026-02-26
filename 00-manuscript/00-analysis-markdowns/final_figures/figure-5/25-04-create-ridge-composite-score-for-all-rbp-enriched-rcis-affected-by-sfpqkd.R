library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- here("final-figures", "figure-5")
save_dir <- file.path(analysis_dir, "rdata_files")
input_dir <- file.path(analysis_dir,"input_files")
plot_dir <- file.path(analysis_dir, "plots")
output_dir <- file.path(analysis_dir, "output_files")
grima_sfpq_rci_expr_df_list <- readRDS(file.path(save_dir, "06-01-save-grima-blood-pbmc-sfpq-affected-rci-expression-matrix.rds"))
selected_rbp_repeat_df <- readRDS("30-sequence-analysis-of-introns/rdata_files/10-02-all-repeats-with-increased-enrichment-for-top10-proteins.rds")
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
# source("../../useful-functions/plot-utils.R")
# source("../../useful-functions/basepair_coverage_functions.R")
# source("../src/colours.R")
source("48-test-SFPQ-affected-RCI-expression-in-survival-models/scripts/src/sources.R")
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
# How many RBPs have Rloops
# ==========================================================================
write.csv(selected_rbp_repeat_df, file = file.path(input_dir, "10-02-enrichment-of-top10-rbp-binding-sites-in-repeat-proximal-regions.csv"))
# selected_rbp_repeat_df$rloop <- selected_rbp_repeat_df$hash_id %in% sfpq_affected_rloop_rcis
# table(unique(selected_rbp_repeat_df[,c("hash_id", "rloop")])[,c("rloop")])
# FALSE  TRUE
#     4   810
# table(unique(selected_rbp_repeat_df[,c("rbp", "hash_id", "rloop")])[,c("rbp", "rloop")])
#          rloop
# rbp       FALSE TRUE
#   AATF        0   32
#   CPEB4       1  243
#   CPSF6       1  243
#   CSTF2       0   17
#   DDX21       0   17
#   DDX6        0  170
#   DROSHA      1  186
#   EIF3D       2  233
#   FAM120A     0  183
#   GARS        0   14
#   GEMIN5      1  152
#   GPKOW       3  220
#   GTF2F1      1  167
#   METTL1      3  283
#   MORC2       0   17
#   MTPAP       2  269
#   NIP7        1  187
#   NIPBL       0   32
#   NOL12       0  178
#   NPM1        3  261
#   PCBP1       0  201
#   PCBP2       0  194
#   PUM1        1  228
#   RBM22       0   27
#   RPS11       1  164
#   RPS6        0   32
#   SF3B1       0   27
#   STAU2       0   22
#   TROVE2      2  235
#   ZC3H8       0   16

# ==========================================================================
# Function
# ==========================================================================

prepare_progression_rsf_df <- function(df){
    df$time <- df$disease_duration
    df$status <- 1
    df$status[df$Deceased != 1] <- 0
    return(df)
}

# ==========================================================================
# Load data
# ==========================================================================

patient_df <- grima_sfpq_rci_expr_df_list$patient_df
sfpq_rci_expr_matrix <- grima_sfpq_rci_expr_df_list$sfpq_affected_rcis_matrix
sfpq_affected_rcis <- grima_sfpq_rci_expr_df_list$sfpq_affected_rcis
sfpq_affected_rloop_rcis <- grima_sfpq_rci_expr_df_list$sfpq_affected_rloop_containing_rcis


# ==========================================================================
# Create surv data
# ==========================================================================

patient_surv_df <- patient_df
patient_surv_df$disease_duration <- as.numeric(patient_surv_df$Disease_duration_months)
patient_surv_df$Disease_status <- factor(patient_surv_df$Disease_status, levels = c("Control", "sALS"))
patient_surv_df$age_at_collection <- as.numeric(patient_surv_df$Age_at_collection)
patient_surv_df$age_at_onset <- as.numeric(patient_surv_df$Age_of_onset)

# ==========================================================================
# Prepare RSF df
# ==========================================================================

non_na_patient_surv_df <- patient_surv_df[!is.na(patient_surv_df$age_at_collection),]
als_surv_df <- non_na_patient_surv_df[non_na_patient_surv_df$Disease_status == "sALS",]
als_rsf_surv_df <- prepare_progression_rsf_df(als_surv_df)

# ==========================================================================
# Create disease duration data
# ==========================================================================

# remove rcis with too many NAs

abundant_expr_rci_matrix <- sfpq_rci_expr_matrix[rowSums(sfpq_rci_expr_matrix > 0) > (ncol(sfpq_rci_expr_matrix) * 0.8),]

selected_sfpq_affected_rcis <- sfpq_affected_rcis[sfpq_affected_rcis %in% rownames(abundant_expr_rci_matrix)]
selected_sfpq_affected_rloop_rcis <- sfpq_affected_rloop_rcis[sfpq_affected_rloop_rcis %in% rownames(abundant_expr_rci_matrix)]

hash_id_freq <- table(unique(selected_rbp_repeat_df[,c("hash_id", "repeat_class")])$hash_id)
rbp_hash_id_freq_df <- data.frame(hash_id = names(hash_id_freq),
	frequency = as.vector(hash_id_freq))
rownames(rbp_hash_id_freq_df) <- NULL

selected_rbp_enriched_rcis <- lapply(split(selected_rbp_repeat_df, paste(selected_rbp_repeat_df$rbp, selected_rbp_repeat_df$repeat_class)), function(x){
	unique(x$hash_id[x$hash_id %in% rownames(abundant_expr_rci_matrix) & 
        x$hash_id %in% rbp_hash_id_freq_df$hash_id[rbp_hash_id_freq_df$frequency == 1]])
})

selected_rbp_enriched_rloop_containing_rcis <- lapply(split(selected_rbp_repeat_df, paste(selected_rbp_repeat_df$rbp, selected_rbp_repeat_df$repeat_class)), function(x){
	x$hash_id[x$hash_id %in% rownames(abundant_expr_rci_matrix) & x$hash_id %in% selected_sfpq_affected_rloop_rcis]
})

# range(vapply(selected_rbp_enriched_rloop_containing_rcis, length, 1) / vapply(selected_rbp_enriched_rcis, length, 1))
# [1] 1 1

selected_rbp_enriched_rcis <- selected_rbp_enriched_rcis[vapply(selected_rbp_enriched_rcis, length, 1) > 3]

all_selected_rcis <- unique(do.call(c, selected_rbp_enriched_rcis))
abundant_single_repeat_class_selected_rbp_repeat_df <- selected_rbp_repeat_df[selected_rbp_repeat_df$hash_id %in% all_selected_rcis,]
write.csv(abundant_single_repeat_class_selected_rbp_repeat_df, file = file.path(input_dir, "10-03-selected-rcis-with-enrichment-of-top10-rbp-binding-sites-in-repeat-proximal-regions.csv"))

# disease_duration_model_outputs <- do.call(rbind, lapply(selected_sfpq_affected_rloop_rcis, function(x){
#     als_rsf_surv_df$scaled_var <- scale(abundant_expr_rci_matrix[x,als_rsf_surv_df$column_id])
#     model <- coxph(Surv(time, status) ~ scaled_var, 
#         data=als_rsf_surv_df, x=TRUE)
#     return_df <- as.data.frame(tidy(model, exponentiate = TRUE, conf.int = TRUE))
#     cox_proportionality <- cox.zph(model)
#     return_df$cox_zph <- cox_proportionality$table[1,3]
#     return_df$term[return_df$term == "scaled_var"] <- x
#     return(return_df)
# }))
# disease_duration_model_outputs[disease_duration_model_outputs$p.value < 0.05,]


# Construct elastic net score
require(glmnet)
rci_matrix_t <- t(abundant_expr_rci_matrix[,als_rsf_surv_df$column_id])
X <- rci_matrix_t
y <- Surv(als_rsf_surv_df$time, als_rsf_surv_df$status)


all_models <- lapply(selected_rbp_enriched_rcis, function(x){
	set.seed(123)
	foldid <- sample(rep(1:10, length.out=nrow(X)))
	X2 <- X[,x]
	cv.glmnet(X2, y, family="cox", alpha=0, nfolds = 10, foldid=foldid)
})



scores <- lapply(names(all_models), function(x){
    cur_data <- X[,selected_rbp_enriched_rcis[[x]]]
    as.numeric(predict(all_models[[x]], newx = cur_data, s="lambda.min", type="link"))
})
names(scores) <- names(all_models)

#
library(rms)
fit_baseline <- cph(Surv(time,status) ~ age_at_onset + Sex,
       data=als_rsf_surv_df,
       x=TRUE, y=TRUE, surv=TRUE)
validate(fit_baseline, method="boot", B=300)

#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3714   0.3773 0.3638   0.0135          0.3579  0.2274 0.4987 300
# R2        0.2514   0.2708 0.2408   0.0300          0.2214  0.0659 0.3768 300
# Slope     1.0000   1.0000 0.9415   0.0585          0.9415  0.4853 1.6190 300
# D         0.0423   0.0473 0.0402   0.0071          0.0352  0.0005 0.0660 300
# U        -0.0032  -0.0032 0.0030  -0.0062          0.0030 -0.0027 0.0241 300
# Q         0.0455   0.0505 0.0372   0.0133          0.0322 -0.0197 0.0634 300
# g         0.8000   0.8532 0.7710   0.0821          0.7179  0.3084 1.0681 300


all_cox_validations <- lapply(scores, function(x){
    als_rsf_surv_df$score <- x
    fit <- cph(Surv(time,status) ~ age_at_onset + Sex + score,
           data=als_rsf_surv_df,
           x=TRUE, y=TRUE, surv=TRUE)
    validate(fit, method="boot", B=300)
})
# $`AATF quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4192   0.4276 0.4086   0.0190          0.4002  0.2706 0.5324 300
# R2        0.3112   0.3360 0.2956   0.0404          0.2708  0.1161 0.4204 300
# Slope     1.0000   1.0000 0.8937   0.1063          0.8937  0.4633 1.3562 300
# D         0.0549   0.0617 0.0516   0.0101          0.0448  0.0056 0.0782 300
# U        -0.0032  -0.0032 0.0054  -0.0085          0.0054 -0.0051 0.0578 300
# Q         0.0581   0.0649 0.0462   0.0187          0.0394 -0.0364 0.0795 300
# g         0.8509   0.9412 0.8242   0.1169          0.7340  0.2860 1.0770 300

# $`CPEB4 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.5580   0.5610 0.5502   0.0108          0.5472  0.4602 0.6372 300
# R2        0.5294   0.5352 0.5175   0.0177          0.5117  0.4021 0.6368 300
# Slope     1.0000   1.0000 0.9438   0.0562          0.9438  0.6669 1.2492 300
# D         0.1126   0.1161 0.1089   0.0072          0.1054  0.0643 0.1435 300
# U        -0.0032  -0.0032 0.0020  -0.0052          0.0020 -0.0019 0.0150 300
# Q         0.1158   0.1192 0.1068   0.0124          0.1034  0.0544 0.1428 300
# g         1.5024   1.5506 1.4578   0.0928          1.4096  0.9658 1.8166 300

# $`CPSF6 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.5487   0.5552 0.5423   0.0129          0.5358  0.4363 0.6411 300
# R2        0.5359   0.5450 0.5234   0.0216          0.5143  0.3896 0.6449 300
# Slope     1.0000   1.0000 0.9357   0.0643          0.9357  0.6382 1.2436 300
# D         0.1147   0.1192 0.1107   0.0084          0.1063  0.0586 0.1493 300
# U        -0.0032  -0.0032 0.0026  -0.0058          0.0026 -0.0025 0.0219 300
# Q         0.1179   0.1223 0.1081   0.0142          0.1037  0.0416 0.1486 300
# g         1.4900   1.5485 1.4448   0.1036          1.3864  0.8375 1.8222 300

# $`DDX6 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4863   0.4923 0.4765   0.0158          0.4705  0.3510 0.5918 300
# R2        0.4272   0.4439 0.4145   0.0294          0.3978  0.2554 0.5420 300
# Slope     1.0000   1.0000 0.9308   0.0692          0.9308  0.6228 1.2928 300
# D         0.0829   0.0889 0.0796   0.0093          0.0735  0.0296 0.1105 300
# U        -0.0032  -0.0032 0.0026  -0.0058          0.0026 -0.0023 0.0295 300
# Q         0.0860   0.0921 0.0769   0.0151          0.0709  0.0114 0.1103 300
# g         1.1869   1.2510 1.1526   0.0984          1.0885  0.6546 1.4624 300

# $`DROSHA hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4937   0.5027 0.4805   0.0223          0.4714  0.3667 0.5793 300
# R2        0.4116   0.4347 0.3962   0.0385          0.3731  0.2314 0.5083 300
# Slope     1.0000   1.0000 0.9085   0.0915          0.9085  0.5614 1.2973 300
# D         0.0788   0.0861 0.0749   0.0112          0.0676  0.0245 0.1023 300
# U        -0.0032  -0.0032 0.0034  -0.0066          0.0034 -0.0031 0.0270 300
# Q         0.0819   0.0893 0.0715   0.0178          0.0642  0.0043 0.1026 300
# g         1.1363   1.2192 1.0933   0.1258          1.0105  0.5363 1.3739 300

# $`EIF3D trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4914   0.4960 0.4781   0.0180          0.4734  0.3505 0.6055 300
# R2        0.4319   0.4525 0.4187   0.0337          0.3981  0.2240 0.5557 300
# Slope     1.0000   1.0000 0.9195   0.0805          0.9195  0.5590 1.3593 300
# D         0.0841   0.0916 0.0807   0.0110          0.0731  0.0195 0.1166 300
# U        -0.0032  -0.0032 0.0037  -0.0069          0.0037 -0.0033 0.0296 300
# Q         0.0873   0.0948 0.0770   0.0178          0.0694 -0.0024 0.1167 300
# g         1.2277   1.3102 1.1865   0.1237          1.1039  0.5545 1.5646 300

# $`FAM120A hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4853   0.4852 0.4715   0.0137          0.4716  0.3642 0.5797 300
# R2        0.3958   0.4113 0.3797   0.0316          0.3642  0.2211 0.5020 300
# Slope     1.0000   1.0000 0.9301   0.0699          0.9301  0.5948 1.3442 300
# D         0.0748   0.0797 0.0708   0.0089          0.0659  0.0253 0.0993 300
# U        -0.0032  -0.0032 0.0026  -0.0058          0.0026 -0.0023 0.0195 300
# Q         0.0779   0.0829 0.0682   0.0147          0.0633  0.0114 0.0994 300
# g         1.1021   1.1553 1.0596   0.0957          1.0064  0.5816 1.3693 300

# $`GEMIN5 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.5049   0.5069 0.4912   0.0157          0.4892  0.3796 0.6013 300
# R2        0.4147   0.4362 0.4002   0.0360          0.3787  0.2192 0.5257 300
# Slope     1.0000   1.0000 0.9251   0.0749          0.9251  0.5040 1.3732 300
# D         0.0796   0.0871 0.0759   0.0111          0.0685  0.0181 0.1066 300
# U        -0.0032  -0.0032 0.0049  -0.0081          0.0049 -0.0048 0.0548 300
# Q         0.0828   0.0902 0.0710   0.0192          0.0636 -0.0224 0.1069 300
# g         1.0678   1.1663 1.0389   0.1274          0.9404  0.3281 1.3376 300

# $`GPKOW trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.5068   0.5038 0.4936   0.0103          0.4965  0.3897 0.6178 300
# R2        0.4497   0.4591 0.4360   0.0231          0.4266  0.2837 0.5751 300
# Slope     1.0000   1.0000 0.9368   0.0632          0.9368  0.5701 1.3494 300
# D         0.0889   0.0931 0.0852   0.0078          0.0811  0.0362 0.1206 300
# U        -0.0032  -0.0032 0.0033  -0.0065          0.0033 -0.0031 0.0301 300
# Q         0.0921   0.0962 0.0819   0.0143          0.0777  0.0156 0.1199 300
# g         1.2490   1.3095 1.2086   0.1009          1.1480  0.6294 1.5638 300

# $`GTF2F1 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4974   0.4995 0.4867   0.0128          0.4847  0.3646 0.6021 300
# R2        0.4217   0.4301 0.4072   0.0229          0.3987  0.2552 0.5414 300
# Slope     1.0000   1.0000 0.9438   0.0562          0.9438  0.6038 1.3671 300
# D         0.0814   0.0849 0.0777   0.0072          0.0742  0.0305 0.1111 300
# U        -0.0032  -0.0032 0.0027  -0.0058          0.0027 -0.0024 0.0251 300
# Q         0.0846   0.0881 0.0750   0.0130          0.0715  0.0125 0.1102 300
# g         1.1664   1.2166 1.1285   0.0881          1.0783  0.5854 1.4656 300

# $`METTL1 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4183   0.4352 0.4127   0.0224          0.3958  0.2853 0.5203 300
# R2        0.3056   0.3375 0.2916   0.0459          0.2597  0.1251 0.3990 300
# Slope     1.0000   1.0000 0.8989   0.1011          0.8989  0.5110 1.3698 300
# D         0.0537   0.0618 0.0507   0.0111          0.0426  0.0091 0.0736 300
# U        -0.0032  -0.0032 0.0036  -0.0068          0.0036 -0.0034 0.0290 300
# Q         0.0569   0.0650 0.0471   0.0179          0.0390 -0.0117 0.0741 300
# g         0.8814   0.9760 0.8493   0.1267          0.7547  0.3581 1.1047 300

# $`METTL1 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4811   0.4859 0.4717   0.0142          0.4669  0.3647 0.5842 300
# R2        0.4063   0.4217 0.3897   0.0320          0.3743  0.2290 0.5181 300
# Slope     1.0000   1.0000 0.9138   0.0862          0.9138  0.5368 1.3394 300
# D         0.0774   0.0829 0.0733   0.0096          0.0678  0.0259 0.1045 300
# U        -0.0032  -0.0032 0.0041  -0.0073          0.0041 -0.0040 0.0365 300
# Q         0.0806   0.0861 0.0692   0.0169          0.0637 -0.0028 0.1037 300
# g         1.0749   1.1460 1.0308   0.1151          0.9598  0.4944 1.3278 300

# $`MTPAP quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3940   0.4004 0.3839   0.0166          0.3775  0.2526 0.4989 300
# R2        0.2904   0.3107 0.2739   0.0369          0.2536  0.1067 0.3929 300
# Slope     1.0000   1.0000 0.9118   0.0882          0.9118  0.5310 1.3864 300
# D         0.0504   0.0557 0.0470   0.0087          0.0417  0.0063 0.0710 300
# U        -0.0032  -0.0032 0.0029  -0.0061          0.0029 -0.0029 0.0206 300
# Q         0.0536   0.0589 0.0440   0.0148          0.0387 -0.0102 0.0709 300
# g         0.8483   0.9105 0.8126   0.0979          0.7504  0.3904 1.0570 300

# $`MTPAP trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4779   0.4900 0.4737   0.0163          0.4616  0.3582 0.5639 300
# R2        0.4242   0.4390 0.4108   0.0282          0.3959  0.2498 0.5363 300
# Slope     1.0000   1.0000 0.9407   0.0593          0.9407  0.5599 1.3773 300
# D         0.0820   0.0874 0.0786   0.0088          0.0732  0.0283 0.1099 300
# U        -0.0032  -0.0032 0.0031  -0.0063          0.0032 -0.0030 0.0260 300
# Q         0.0852   0.0906 0.0755   0.0151          0.0701  0.0087 0.1092 300
# g         1.1971   1.2591 1.1588   0.1004          1.0967  0.5726 1.5186 300

# $`NIP7 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4835   0.4909 0.4745   0.0164          0.4671  0.3562 0.5905 300
# R2        0.4077   0.4280 0.3938   0.0342          0.3735  0.2344 0.5329 300
# Slope     1.0000   1.0000 0.9166   0.0834          0.9166  0.5790 1.3205 300
# D         0.0778   0.0844 0.0743   0.0101          0.0677  0.0274 0.1085 300
# U        -0.0032  -0.0032 0.0031  -0.0063          0.0031 -0.0030 0.0209 300
# Q         0.0809   0.0876 0.0712   0.0164          0.0645  0.0107 0.1074 300
# g         1.1545   1.2295 1.1158   0.1137          1.0408  0.5765 1.4525 300

# $`NIPBL quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4280   0.4280 0.4117   0.0164          0.4117  0.2892 0.5342 300
# R2        0.3109   0.3307 0.2961   0.0346          0.2763  0.1237 0.4133 300
# Slope     1.0000   1.0000 0.9073   0.0927          0.9073  0.5254 1.3564 300
# D         0.0549   0.0603 0.0517   0.0087          0.0462  0.0077 0.0761 300
# U        -0.0032  -0.0032 0.0036  -0.0068          0.0036 -0.0035 0.0348 300
# Q         0.0580   0.0635 0.0481   0.0154          0.0426 -0.0204 0.0755 300
# g         0.8732   0.9441 0.8448   0.0993          0.7738  0.3586 1.0691 300

# $`NOL12 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4737   0.4793 0.4644   0.0148          0.4588  0.3413 0.5811 300
# R2        0.3931   0.4088 0.3790   0.0298          0.3633  0.2147 0.5043 300
# Slope     1.0000   1.0000 0.9173   0.0827          0.9173  0.5804 1.3276 300
# D         0.0741   0.0795 0.0707   0.0088          0.0653  0.0232 0.1008 300
# U        -0.0032  -0.0032 0.0030  -0.0062          0.0030 -0.0028 0.0282 300
# Q         0.0772   0.0827 0.0676   0.0150          0.0622  0.0023 0.1002 300
# g         1.1277   1.2007 1.0872   0.1135          1.0142  0.5428 1.4045 300

# $`NPM1 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4178   0.4106 0.4016   0.0090          0.4088  0.2813 0.5449 300
# R2        0.2955   0.3077 0.2811   0.0266          0.2689  0.1137 0.4095 300
# Slope     1.0000   1.0000 0.9290   0.0710          0.9290  0.5247 1.4231 300
# D         0.0515   0.0552 0.0485   0.0068          0.0447  0.0070 0.0745 300
# U        -0.0032  -0.0032 0.0032  -0.0064          0.0032 -0.0032 0.0280 300
# Q         0.0547   0.0584 0.0452   0.0132          0.0415 -0.0150 0.0738 300
# g         0.8505   0.9019 0.8181   0.0837          0.7668  0.3438 1.0690 300

# $`NPM1 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4783   0.4857 0.4724   0.0132          0.4651  0.3490 0.5838 300
# R2        0.4319   0.4442 0.4154   0.0289          0.4030  0.2426 0.5595 300
# Slope     1.0000   1.0000 0.9285   0.0715          0.9285  0.5731 1.3308 300
# D         0.0841   0.0891 0.0798   0.0093          0.0748  0.0251 0.1162 300
# U        -0.0032  -0.0032 0.0034  -0.0066          0.0034 -0.0030 0.0372 300
# Q         0.0873   0.0922 0.0764   0.0158          0.0714 -0.0026 0.1148 300
# g         1.2249   1.2767 1.1709   0.1058          1.1191  0.5884 1.5434 300

# $`PCBP1 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4611   0.4709 0.4527   0.0182          0.4429  0.3297 0.5713 300
# R2        0.3860   0.4028 0.3717   0.0311          0.3549  0.2082 0.5074 300
# Slope     1.0000   1.0000 0.9308   0.0692          0.9308  0.6155 1.3684 300
# D         0.0723   0.0782 0.0689   0.0093          0.0630  0.0205 0.1000 300
# U        -0.0032  -0.0032 0.0025  -0.0057          0.0025 -0.0022 0.0176 300
# Q         0.0755   0.0814 0.0664   0.0150          0.0605  0.0063 0.0983 300
# g         1.0813   1.1470 1.0448   0.1022          0.9791  0.5427 1.3686 300

# $`PCBP1 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3978   0.4087 0.3877   0.0210          0.3768  0.2615 0.5001 300
# R2        0.2900   0.3140 0.2739   0.0401          0.2499  0.0999 0.3922 300
# Slope     1.0000   1.0000 0.9014   0.0986          0.9014  0.5224 1.3544 300
# D         0.0503   0.0566 0.0470   0.0096          0.0407  0.0035 0.0707 300
# U        -0.0032  -0.0032 0.0032  -0.0064          0.0032 -0.0029 0.0316 300
# Q         0.0535   0.0598 0.0437   0.0160          0.0374 -0.0209 0.0705 300
# g         0.8540   0.9249 0.8173   0.1076          0.7464  0.3355 1.0685 300

# $`PCBP2 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4849   0.4964 0.4767   0.0197          0.4652  0.3517 0.5881 300
# R2        0.4132   0.4347 0.3996   0.0351          0.3782  0.2295 0.5270 300
# Slope     1.0000   1.0000 0.9187   0.0813          0.9187  0.5632 1.3129 300
# D         0.0792   0.0865 0.0758   0.0107          0.0685  0.0228 0.1079 300
# U        -0.0032  -0.0032 0.0035  -0.0066          0.0035 -0.0033 0.0288 300
# Q         0.0824   0.0896 0.0723   0.0173          0.0650  0.0008 0.1081 300
# g         1.1584   1.2417 1.1208   0.1209          1.0375  0.5189 1.4490 300

# $`PCBP2 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3926   0.4002 0.3779   0.0224          0.3703  0.2599 0.4924 300
# R2        0.2806   0.3038 0.2613   0.0425          0.2381  0.0991 0.3731 300
# Slope     1.0000   1.0000 0.8920   0.1080          0.8920  0.5018 1.3391 300
# D         0.0483   0.0541 0.0444   0.0097          0.0386  0.0055 0.0664 300
# U        -0.0032  -0.0032 0.0035  -0.0067          0.0035 -0.0031 0.0322 300
# Q         0.0515   0.0573 0.0409   0.0164          0.0351 -0.0183 0.0668 300
# g         0.8469   0.9097 0.7955   0.1142          0.7327  0.3572 1.0373 300

# $`PUM1 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.5016   0.5063 0.4951   0.0112          0.4904  0.3865 0.6066 300
# R2        0.4768   0.4867 0.4638   0.0229          0.4540  0.3139 0.5879 300
# Slope     1.0000   1.0000 0.9330   0.0670          0.9330  0.6249 1.2930 300
# D         0.0966   0.1010 0.0929   0.0081          0.0884  0.0414 0.1281 300
# U        -0.0032  -0.0032 0.0025  -0.0057          0.0025 -0.0024 0.0197 300
# Q         0.0997   0.1042 0.0904   0.0138          0.0859  0.0273 0.1280 300
# g         1.3617   1.4186 1.3168   0.1018          1.2599  0.7383 1.6875 300

# $`RBM22 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4145   0.4167 0.3939   0.0229          0.3917  0.2705 0.5144 300
# R2        0.2906   0.3194 0.2730   0.0464          0.2442  0.0874 0.3843 300
# Slope     1.0000   1.0000 0.8890   0.1110          0.8890  0.4657 1.3744 300
# D         0.0505   0.0577 0.0468   0.0110          0.0395  0.0012 0.0692 300
# U        -0.0032  -0.0032 0.0047  -0.0079          0.0047 -0.0046 0.0390 300
# Q         0.0536   0.0609 0.0421   0.0188          0.0348 -0.0262 0.0716 300
# g         0.8325   0.9230 0.7970   0.1260          0.7065  0.2773 1.0289 300

# $`RPS11 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4690   0.4743 0.4510   0.0233          0.4457  0.3355 0.5560 300
# R2        0.3610   0.3931 0.3457   0.0474          0.3136  0.1567 0.4498 300
# Slope     1.0000   1.0000 0.8983   0.1017          0.8983  0.4931 1.3243 300
# D         0.0663   0.0752 0.0627   0.0125          0.0538  0.0115 0.0863 300
# U        -0.0032  -0.0032 0.0045  -0.0077          0.0045 -0.0044 0.0414 300
# Q         0.0694   0.0784 0.0583   0.0201          0.0493 -0.0173 0.0882 300
# g         0.9810   1.0901 0.9470   0.1431          0.8379  0.3570 1.2008 300

# $`RPS6 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.4159   0.4239 0.4036   0.0203          0.3956  0.2660 0.5191 300
# R2        0.2989   0.3239 0.2829   0.0410          0.2578  0.0968 0.3965 300
# Slope     1.0000   1.0000 0.9119   0.0881          0.9119  0.5170 1.3388 300
# D         0.0522   0.0588 0.0488   0.0100          0.0422  0.0028 0.0731 300
# U        -0.0032  -0.0032 0.0032  -0.0064          0.0032 -0.0031 0.0346 300
# Q         0.0554   0.0620 0.0456   0.0163          0.0390 -0.0226 0.0738 300
# g         0.8951   0.9595 0.8555   0.1039          0.7912  0.3481 1.1039 300

# $`SF3B1 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3987   0.4167 0.3916   0.0252          0.3735  0.2489 0.4902 300
# R2        0.2918   0.3203 0.2746   0.0457          0.2461  0.0975 0.3862 300
# Slope     1.0000   1.0000 0.8922   0.1078          0.8922  0.5176 1.3555 300
# D         0.0507   0.0579 0.0471   0.0108          0.0399  0.0036 0.0698 300
# U        -0.0032  -0.0032 0.0034  -0.0066          0.0034 -0.0032 0.0309 300
# Q         0.0539   0.0611 0.0437   0.0174          0.0365 -0.0198 0.0702 300
# g         0.8866   0.9633 0.8422   0.1211          0.7656  0.3695 1.0925 300

# $`STAU2 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3871   0.3977 0.3781   0.0196          0.3675  0.2506 0.5025 300
# R2        0.2883   0.3110 0.2717   0.0392          0.2491  0.1000 0.4124 300
# Slope     1.0000   1.0000 0.9167   0.0833          0.9167  0.5155 1.5174 300
# D         0.0500   0.0561 0.0465   0.0095          0.0404  0.0049 0.0745 300
# U        -0.0032  -0.0032 0.0033  -0.0065          0.0033 -0.0031 0.0238 300
# Q         0.0531   0.0592 0.0432   0.0160          0.0371 -0.0140 0.0731 300
# g         0.8672   0.9306 0.8234   0.1072          0.7599  0.3534 1.1250 300

# $`TROVE2 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.5170   0.5174 0.5041   0.0133          0.5037  0.3930 0.6189 300
# R2        0.4574   0.4708 0.4459   0.0249          0.4326  0.2799 0.5773 300
# Slope     1.0000   1.0000 0.9334   0.0666          0.9334  0.5708 1.3376 300
# D         0.0911   0.0967 0.0879   0.0088          0.0823  0.0335 0.1233 300
# U        -0.0032  -0.0032 0.0036  -0.0067          0.0036 -0.0033 0.0309 300
# Q         0.0942   0.0998 0.0843   0.0155          0.0787  0.0125 0.1232 300
# g         1.2438   1.3145 1.2111   0.1034          1.1404  0.6086 1.5616 300

baseline_value_df <- extractC(validate(fit_baseline, method="boot", B=300))
baseline_value_df$rbp <- "baseline"
baseline_value_df$repeat_class <- "none"

all_rbp_cindex_df <- do.call(rbind, lapply(names(all_cox_validations), function(x){
    group_split <- str_split(x, pattern = " ", simplify = T)
    cindex_df <- extractC(all_cox_validations[[x]])
    cindex_df$rbp <- group_split[1]
    cindex_df$repeat_class <- group_split[2]
    return(cindex_df)
}))

all_cindex_df <- rbind(baseline_value_df, all_rbp_cindex_df)


all_results_df <- do.call(rbind, lapply(names(scores), function(x){
	als_rsf_surv_df$ridge_score <- scale(scores[[x]])

	fit_baseline <- coxph(Surv(time, status) ~ age_at_onset + Sex,
	    data = als_rsf_surv_df, x = TRUE)
	fit_ridge <- coxph(Surv(time, status) ~ age_at_onset + Sex + ridge_score,
	    data = als_rsf_surv_df, x = TRUE)
	label_split <- str_split(x, " ", simplify = T)
	return_df <- extract_cox_report(fit_baseline, 
		fit_ridge, 
		term = "ridge_score",
		label = x)
	return_df$rbp <- label_split[1]
	return_df$repeat_class <- label_split[2]
	return(return_df)
}))
# als_rsf_surv_df$ridge_score <- scores$ridge

return_list <- list(coxph_results = all_results_df,
	cindex_list = all_cindex_df,
	rbp_fits = all_models,
    all_rbp_scores = scores, 
    surv_df = als_rsf_surv_df,
	selected_rbp_enriched_rcis = selected_rbp_enriched_rcis)


saveRDS(object = return_list, file = file.path(save_dir, "06-03-save-all-ridge-coefficient-models-from-grima-et-al.rds"))

