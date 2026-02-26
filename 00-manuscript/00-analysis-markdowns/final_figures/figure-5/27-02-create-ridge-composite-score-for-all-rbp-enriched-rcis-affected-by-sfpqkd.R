library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- here("final-figures", "figure-5")
save_dir <- file.path(analysis_dir, "rdata_files")
input_dir <- file.path(analysis_dir,"input_files")
plot_dir <- file.path(analysis_dir, "plots")
output_dir <- file.path(analysis_dir, "output_files")
answerals_sfpq_rci_expr_df_list <- readRDS(file.path(save_dir, "07-01-save-answerals-iPSCs-sfpq-affected-rci-expression-matrix.rds"))
selected_rbp_repeat_df <- readRDS("30-sequence-analysis-of-introns/rdata_files/10-02-all-repeats-with-increased-enrichment-for-top10-proteins.rds")
model_validation_list <- readRDS(file.path(save_dir, "06-03-save-all-ridge-coefficient-models-from-grima-et-al.rds"))
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
# write.csv(selected_rbp_repeat_df, file = file.path(input_dir, "10-02-enrichment-of-top10-rbp-binding-sites-in-repeat-proximal-regions.csv"))
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
    df$status[!df$subject_died] <- 0
    return(df)
}

# ==========================================================================
# Load data
# ==========================================================================

patient_df <- answerals_sfpq_rci_expr_df_list$patient_df
sfpq_rci_expr_matrix <- answerals_sfpq_rci_expr_df_list$sfpq_affected_rcis_matrix
sfpq_affected_rcis <- answerals_sfpq_rci_expr_df_list$sfpq_affected_rcis
sfpq_affected_rloop_rcis <- answerals_sfpq_rci_expr_df_list$sfpq_affected_rloop_containing_rcis

# ==========================================================================
# Create surv data
# ==========================================================================

patient_surv_df <- patient_df

# ==========================================================================
# Prepare RSF df
# ==========================================================================

als_surv_df <- patient_surv_df[
    patient_surv_df$mnd_group == "ALS" & 
    !is.na(patient_surv_df$nfl) & 
    patient_surv_df$disease_duration > 0 & 
    patient_surv_df$subject_died,]
als_rsf_surv_df <- prepare_progression_rsf_df(als_surv_df)

# ==========================================================================
# Create disease duration data
# ==========================================================================

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

selected_rbp_enriched_rcis <- selected_rbp_enriched_rcis[vapply(selected_rbp_enriched_rcis, length, 1) > 3]

all_selected_rcis <- unique(do.call(c, selected_rbp_enriched_rcis))
abundant_single_repeat_class_selected_rbp_repeat_df <- selected_rbp_repeat_df[selected_rbp_repeat_df$hash_id %in% all_selected_rcis,]

write.csv(abundant_single_repeat_class_selected_rbp_repeat_df, file = file.path(input_dir, "10-04-iPSC-MNs-selected-rcis-with-enrichment-of-top10-rbp-binding-sites-in-repeat-proximal-regions.csv"))



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


als_rsf_surv_df$sex <- factor(als_rsf_surv_df$sex, levels = c("Male", "Female"))

library(rms)
fit_baseline <- cph(Surv(time,status) ~ nfl,
       data=als_rsf_surv_df,
       x=TRUE, y=TRUE, surv=TRUE)
validate(fit_baseline, method="boot", B=300)

#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3034   0.3005 0.3034  -0.0029          0.3063  0.2208 0.3952 300
# R2        0.2459   0.2488 0.2459   0.0029          0.2431  0.1502 0.3359 300
# Slope     1.0000   1.0000 1.0009  -0.0009          1.0009  0.7249 1.3217 300
# D         0.0316   0.0323 0.0316   0.0007          0.0309  0.0162 0.0445 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0042 300
# Q         0.0327   0.0333 0.0310   0.0023          0.0304  0.0135 0.0435 300
# g         0.6876   0.6974 0.6876   0.0098          0.6778  0.4811 0.8585 300

cor_test_df <- do.call(rbind, lapply(names(scores), function(x){
    y1 <- str_split(x, pattern = " ", simplify = T)
    als_rsf_surv_df$score <- scores[[x]]
    test_res <- cor.test(als_rsf_surv_df$nfl, als_rsf_surv_df$score)
    data.frame(x1 = "nfl", rbp = y1[1], 
        repeat_class = y1[2], pval = test_res$p.value, 
        cor = test_res$estimate)
}))
cor_test_df <- cor_test_df[cor_test_df$pval < 0.05,]

#        x1     rbp    repeat_class         pval       cor
# cor1  nfl   CPEB4   trinucleotide 0.0094562039 0.1762131
# cor2  nfl   CPSF6   trinucleotide 0.0039703712 0.1952369
# cor4  nfl   DDX21 pentanucleotide 0.0002651704 0.2457659
# cor5  nfl    DDX6  hexanucleotide 0.0021726092 0.2075139
# cor6  nfl  DROSHA  hexanucleotide 0.0067362900 0.1838579
# cor7  nfl   EIF3D   trinucleotide 0.0034170887 0.1983584
# cor8  nfl FAM120A  hexanucleotide 0.0098946893 0.1751692
# cor11 nfl  GEMIN5   trinucleotide 0.0480718635 0.1346677
# cor12 nfl   GPKOW   trinucleotide 0.0026627843 0.2034484
# cor13 nfl  GTF2F1  hexanucleotide 0.0122961514 0.1700875
# cor14 nfl  METTL1  hexanucleotide 0.0385547554 0.1408859
# cor16 nfl  METTL1   trinucleotide 0.0025601169 0.2042401
# cor19 nfl   MTPAP   trinucleotide 0.0046615695 0.1918471
# cor20 nfl    NIP7  hexanucleotide 0.0028784478 0.2018719
# cor21 nfl    NIP7 pentanucleotide 0.0007577065 0.2274585
# cor23 nfl   NOL12  hexanucleotide 0.0012914750 0.2175816
# cor24 nfl   NOL12 pentanucleotide 0.0077966849 0.1805980
# cor26 nfl    NPM1   trinucleotide 0.0031331632 0.2001423
# cor27 nfl   PCBP1  hexanucleotide 0.0118729632 0.1709154
# cor29 nfl   PCBP2  hexanucleotide 0.0129360372 0.1688820
# cor31 nfl    PUM1   trinucleotide 0.0037246213 0.1965715
# cor33 nfl   RPS11  hexanucleotide 0.0036321758 0.1970942
# cor34 nfl   RPS11 pentanucleotide 0.0219981203 0.1557932
# cor38 nfl  TROVE2   trinucleotide 0.0300901883 0.1476221
# cor39 nfl   ZC3H8 pentanucleotide 0.0040570652 0.1947838


all_cox_validations <- lapply(scores, function(x){
    als_rsf_surv_df$score <- x
    fit <- cph(Surv(time,status) ~ nfl + score,
           data=als_rsf_surv_df,
           x=TRUE, y=TRUE, surv=TRUE)
    validate(fit, method="boot", B=300)
})
# $`AATF quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3285   0.3320 0.3247   0.0073          0.3212  0.2399 0.4033 300
# R2        0.2746   0.2821 0.2707   0.0114          0.2632  0.1717 0.3445 300
# Slope     1.0000   1.0000 0.9768   0.0232          0.9768  0.7402 1.2675 300
# D         0.0360   0.0374 0.0354   0.0020          0.0340  0.0189 0.0466 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0036 300
# Q         0.0371   0.0385 0.0348   0.0036          0.0334  0.0163 0.0460 300
# g         0.7445   0.7616 0.7374   0.0242          0.7204  0.5367 0.8807 300

# $`CPEB4 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3501   0.3500 0.3489   0.0011          0.3490  0.2716 0.4311 300
# R2        0.3142   0.3151 0.3106   0.0045          0.3097  0.2170 0.3950 300
# Slope     1.0000   1.0000 0.9841   0.0159          0.9841  0.7468 1.2617 300
# D         0.0424   0.0428 0.0418   0.0010          0.0414  0.0254 0.0552 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0006 0.0041 300
# Q         0.0435   0.0439 0.0412   0.0026          0.0408  0.0226 0.0543 300
# g         0.8099   0.8161 0.8030   0.0131          0.7968  0.6117 0.9679 300

# $`CPSF6 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3804   0.3819 0.3789   0.0030          0.3775  0.2969 0.4599 300
# R2        0.3558   0.3609 0.3532   0.0077          0.3481  0.2498 0.4444 300
# Slope     1.0000   1.0000 0.9823   0.0177          0.9823  0.7544 1.2505 300
# D         0.0495   0.0508 0.0491   0.0017          0.0478  0.0293 0.0643 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0005 0.0049 300
# Q         0.0506   0.0518 0.0484   0.0034          0.0472  0.0258 0.0634 300
# g         0.8921   0.9102 0.8866   0.0237          0.8684  0.6482 1.0739 300

# $`CSTF2 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3201   0.3267 0.3188   0.0079          0.3122  0.2333 0.3935 300
# R2        0.2640   0.2727 0.2610   0.0117          0.2523  0.1676 0.3355 300
# Slope     1.0000   1.0000 0.9724   0.0276          0.9724  0.7483 1.2408 300
# D         0.0344   0.0359 0.0339   0.0020          0.0324  0.0186 0.0450 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0004 0.0032 300
# Q         0.0354   0.0370 0.0334   0.0036          0.0318  0.0161 0.0444 300
# g         0.7256   0.7447 0.7210   0.0236          0.7019  0.5286 0.8672 300

# $`DDX21 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3063   0.3040 0.3023   0.0017          0.3046  0.2232 0.3852 300
# R2        0.2476   0.2495 0.2450   0.0045          0.2432  0.1525 0.3271 300
# Slope     1.0000   1.0000 0.9944   0.0056          0.9944  0.7273 1.2859 300
# D         0.0319   0.0323 0.0315   0.0009          0.0310  0.0165 0.0433 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0053 300
# Q         0.0329   0.0334 0.0309   0.0025          0.0304  0.0133 0.0425 300
# g         0.6918   0.6992 0.6891   0.0101          0.6817  0.4857 0.8462 300

# $`DDX6 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3649   0.3673 0.3637   0.0036          0.3613  0.2836 0.4479 300
# R2        0.3381   0.3421 0.3355   0.0066          0.3316  0.2457 0.4321 300
# Slope     1.0000   1.0000 0.9834   0.0166          0.9834  0.7530 1.2960 300
# D         0.0464   0.0474 0.0460   0.0014          0.0450  0.0295 0.0615 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0005 0.0037 300
# Q         0.0475   0.0485 0.0454   0.0031          0.0444  0.0274 0.0603 300
# g         0.8585   0.8730 0.8527   0.0203          0.8382  0.6589 1.0406 300

# $`DROSHA hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3354   0.3364 0.3342   0.0022          0.3332  0.2593 0.4159 300
# R2        0.2956   0.2990 0.2926   0.0064          0.2892  0.2061 0.3786 300
# Slope     1.0000   1.0000 0.9867   0.0133          0.9867  0.7619 1.2888 300
# D         0.0394   0.0401 0.0389   0.0013          0.0381  0.0241 0.0520 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0039 300
# Q         0.0404   0.0412 0.0383   0.0028          0.0376  0.0216 0.0508 300
# g         0.7782   0.7891 0.7722   0.0168          0.7613  0.5863 0.9370 300

# $`EIF3D trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3418   0.3397 0.3424  -0.0027          0.3445  0.2612 0.4331 300
# R2        0.3125   0.3095 0.3092   0.0003          0.3122  0.2151 0.4088 300
# Slope     1.0000   1.0000 1.0031  -0.0031          1.0031  0.7748 1.2978 300
# D         0.0421   0.0419 0.0416   0.0003          0.0418  0.0252 0.0572 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0004 0.0038 300
# Q         0.0432   0.0430 0.0411   0.0019          0.0413  0.0231 0.0560 300
# g         0.8207   0.8193 0.8138   0.0055          0.8152  0.6124 1.0021 300

# $`FAM120A hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3421   0.3420 0.3390   0.0031          0.3390  0.2572 0.4244 300
# R2        0.3076   0.3111 0.3043   0.0069          0.3008  0.2059 0.3907 300
# Slope     1.0000   1.0000 0.9868   0.0132          0.9868  0.7471 1.2982 300
# D         0.0413   0.0422 0.0408   0.0014          0.0399  0.0234 0.0543 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0005 0.0037 300
# Q         0.0424   0.0432 0.0402   0.0031          0.0393  0.0206 0.0531 300
# g         0.7997   0.8124 0.7931   0.0193          0.7804  0.5836 0.9570 300

# $`GARS pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3057   0.3050 0.3032   0.0017          0.3039  0.2286 0.3896 300
# R2        0.2509   0.2546 0.2480   0.0066          0.2443  0.1573 0.3286 300
# Slope     1.0000   1.0000 0.9845   0.0155          0.9845  0.7211 1.2856 300
# D         0.0324   0.0331 0.0319   0.0012          0.0312  0.0173 0.0437 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0042 300
# Q         0.0334   0.0342 0.0313   0.0028          0.0306  0.0148 0.0430 300
# g         0.7023   0.7118 0.6965   0.0152          0.6870  0.4973 0.8477 300

# $`GEMIN5 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3203   0.3249 0.3183   0.0066          0.3138  0.2207 0.4089 300
# R2        0.2648   0.2741 0.2612   0.0129          0.2520  0.1460 0.3488 300
# Slope     1.0000   1.0000 0.9753   0.0247          0.9753  0.7259 1.3069 300
# D         0.0345   0.0362 0.0339   0.0023          0.0322  0.0147 0.0469 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0006 0.0042 300
# Q         0.0355   0.0373 0.0333   0.0040          0.0316  0.0120 0.0463 300
# g         0.7264   0.7489 0.7213   0.0276          0.6989  0.4992 0.8821 300

# $`GEMIN5 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3526   0.3579 0.3531   0.0048          0.3478  0.2709 0.4208 300
# R2        0.3151   0.3195 0.3119   0.0076          0.3076  0.2128 0.3867 300
# Slope     1.0000   1.0000 0.9846   0.0154          0.9846  0.7530 1.2444 300
# D         0.0426   0.0435 0.0420   0.0015          0.0411  0.0246 0.0541 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0035 300
# Q         0.0436   0.0446 0.0415   0.0031          0.0405  0.0219 0.0533 300
# g         0.8090   0.8225 0.8027   0.0198          0.7893  0.5945 0.9485 300

# $`GPKOW trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3654   0.3671 0.3647   0.0024          0.3630  0.2775 0.4466 300
# R2        0.3452   0.3467 0.3414   0.0053          0.3399  0.2412 0.4392 300
# Slope     1.0000   1.0000 0.9866   0.0134          0.9866  0.7446 1.2554 300
# D         0.0477   0.0483 0.0470   0.0013          0.0464  0.0284 0.0631 300
# U        -0.0011  -0.0011 0.0007  -0.0017          0.0007 -0.0006 0.0054 300
# Q         0.0487   0.0493 0.0464   0.0030          0.0457  0.0248 0.0621 300
# g         0.8718   0.8828 0.8639   0.0189          0.8529  0.6395 1.0532 300

# $`GTF2F1 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3413   0.3445 0.3378   0.0068          0.3345  0.2562 0.4273 300
# R2        0.3067   0.3141 0.3034   0.0108          0.2960  0.2031 0.3906 300
# Slope     1.0000   1.0000 0.9805   0.0195          0.9805  0.7522 1.2729 300
# D         0.0412   0.0427 0.0406   0.0020          0.0391  0.0232 0.0543 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0042 300
# Q         0.0422   0.0437 0.0400   0.0037          0.0386  0.0207 0.0534 300
# g         0.7994   0.8184 0.7932   0.0252          0.7742  0.5746 0.9521 300

# $`METTL1 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3370   0.3401 0.3344   0.0057          0.3313  0.2521 0.4104 300
# R2        0.3027   0.3090 0.2996   0.0094          0.2932  0.2049 0.3829 300
# Slope     1.0000   1.0000 0.9849   0.0151          0.9849  0.7513 1.2541 300
# D         0.0405   0.0418 0.0400   0.0018          0.0387  0.0236 0.0530 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0036 300
# Q         0.0416   0.0429 0.0395   0.0034          0.0382  0.0213 0.0523 300
# g         0.8081   0.8237 0.8014   0.0223          0.7858  0.5989 0.9671 300

# $`METTL1 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3167   0.3197 0.3127   0.0070          0.3097  0.2253 0.3983 300
# R2        0.2593   0.2702 0.2563   0.0140          0.2453  0.1467 0.3330 300
# Slope     1.0000   1.0000 0.9721   0.0279          0.9721  0.7227 1.2511 300
# D         0.0336   0.0356 0.0332   0.0024          0.0312  0.0152 0.0446 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0006 0.0042 300
# Q         0.0347   0.0366 0.0326   0.0040          0.0307  0.0120 0.0440 300
# g         0.7135   0.7353 0.7088   0.0264          0.6871  0.4896 0.8535 300

# $`METTL1 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3559   0.3614 0.3546   0.0068          0.3492  0.2807 0.4339 300
# R2        0.3330   0.3411 0.3299   0.0112          0.3218  0.2339 0.4146 300
# Slope     1.0000   1.0000 0.9767   0.0233          0.9767  0.7485 1.2241 300
# D         0.0456   0.0472 0.0450   0.0022          0.0434  0.0274 0.0587 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0005 0.0052 300
# Q         0.0466   0.0483 0.0444   0.0039          0.0428  0.0238 0.0578 300
# g         0.8509   0.8731 0.8445   0.0286          0.8224  0.6248 1.0056 300

# $`MORC2 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3257   0.3278 0.3210   0.0068          0.3189  0.2344 0.4115 300
# R2        0.2645   0.2710 0.2613   0.0097          0.2548  0.1606 0.3446 300
# Slope     1.0000   1.0000 0.9826   0.0174          0.9826  0.7298 1.2722 300
# D         0.0344   0.0357 0.0339   0.0017          0.0327  0.0175 0.0464 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0039 300
# Q         0.0355   0.0367 0.0334   0.0034          0.0321  0.0151 0.0458 300
# g         0.7223   0.7371 0.7178   0.0193          0.7031  0.5169 0.8763 300

# $`MTPAP quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3260   0.3298 0.3221   0.0077          0.3183  0.2401 0.3967 300
# R2        0.2706   0.2795 0.2665   0.0130          0.2576  0.1728 0.3376 300
# Slope     1.0000   1.0000 0.9716   0.0284          0.9716  0.7387 1.2001 300
# D         0.0354   0.0370 0.0348   0.0022          0.0332  0.0193 0.0455 300
# U        -0.0011  -0.0011 0.0005  -0.0015          0.0005 -0.0004 0.0037 300
# Q         0.0364   0.0380 0.0343   0.0038          0.0327  0.0169 0.0451 300
# g         0.7364   0.7562 0.7295   0.0267          0.7096  0.5353 0.8702 300

# $`MTPAP trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3569   0.3587 0.3545   0.0042          0.3527  0.2746 0.4303 300
# R2        0.3271   0.3313 0.3238   0.0076          0.3195  0.2349 0.4075 300
# Slope     1.0000   1.0000 0.9816   0.0184          0.9816  0.7853 1.2473 300
# D         0.0446   0.0455 0.0440   0.0015          0.0430  0.0282 0.0576 300
# U        -0.0011  -0.0011 0.0005  -0.0015          0.0005 -0.0004 0.0029 300
# Q         0.0456   0.0466 0.0435   0.0031          0.0426  0.0264 0.0569 300
# g         0.8427   0.8566 0.8360   0.0206          0.8222  0.6494 0.9959 300

# $`NIP7 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3482   0.3519 0.3468   0.0052          0.3430  0.2622 0.4297 300
# R2        0.3142   0.3205 0.3111   0.0094          0.3048  0.2126 0.3977 300
# Slope     1.0000   1.0000 0.9811   0.0189          0.9811  0.7584 1.2357 300
# D         0.0424   0.0437 0.0419   0.0018          0.0406  0.0244 0.0557 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0036 300
# Q         0.0435   0.0448 0.0414   0.0034          0.0401  0.0217 0.0548 300
# g         0.8203   0.8356 0.8152   0.0204          0.7999  0.6023 0.9806 300

# $`NIP7 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3099   0.3136 0.3086   0.0050          0.3049  0.2215 0.3925 300
# R2        0.2531   0.2596 0.2509   0.0087          0.2444  0.1523 0.3368 300
# Slope     1.0000   1.0000 0.9752   0.0248          0.9752  0.7292 1.3211 300
# D         0.0327   0.0339 0.0324   0.0016          0.0311  0.0163 0.0449 300
# U        -0.0011  -0.0011 0.0007  -0.0017          0.0007 -0.0006 0.0051 300
# Q         0.0337   0.0350 0.0317   0.0033          0.0305  0.0128 0.0439 300
# g         0.7085   0.7260 0.7040   0.0220          0.6865  0.4932 0.8681 300

# $`NIPBL quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3348   0.3344 0.3300   0.0044          0.3304  0.2511 0.4135 300
# R2        0.2797   0.2841 0.2758   0.0084          0.2713  0.1787 0.3559 300
# Slope     1.0000   1.0000 0.9834   0.0166          0.9834  0.7180 1.2587 300
# D         0.0368   0.0377 0.0362   0.0015          0.0353  0.0199 0.0484 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0006 0.0051 300
# Q         0.0379   0.0388 0.0356   0.0032          0.0346  0.0162 0.0476 300
# g         0.7504   0.7629 0.7435   0.0195          0.7309  0.5230 0.8994 300

# $`NOL12 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3486   0.3555 0.3468   0.0087          0.3399  0.2561 0.4291 300
# R2        0.3278   0.3363 0.3246   0.0117          0.3161  0.2236 0.4172 300
# Slope     1.0000   1.0000 0.9750   0.0250          0.9750  0.7491 1.2692 300
# D         0.0447   0.0464 0.0441   0.0023          0.0424  0.0259 0.0591 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0005 0.0042 300
# Q         0.0457   0.0475 0.0435   0.0040          0.0417  0.0231 0.0583 300
# g         0.8510   0.8734 0.8447   0.0288          0.8222  0.6267 1.0286 300

# $`NOL12 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3194   0.3156 0.3175  -0.0019          0.3213  0.2353 0.4077 300
# R2        0.2624   0.2635 0.2597   0.0039          0.2585  0.1610 0.3445 300
# Slope     1.0000   1.0000 0.9929   0.0071          0.9929  0.7319 1.3109 300
# D         0.0341   0.0345 0.0337   0.0008          0.0333  0.0177 0.0463 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0039 300
# Q         0.0352   0.0356 0.0331   0.0024          0.0327  0.0152 0.0454 300
# g         0.7227   0.7295 0.7182   0.0112          0.7115  0.5162 0.8842 300

# $`NPM1 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3386   0.3401 0.3354   0.0047          0.3339  0.2588 0.4122 300
# R2        0.2900   0.2953 0.2856   0.0097          0.2803  0.1940 0.3593 300
# Slope     1.0000   1.0000 0.9757   0.0243          0.9757  0.7475 1.2122 300
# D         0.0385   0.0395 0.0378   0.0018          0.0367  0.0222 0.0492 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0043 300
# Q         0.0395   0.0406 0.0372   0.0033          0.0362  0.0196 0.0488 300
# g         0.7586   0.7730 0.7511   0.0219          0.7367  0.5588 0.8902 300

# $`NPM1 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3591   0.3619 0.3580   0.0039          0.3552  0.2710 0.4346 300
# R2        0.3337   0.3364 0.3302   0.0061          0.3276  0.2286 0.4192 300
# Slope     1.0000   1.0000 0.9873   0.0127          0.9873  0.7382 1.2755 300
# D         0.0457   0.0464 0.0451   0.0013          0.0443  0.0265 0.0596 300
# U        -0.0011  -0.0011 0.0007  -0.0017          0.0007 -0.0006 0.0051 300
# Q         0.0467   0.0475 0.0444   0.0031          0.0437  0.0231 0.0587 300
# g         0.8509   0.8604 0.8435   0.0169          0.8341  0.6284 1.0195 300

# $`PCBP1 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3306   0.3245 0.3271  -0.0026          0.3332  0.2516 0.4129 300
# R2        0.2961   0.2931 0.2926   0.0005          0.2956  0.2058 0.3872 300
# Slope     1.0000   1.0000 0.9995   0.0005          0.9995  0.7422 1.3029 300
# D         0.0394   0.0392 0.0389   0.0003          0.0391  0.0241 0.0536 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0037 300
# Q         0.0405   0.0403 0.0383   0.0019          0.0385  0.0219 0.0525 300
# g         0.7809   0.7790 0.7741   0.0049          0.7761  0.5906 0.9613 300

# $`PCBP1 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3336   0.3321 0.3297   0.0024          0.3312  0.2515 0.4134 300
# R2        0.2799   0.2821 0.2759   0.0063          0.2736  0.1847 0.3633 300
# Slope     1.0000   1.0000 0.9899   0.0101          0.9899  0.7331 1.2668 300
# D         0.0368   0.0374 0.0362   0.0012          0.0356  0.0210 0.0495 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0005 0.0064 300
# Q         0.0379   0.0385 0.0356   0.0029          0.0350  0.0174 0.0487 300
# g         0.7472   0.7569 0.7407   0.0162          0.7311  0.5441 0.9039 300

# $`PCBP2 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3493   0.3516 0.3468   0.0048          0.3445  0.2662 0.4257 300
# R2        0.3187   0.3224 0.3153   0.0071          0.3115  0.2211 0.3983 300
# Slope     1.0000   1.0000 0.9802   0.0198          0.9802  0.7616 1.2308 300
# D         0.0431   0.0440 0.0426   0.0014          0.0417  0.0259 0.0558 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0036 300
# Q         0.0442   0.0451 0.0421   0.0030          0.0412  0.0236 0.0552 300
# g         0.8256   0.8372 0.8188   0.0184          0.8072  0.6263 0.9812 300

# $`PCBP2 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3378   0.3385 0.3345   0.0040          0.3338  0.2641 0.4130 300
# R2        0.2899   0.2937 0.2858   0.0079          0.2820  0.2027 0.3670 300
# Slope     1.0000   1.0000 0.9812   0.0188          0.9812  0.7489 1.2574 300
# D         0.0384   0.0392 0.0378   0.0015          0.0370  0.0238 0.0502 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0034 300
# Q         0.0395   0.0403 0.0373   0.0030          0.0364  0.0218 0.0495 300
# g         0.7681   0.7793 0.7607   0.0186          0.7495  0.5799 0.9155 300

# $`PUM1 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3594   0.3569 0.3571  -0.0001          0.3595  0.2758 0.4422 300
# R2        0.3181   0.3177 0.3143   0.0033          0.3148  0.2165 0.4048 300
# Slope     1.0000   1.0000 0.9875   0.0125          0.9875  0.7412 1.2394 300
# D         0.0431   0.0432 0.0424   0.0008          0.0422  0.0251 0.0567 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0044 300
# Q         0.0441   0.0443 0.0419   0.0024          0.0417  0.0220 0.0560 300
# g         0.8133   0.8178 0.8067   0.0111          0.8022  0.6015 0.9792 300

# $`RBM22 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3351   0.3396 0.3332   0.0064          0.3286  0.2452 0.4160 300
# R2        0.2908   0.2966 0.2868   0.0098          0.2810  0.1859 0.3717 300
# Slope     1.0000   1.0000 0.9778   0.0222          0.9778  0.7289 1.2714 300
# D         0.0386   0.0398 0.0380   0.0018          0.0368  0.0207 0.0510 300
# U        -0.0011  -0.0011 0.0006  -0.0017          0.0006 -0.0006 0.0046 300
# Q         0.0396   0.0408 0.0373   0.0035          0.0361  0.0178 0.0504 300
# g         0.7735   0.7891 0.7657   0.0234          0.7501  0.5496 0.9293 300

# $`RPS11 hexanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3392   0.3367 0.3369  -0.0002          0.3394  0.2603 0.4265 300
# R2        0.2983   0.2984 0.2951   0.0033          0.2950  0.2050 0.3822 300
# Slope     1.0000   1.0000 0.9937   0.0063          0.9937  0.7413 1.2917 300
# D         0.0398   0.0400 0.0393   0.0008          0.0390  0.0238 0.0528 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0041 300
# Q         0.0408   0.0411 0.0387   0.0024          0.0384  0.0211 0.0515 300
# g         0.7884   0.7923 0.7827   0.0097          0.7787  0.5947 0.9538 300

# $`RPS11 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3060   0.3095 0.3058   0.0037          0.3023  0.2197 0.3927 300
# R2        0.2533   0.2582 0.2506   0.0076          0.2457  0.1595 0.3370 300
# Slope     1.0000   1.0000 0.9808   0.0192          0.9808  0.7340 1.2833 300
# D         0.0327   0.0337 0.0323   0.0014          0.0313  0.0171 0.0446 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0052 300
# Q         0.0338   0.0347 0.0317   0.0030          0.0308  0.0139 0.0438 300
# g         0.7097   0.7216 0.7051   0.0164          0.6933  0.5157 0.8773 300

# $`RPS6 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3300   0.3275 0.3251   0.0024          0.3276  0.2459 0.4044 300
# R2        0.2738   0.2756 0.2693   0.0063          0.2675  0.1771 0.3514 300
# Slope     1.0000   1.0000 0.9847   0.0153          0.9847  0.7314 1.2763 300
# D         0.0359   0.0364 0.0352   0.0012          0.0347  0.0198 0.0474 300
# U        -0.0011  -0.0011 0.0006  -0.0016          0.0006 -0.0005 0.0046 300
# Q         0.0369   0.0374 0.0346   0.0028          0.0341  0.0168 0.0465 300
# g         0.7429   0.7487 0.7349   0.0139          0.7290  0.5368 0.8996 300

# $`SF3B1 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3329   0.3323 0.3304   0.0019          0.3310  0.2546 0.4150 300
# R2        0.2785   0.2788 0.2750   0.0039          0.2746  0.1856 0.3591 300
# Slope     1.0000   1.0000 0.9883   0.0117          0.9883  0.7407 1.2755 300
# D         0.0366   0.0369 0.0361   0.0008          0.0358  0.0213 0.0489 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0036 300
# Q         0.0377   0.0379 0.0355   0.0024          0.0353  0.0188 0.0479 300
# g         0.7443   0.7477 0.7385   0.0091          0.7352  0.5436 0.8957 300

# $`STAU2 quadnucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3257   0.3299 0.3230   0.0070          0.3187  0.2377 0.3985 300
# R2        0.2788   0.2872 0.2760   0.0112          0.2676  0.1824 0.3479 300
# Slope     1.0000   1.0000 0.9779   0.0221          0.9779  0.7491 1.2557 300
# D         0.0367   0.0382 0.0362   0.0020          0.0347  0.0206 0.0472 300
# U        -0.0011  -0.0011 0.0005  -0.0016          0.0005 -0.0005 0.0034 300
# Q         0.0377   0.0393 0.0357   0.0036          0.0342  0.0184 0.0467 300
# g         0.7549   0.7742 0.7499   0.0243          0.7306  0.5573 0.8895 300

# $`TROVE2 trinucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3349   0.3371 0.3319   0.0052          0.3297  0.2510 0.4159 300
# R2        0.2920   0.2969 0.2883   0.0086          0.2835  0.1930 0.3746 300
# Slope     1.0000   1.0000 0.9824   0.0176          0.9824  0.7695 1.2454 300
# D         0.0388   0.0398 0.0382   0.0016          0.0372  0.0219 0.0514 300
# U        -0.0011  -0.0011 0.0005  -0.0015          0.0005 -0.0004 0.0032 300
# Q         0.0398   0.0409 0.0377   0.0031          0.0367  0.0199 0.0508 300
# g         0.7776   0.7877 0.7698   0.0179          0.7596  0.5731 0.9342 300

# $`ZC3H8 pentanucleotide`
#       index.orig training   test optimism index.corrected   Lower  Upper   n
# Dxy       0.3215   0.3186 0.3187   0.0000          0.3216  0.2402 0.4052 300
# R2        0.2606   0.2623 0.2572   0.0050          0.2556  0.1628 0.3459 300
# Slope     1.0000   1.0000 0.9821   0.0179          0.9821  0.7044 1.3008 300
# D         0.0338   0.0343 0.0333   0.0010          0.0328  0.0176 0.0461 300
# U        -0.0011  -0.0011 0.0007  -0.0017          0.0007 -0.0006 0.0045 300
# Q         0.0349   0.0354 0.0326   0.0027          0.0322  0.0148 0.0456 300
# g         0.7201   0.7272 0.7150   0.0122          0.7080  0.5235 0.8883 300

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

	fit_baseline <- coxph(Surv(time, status) ~ nfl,
	    data = als_rsf_surv_df, x = TRUE)
	fit_ridge <- coxph(Surv(time, status) ~ nfl + ridge_score,
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

saveRDS(object = return_list, file = file.path(save_dir, "07-02-save-all-ridge-coefficient-models-from-answerals.rds"))
