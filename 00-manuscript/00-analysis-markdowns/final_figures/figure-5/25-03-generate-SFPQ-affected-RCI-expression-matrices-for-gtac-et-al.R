library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- here("final-figures", "figure-5")
save_dir <- file.path(analysis_dir, "rdata_files")
input_dir <- file.path(analysis_dir,"input_files")
plot_dir <- file.path(analysis_dir, "plots")
output_dir <- file.path(analysis_dir, "output_files")
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

# ==========================================================================
# Load data
# ==========================================================================


#
# Load repeats and introns containing repeats
# --------------------------------------------------------------------------
repeat_ranges_list <- readRDS("final-figures/figure-1/rdata_files/00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds")
repeat_ranges <- unique(repeat_ranges_list$repeat_ranges)
feature_ranges_with_repeats <- repeat_ranges_list$repeat_containing_ranges

#
# Load total feature ranges
# --------------------------------------------------------------------------

feature_ranges <- readRDS("05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")

#
# Load significant events and repeat containing introns
# --------------------------------------------------------------------------

sig_result_list <- readRDS("20-import-significant-events/rdata_files/02-significantly-differentially-expressed-events.rds")
sig_result_df <- sig_result_list$significant_events

repeat_containing_introns_df <- sig_result_df[sig_result_df$status == "repeat_containing" & 
    sig_result_df$region == "intron",]
repeat_containing_introns_df$cohort <- "public-data"
repeat_containing_introns_df$cohort[grepl(x = repeat_containing_introns_df$study, pattern = "nygc")] <- "nygc"
repeat_containing_introns_df$cohort[grepl(x = repeat_containing_introns_df$study, pattern = "answerals")] <- "answerals"

repeat_containing_introns_df$tissue_type <- "iPSCs"
repeat_containing_introns_df$tissue_type[grepl(x = repeat_containing_introns_df$study, pattern = "nygc")] <- "postmortem"

#
# Load Grima et al data
# --------------------------------------------------------------------------

gtac_object_list <- readRDS("38-polyA-blood-pbmcs/rdata_files/01-list_of_counts_stat_samplesheet_from_polyA_pbmcs.rds")
patient_df <- gtac_object_list$sample_df

#
# Load SFPQ KD penetrant RCIs
# --------------------------------------------------------------------------

sfpq_kd_rcis_df <- readRDS(file.path("40-giulia-rbp-kd/rdata_files/04-06-introns-affected-by-SFPQKD-with-rloop-density-and-rci-penetrance.rds"))
sfpq_kd_rcis <- rownames(sfpq_kd_rcis_df[sfpq_kd_rcis_df$rci_status == "RCI",])
rloop_containing_sfpq_kd_rcis <- rownames(sfpq_kd_rcis_df)[sfpq_kd_rcis_df$rci_status == "RCI" & 
    sfpq_kd_rcis_df$rloop_status == "R-loop containing"]

# ==========================================================================
# Functions
# ==========================================================================

prepare_progression_rsf_df <- function(df){
    df$time <- df$disease_duration
    df$status <- 1
    df$status[df$Deceased != 1] <- 0
    return(df)
}

# ==========================================================================
# Generate TPM matrix
# ==========================================================================

counts_matrix <- gtac_object_list$counts
temp_matrix <- counts_matrix / (width(feature_ranges)[match(rownames(counts_matrix), feature_ranges$hash_id)] / 1e3)
tpm_matrix <- t(t(temp_matrix) / (colSums(temp_matrix) / 1e6))

keep_hash_ids <- feature_ranges_with_repeats$hash_id[match(unique(sfpq_kd_rcis), feature_ranges_with_repeats$hash_id)]
keep_tpm_matrix <- tpm_matrix[keep_hash_ids,colnames(tpm_matrix) %in% patient_df$column_id]

# ==========================================================================
# Construct DESeq2 object
# ==========================================================================
require(DESeq2)

keep_patient_df <- patient_df[patient_df$AFFECTION_STATUS == "1",]
keep_patient_df$sex <- factor(keep_patient_df$sex, levels = c("female", "male"))

deseq_object <- DESeqDataSetFromMatrix(countData = round(counts_matrix[keep_hash_ids,])[,keep_patient_df$column_id],
    colData = keep_patient_df,
    design = ~ sex)

diff_object <- DESeq(deseq_object)
normalized_counts <- counts(diff_object, normalized=TRUE)

# ==========================================================================
# Generate Penetrance matrix
# ==========================================================================

patient_penetrance_matrix <- do.call(rbind, lapply(seq_len(nrow(keep_patient_df)), function(x){
    temp_df <- keep_patient_df[x,]
    current_columns <- temp_df$column_id
    current_tpm_vector <- Reduce("+", lapply(current_columns, function(x){
        normalized_counts[,x]
    })) / length(current_columns)
    current_threshold <- 0.1
    a_vector <- rep(0, length(unique(sfpq_kd_rcis)))
    names(a_vector) <- unique(sfpq_kd_rcis)
    a_vector[names(current_tpm_vector[current_tpm_vector > current_threshold])] <- 1
    return(a_vector)
}))

# ==========================================================================
# Save data
# ==========================================================================

sfpq_affected_rcis_matrix <- normalized_counts[sfpq_kd_rcis,]
keep_patient_df$tissue_group <- "blood"
keep_patient_df$penetrance_row <- match(keep_patient_df$column_id, colnames(sfpq_affected_rcis_matrix))

saveRDS(object = list(patient_df = keep_patient_df, 
    penetrance_matrix = patient_penetrance_matrix,
    sfpq_affected_rcis_matrix = sfpq_affected_rcis_matrix,
    sfpq_affected_rcis = sfpq_kd_rcis,
    sfpq_affected_rloop_containing_rcis = rloop_containing_sfpq_kd_rcis),
    file = file.path(save_dir, "06-03-save-gtac-blood-pbmc-sfpq-affected-rci-expression-matrix.rds"))