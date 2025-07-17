data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "40-giulia-rbp-kd")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
samplesheet_df <- read.csv("/nemo/lab/patanir/home/users/ziffo/patani-collab/motor-neuron-rbp-knockdown/sample-details/samplesheet.csv")
samples_param_df <- read.csv(file.path(analysis_input_dir, "01-01-public-data-params-dataframe.txt"))
results_list <- readRDS(file.path(analysis_save_dir, "01-02-list_of_counts_stat_samplesheet_from_answerals.rds"))
sig_result_list <- readRDS(file.path("../20-import-significant-events/rdata_files/02-significantly-differentially-expressed-events.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
require("DESeq2")
source(file.path(analysis_dir, "scripts", "src", "experiment_matrix_functions.R"))

sample_df <- results_list$sample_df
counts_matrix <- results_list$counts
min_count <- 1

sig_result_df <- sig_result_list$significant_events
repeat_containing_introns_df <- sig_result_df[sig_result_df$status == "repeat_containing" & 
    sig_result_df$region == "intron",]
repeat_containing_introns_df$hash_id <- feature_ranges$hash_id[match(repeat_containing_introns_df$feature_ids, feature_ranges$feature_id)]



experiment_df <- samplesheet_df[match(sample_df$sample_name, samplesheet_df$sample),]
experiment_df <- experiment_df[experiment_df$fraction == "whole",]
experiment_df <- experiment_df[,c("sample", "fraction", "cellline", "condition")]
experiment_df$condition <- factor(experiment_df$condition, levels = c("untreated", "scramble", "sfpqkd", "fuskd", "tdp43kd"))


wholecell_count_matrix <- counts_matrix[,sample_df$column_id[match(experiment_df$sample, sample_df$sample_name)]]
colnames(wholecell_count_matrix) <- experiment_df$sample



wholecell_count_matrix <- wholecell_count_matrix[rowMeans(wholecell_count_matrix) > min_count,]
deseq_object <- DESeqDataSetFromMatrix(countData = round(wholecell_count_matrix), 
    colData = experiment_df, 
    design = ~ 0 + cellline + condition)
differential_object <- DESeq(deseq_object)


# col_data <- colData(differential_object)
model_matrix <- model.matrix(design(differential_object), colData(differential_object))

sfpqkd_vs_scramble <- colMeans(model_matrix[differential_object$condition == "sfpqkd",]) - colMeans(model_matrix[differential_object$condition == "scramble",])
fuskd_vs_scramble <- colMeans(model_matrix[differential_object$condition == "fuskd",]) - colMeans(model_matrix[differential_object$condition == "scramble",])
tdp43kd_vs_scramble <- colMeans(model_matrix[differential_object$condition == "tdp43kd",]) - colMeans(model_matrix[differential_object$condition == "scramble",])

contrast_list <- list(
    "sfpqkd_vs_scramble" = sfpqkd_vs_scramble,
    "fuskd_vs_scramble" = fuskd_vs_scramble,
    "tdp43kd_vs_scramble" = tdp43kd_vs_scramble)

all_contrasts_list <- lapply(contrast_list, function(a_contrast){
    current_results <- results(object = differential_object, contrast = a_contrast)
    return(current_results)
})

lapply(all_contrasts_list, function(x){
    nrow(x[x$padj < 0.05 & !is.na(x$padj),])
})
# $sfpqkd_vs_scramble
# [1] 44788

# $fuskd_vs_scramble
# [1] 898

# $tdp43kd_vs_scramble
# [1] 1576

lapply(all_contrasts_list, function(x){
    nrow(x[x$padj < 0.05 & !is.na(x$padj) & rownames(x) %in% repeat_containing_introns_df$hash_id,])
})
# $sfpqkd_vs_scramble
# [1] 1272

# $fuskd_vs_scramble
# [1] 30

# $tdp43kd_vs_scramble
# [1] 28


lapply(all_contrasts_list, function(x){
    nrow(x[x$padj < 0.05 & !is.na(x$padj) & rownames(x) %in% repeat_containing_introns_df$hash_id,]) / nrow(x[x$padj < 0.05 & !is.na(x$padj),])
})


lapply(all_contrasts_list, function(x){
    table(x[x$padj < 0.05 & !is.na(x$padj) & rownames(x) %in% repeat_containing_introns_df$hash_id,"log2FoldChange"] > 0)
})


saveRDS(object = list(contrasts = all_contrasts_list, experiment_matrix = colData(differential_object), count_matrix = wholecell_count_matrix), file = file.path(analysis_save_dir, "01-03-list-of-contrasts-for-sfpqkd-fuskd-and-tdp43kd.rds"))

