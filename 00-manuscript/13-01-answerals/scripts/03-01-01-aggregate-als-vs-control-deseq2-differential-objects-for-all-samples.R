# ==========================================================================
# Intial script setup
# ==========================================================================
data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "13-01-answerals")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
isPairedEnd <- T
strand_specificity <- 2
sampleinfo_df <- read.csv("/nemo/project/proj-luscombn-patani/working/answerals/sample-details/samplesheet.csv")
alignment_dir <- "/nemo/project/proj-luscombn-patani/working/"
samples <- c("answerals" = "answerals")
object_list <- readRDS(file.path(analysis_save_dir, "01-02-list_of_counts_stat_samplesheet_from_answerals.rds"))
# ==========================================================================
# Imports
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
require("DESeq2")
source(file.path(analysis_dir, "scripts", "src", "experiment_matrix_functions.R"))
source(file.path(analysis_dir, "../05-simulated-reads-differential-profiling/scripts", "src", "deseq_simulated_reads_differential_expression_profiling.R"))
# ==========================================================================
# Analysis
# ==========================================================================


# Simple model differential expression 
simple_model_differential_object <- readRDS(file.path(analysis_save_dir, "02-01-01-save-deseq2-simple-differential-object.rds"))


contrast_list <- list()

model_matrix <- model.matrix(design(simple_model_differential_object), colData(simple_model_differential_object))
contrast_list[["sals_vs_control"]] <- colMeans(model_matrix[simple_model_differential_object$condition_group == "sals",]) - colMeans(model_matrix[simple_model_differential_object$condition_group == "ctrl",])
contrast_list[["fals_vs_control"]] <- colMeans(model_matrix[simple_model_differential_object$condition_group == "fals",]) - colMeans(model_matrix[simple_model_differential_object$condition_group == "ctrl",])
contrast_list[["mnd_vs_control"]] <- colMeans(model_matrix[simple_model_differential_object$condition_group == "other_mnd",]) - colMeans(model_matrix[simple_model_differential_object$condition_group == "ctrl",])
contrast_list[["sals_vs_fals"]] <- colMeans(model_matrix[simple_model_differential_object$condition_group == "sals",]) - colMeans(model_matrix[simple_model_differential_object$condition_group == "fals",])


all_results_df <- do.call(rbind, lapply(names(contrast_list), function(y){
    current_results <- results(object = simple_model_differential_object, contrast = contrast_list[[y]])
    filter <- !is.na(current_results$pvalue)
    return_df <- data.frame(study = "answerals",
        contrasts = y, 
        feature_ids = rownames(current_results)[filter],
        logfc = current_results$log2FoldChange[filter],
        pval = current_results$pvalue[filter],
        padj = current_results$padj[filter])
    return(return_df)
}))

saveRDS(object = all_results_df, file = file.path(analysis_save_dir, "03-01-01-answerals-get-als-vs-control-pairwise-comparisons.rds"))
