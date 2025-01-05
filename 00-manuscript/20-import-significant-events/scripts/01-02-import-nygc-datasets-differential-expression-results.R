data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "20-import-significant-events")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
repeat_object_list <- readRDS(file.path("../01-repeat-characterisation/rdata_files", "00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds"))
nygc_results_filepath_df <- readRDS(file.path("../13-nygc-consortium/rdata_files/03-04-differentially-expressed-events-filepaths.rds"))
# object_list <- readRDS(file.path("../13-nygc-consortium/rdata_files/03-02-list_of_counts_stat_samplesheet_from_nygc.rds"))
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
source("../../useful-functions/plot-utils.R")
source("../../useful-functions/basepair_coverage_functions.R")
require("cowplot")
require("ggplot2")
require("GenomicRanges")
require("HiCBricks")
require("stringr")
require("RColorBrewer")
require("rslurm")
require("dplyr")
require("rslurm")

repeat_ranges <- repeat_object_list$repeat_ranges

non_specific_results_filepath_df <- nygc_results_filepath_df[nygc_results_filepath_df$tissue == "non-specific",]
nygc_non_specific_results_df <- do.call(rbind, lapply(seq_len(nrow(non_specific_results_filepath_df)), function(x){
    a_row <- non_specific_results_filepath_df[x,]
    result_object <- readRDS(a_row$file_path)
    filter <- !is.na(result_object$pvalue)
    return_df <- data.frame(study = "nygc_tissue_non_specific",
        contrasts = paste(a_row$condition_1, "vs", a_row$condition_2, sep = "_"),
        feature_ids = rownames(result_object)[filter],
        logfc = result_object$log2FoldChange[filter],
        pval = result_object$pvalue[filter],
        padj = result_object$padj[filter])
}))

tissue_specific_results_filepath_df <- nygc_results_filepath_df[nygc_results_filepath_df$tissue != "non-specific",]
nygc_tissue_specific_results_df <- do.call(rbind, lapply(seq_len(nrow(tissue_specific_results_filepath_df)), function(x){
    a_row <- tissue_specific_results_filepath_df[x,]
    result_object <- readRDS(a_row$file_path)
    filter <- !is.na(result_object$pvalue)
    return_df <- data.frame(study = paste("nygc", a_row$tissue, sep = "_"),
        contrasts = paste(a_row$condition_1, "vs", a_row$condition_2, sep = "_"),
        feature_ids = rownames(result_object)[filter],
        logfc = result_object$log2FoldChange[filter],
        pval = result_object$pvalue[filter],
        padj = result_object$padj[filter])
}))

nygc_result_df <- rbind(nygc_non_specific_results_df, nygc_tissue_specific_results_df)

sig_nygc_results_df <- nygc_result_df[nygc_result_df$pval < 0.05 & !is.na(nygc_result_df$pval),]

replace_feature_ids <- feature_ranges$feature_id[match(sig_nygc_results_df$feature_ids, feature_ranges$hash_id)]
sig_nygc_results_df$feature_ids <- replace_feature_ids

saveRDS(object = sig_nygc_results_df, file = file.path(analysis_save_dir, "01-03-import-nygc-datasets-differential-expression.rds"))