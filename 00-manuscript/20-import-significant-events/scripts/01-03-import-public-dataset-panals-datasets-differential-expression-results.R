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
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
repeat_object_list <- readRDS(file.path("../01-repeat-characterisation/rdata_files", "00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds"))
panals_results_df <- readRDS(file.path("../13-03-public-data/rdata_files/03-01-02-assembled-contrasts-from-pan-als-comparison.rds"))
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
# feature_ranges <- object_list$ranges

replace_feature_ids <- feature_ranges$feature_id[match(panals_results_df$feature_ids, feature_ranges$hash_id)]
panals_results_df$feature_ids <- replace_feature_ids

saveRDS(object = panals_results_df, file = file.path(analysis_save_dir, "01-03-import-panals-public-datasets-differential-expression.rds"))

panals_results_df <- readRDS(file.path(analysis_save_dir, "01-03-import-panals-public-datasets-differential-expression.rds"))
# panals_results_df[grepl(x = panals_results_df$feature_ids, pattern = "ENSG00000186868") & panals_results_df$pval < 0.05,]


# "ENSG00000186868:17:45961605:45961757:+:exon"