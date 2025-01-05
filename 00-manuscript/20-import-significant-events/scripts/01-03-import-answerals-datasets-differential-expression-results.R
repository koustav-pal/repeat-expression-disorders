data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "20-final-frequency-profiles-of-differential-expression")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
repeat_object_list <- readRDS(file.path("../01-repeat-characterisation/rdata_files", "00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds"))
answerals_als_vs_control_results_df <- readRDS(file.path("../13-01-answerals/rdata_files/03-01-01-answerals-get-als-vs-control-pairwise-comparisons.rds"))
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
feature_ranges <- object_list$ranges

replace_feature_ids <- feature_ranges$feature_id[match(answerals_als_vs_control_results_df$feature_ids, feature_ranges$hash_id)]
answerals_als_vs_control_results_df$feature_ids <- replace_feature_ids

saveRDS(object = answerals_als_vs_control_results_df, file = file.path(analysis_save_dir, "01-03-import-answerals-datasets-differential-expression.rds"))