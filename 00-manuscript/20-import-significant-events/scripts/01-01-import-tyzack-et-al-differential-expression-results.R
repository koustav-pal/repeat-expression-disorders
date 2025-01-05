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
object_list <- readRDS(file.path("../05-simulated-reads-differential-profiling/rdata_files/02-03-list_of_featureCounts_matrix_in_exons_and_introns_by_molecule.rds"))
repeat_object_list <- readRDS(file.path("../01-repeat-characterisation/rdata_files", "00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds"))
fractional_counts_job_list <- readRDS("../09-01-QC-compare-fractional-counts-vs-abs-counts/rdata_files/01_submit_jobs_to_perform_fractional_counts_differential_expression_profiling_using_deseq2.rds")
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

slurm_dir <- file.path("../09-01-QC-compare-fractional-counts-vs-abs-counts/rslurm/differential_expression_profiling")
if(!dir.exists(slurm_dir)){
    dir.create(slurm_dir, recursive = T)
}
setwd(slurm_dir)
# file.remove(file.path("_rslurm_deseq_diff_expression", list.files("./_rslurm_deseq_diff_expression/")))
slurm_opts <- list(
    "time" = "72:00:00", 
    "partition" = "cpu",
    "mem-per-cpu" = "64GB")

fractional_counts_job_output <- get_slurm_out(fractional_counts_job_list)
# abs_counts_job_output <- get_slurm_out(abs_counts_job_list)

repeat_ranges <- repeat_object_list$repeat_ranges
feature_ranges <- object_list$ranges

all_de_events_df <- do.call(rbind, fractional_counts_job_output[c("vcp_d0", "vcp_d3", "vcp_d7", "vcp_d14", "vcp_d22", "vcp_d35")])
all_de_events_df$study <- "tyzack-et-al-2021"
all_de_events_df$contrasts <- all_de_events_df$contrast
all_de_events_df$logfc <- all_de_events_df$lfc
all_de_events_df$pval <- all_de_events_df$pvalue
return_df <- all_de_events_df[,c("study", "contrasts", "feature_ids", "logfc", "pval", "padj")] 

saveRDS(object = return_df, file = file.path(analysis_save_dir, "01-01-import-tyzack-et-al-differential-expression.rds"))