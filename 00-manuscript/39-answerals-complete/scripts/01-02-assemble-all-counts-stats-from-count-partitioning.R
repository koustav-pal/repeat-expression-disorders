library(here)
project_dir <- here()
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
analysis_dir <- file.path(project_dir, "39-answerals-complete")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS(here("05-simulated-reads-differential-profiling","rdata_files","01-all-merged-exonic-and-intronic-segments-in-genes.rds"))
samples_param_df <- read.csv(file.path(analysis_input_dir, "01-01-public-data-params-dataframe.txt"))
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
source(here("src", "subread_featurecounts_parsing_functions.R"))
slurm_job_list <- readRDS(file.path(analysis_save_dir, "01-01-submit-jobs-for-exonic-intronic-fractional-counts.rds"))
slurm_dir <- file.path(rslurm_dir, "featurecounts_of_exons_and_introns")

result_list <- assemble_all_counts_from_jobs(slurm_job_list, slurm_dir, samples_param_df)

saveRDS(object = result_list, 
	file = file.path(analysis_save_dir, 
	"01-02-list_of_counts_stat_samplesheet_from_answerals.rds"))