# ==========================================================================
# Intial script setup
# ==========================================================================
data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "13-nygc-consortium")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
samples_param_df <- read.csv(file.path(analysis_input_dir, "01-01-nygc-params-dataframe.txt"))
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
source(file.path(project_dir, "analysis", "src", "subread_featurecounts_parsing_functions.R"))
slurm_job_list <- readRDS(file.path(analysis_save_dir, "03-submit-jobs-for-exonic-intronic-fractional-counts.rds"))
slurm_dir <- file.path(rslurm_dir, "featurecounts_of_exons_and_introns")

result_list <- assemble_all_counts_from_jobs(slurm_job_list, slurm_dir, samples_param_df)

saveRDS(object = result_list, 
	file = file.path(analysis_save_dir, 
	"03-02-list_of_counts_stat_samplesheet_from_nygc.rds"))