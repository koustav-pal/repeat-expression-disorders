data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "13-03-public-data")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
samples_param_df <- read.csv(file.path(analysis_input_dir, "01-01-public-data-params-dataframe.txt"))
results_list <- readRDS(file.path(analysis_save_dir, "02-list_of_counts_stat_samplesheet_from_public_data.rds"))
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

all_functions_list <- return_function_list()

all_save_path_df <- do.call(rbind, lapply(seq_along(all_functions_list), function(x){
    sample_name <- names(all_functions_list)[[x]]
    message(sample_name)
    FUN <- all_functions_list[[sample_name]]
    experiment_df <- FUN(sample_df)$experiment
    save_path <- file.path(analysis_save_dir, paste("02-01", x, "save-deseq2", sample_name, "differential-object.rds", sep = "-"))
    file_path_df <- perform_diff_expr(counts_matrix, 
        sample_df, 
        FUN, 
        sample_name, 
        min_count = 1, 
        study_specific_design = T, 
        save_path = save_path)
    return(file_path_df)
}))

saveRDS(object = all_save_path_df, file = file.path(analysis_save_dir, "02-01-study-specific-differential-expression-profiling-paths.rds"))

