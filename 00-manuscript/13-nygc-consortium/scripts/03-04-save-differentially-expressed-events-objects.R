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
analysis_plot_dir <- file.path(analysis_dir, "plots")
rslurm_dir <- file.path(analysis_dir, "rslurm")
alignment_dir <- file.path(project_dir, "analysis", "Alignments")
# object_list <- readRDS(file.path(analysis_save_dir, "02-03-list_of_featureCounts_matrix_in_exons_and_introns_by_patient.rds"))
isPairedEnd <- T
strand_specificity <- 2
# sampleinfo_df <- read.csv(file.path(analysis_input_dir, "00_filtered_sample_gsm_list.txt"))
# clinical_metadata_df <- read.csv("../0103-NYGC-consortium/input_files/samplesheet.csv")
repeat_object_list <- readRDS(file.path("../01-repeat-characterisation/rdata_files", "00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds"))
# hdf5_matrix <- HDF5Array(file.path(analysis_save_dir, "02-03-featurewise-counts-matrix.hdf5"), "feature_counts")
# ==========================================================================
# Imports
# ==========================================================================
require("GenomicRanges")
require("rslurm")
require("stringr")
require("HiCBricks")
require("PCAtools")
require("matrixStats")
require("HDF5Array")
require("DelayedArray")
require("DESeq2")
source("../09-01-QC-compare-fractional-counts-vs-abs-counts/scripts/src/deseq_functions.R")
source("../../useful-functions/plot-utils.R")

differential_object <- readRDS(file.path(analysis_save_dir, "03-03-save-deseq2-differential-object.rds"))

col_data_df <- colData(differential_object)

model_matrix <- model.matrix(design(differential_object), colData(differential_object))

control_results <- colMeans(model_matrix[differential_object$condition == "ctrl",])
# tissue specific expression changes
conditions <- unique(differential_object$condition)
conditions <- conditions[conditions != "ctrl"]
tissue_non_specific_result_df_list <- lapply(conditions, function(a_condition){
    condition_results <- colMeans(model_matrix[differential_object$condition == a_condition,])
    save_file <- file.path(analysis_save_dir, paste("03-04", a_condition, "vs", "control-differentially-expressed-events.rds", sep = "-"))
    if(!file.exists(save_file)){
        results_object <- results(differential_object, contrast = condition_results - control_results)
        saveRDS(object = results_object, file = save_file)  
    }
    return(data.frame(condition_1 = a_condition, 
        condition_2 = "control", 
        tissue = "non-specific", 
        file_path = save_file))
})
tissue_non_specific_result_df <- do.call(rbind, tissue_non_specific_result_df_list)

tissues <- unique(col_data_df$tissue)
conditions <- unique(differential_object$condition)
conditions <- conditions[conditions != "ctrl"]
tissue_specific_result_df_list <- lapply(unique(col_data_df$tissue), function(a_tissue){
    message(a_tissue)
    result_df_list <- do.call(rbind, lapply(conditions, function(a_condition){
            message(a_condition)
            save_file <- file.path(analysis_save_dir, paste("03-04", a_tissue, a_condition, "vs", "control-differentially-expressed-events.rds", sep = "-"))
            if(!file.exists(save_file)){
                message(save_file)
                condition_results <- colMeans(model_matrix[
                    differential_object$condition == a_condition & 
                    differential_object$tissue == a_tissue,])
                control_results <- colMeans(model_matrix[
                    differential_object$condition == "ctrl" & 
                    differential_object$tissue == a_tissue,])
                results_object <- results(differential_object, contrast = condition_results - control_results)
                saveRDS(object = results_object, file = save_file)
            }
            return(data.frame(condition_1 = a_condition, 
                condition_2 = "control", 
                tissue = a_tissue, 
                file_path = save_file))
        }))
    result_df_list
})
tissue_specific_result_df <- do.call(rbind, tissue_specific_result_df_list)
all_result_df <- rbind(tissue_non_specific_result_df, tissue_specific_result_df)

saveRDS(object = all_result_df, file = file.path(analysis_save_dir, "03-04-differentially-expressed-events-filepaths.rds"))