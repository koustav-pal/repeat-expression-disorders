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
source(file.path(analysis_dir, "scripts", "src", "helper_functions.R"))


sample_df <- results_list$sample_df
counts_matrix <- results_list$counts

all_functions_list <- return_function_list()

save_paths_df <- readRDS(file.path(analysis_save_dir, "02-01-study-specific-differential-expression-profiling-paths.rds"))
all_differential_objects_list <-  lapply(seq_len(nrow(save_paths_df)), function(x){
    sample_name <- save_paths_df$sample_name[x]
    message(sample_name)
    readRDS(save_paths_df$save_path[x])
})
names(all_differential_objects_list) <- save_paths_df$sample_name


contrast_df <- read.table(file.path(analysis_input_dir, "contrast_table.txt"), 
    sep = "\t", header = T)

contrast_df_split <- split(contrast_df, contrast_df$molecule)


all_results_df <- do.call(rbind, lapply(names(all_differential_objects_list), function(x){
    contrast_list <- build_contrasts(all_differential_objects_list[[x]], contrast_df_split[[x]])
    result_df <- do.call(rbind, lapply(names(contrast_list), function(y){
        current_results <- results(object = all_differential_objects_list[[x]], contrast = contrast_list[[y]])
        filter <- !is.na(current_results$pvalue)
        return_df <- data.frame(study = x,
            contrasts = y, 
            feature_ids = rownames(current_results)[filter],
            logfc = current_results$log2FoldChange[filter],
            pval = current_results$pvalue[filter],
            padj = current_results$padj[filter])
        return(return_df)
    }))
    return(result_df)
}))

saveRDS(object = all_results_df, file = file.path(analysis_save_dir, "03-01-01-assembled-contrasts-from-all-study-specific-pairwise-comparisons.rds"))


pan_als_differential_expression_object <- readRDS(file.path(analysis_save_dir, "02-02-01-save-deseq2-all-samples-als-vs-control-differential-object.rds"))
contrast_df <- data.frame(molecule = "pan-als", from_var = "als_group",
    from_value = "als", to_var = "als_group",
    to_value = "ctrl", contrast_name = "diff_expr_als_vs_ctrl",
    group = 1)

contrast_list <- build_contrasts(pan_als_differential_expression_object, contrast_df)
pan_als_results <- results(object = pan_als_differential_expression_object, contrast = contrast_list[[1]])
filter <- !is.na(pan_als_results$pvalue)
pan_als_all_df <- data.frame(study = "pan-als",
    contrasts = "diff_expr_als_vs_ctrl", 
    feature_ids = rownames(pan_als_results)[filter],
    logfc = pan_als_results$log2FoldChange[filter],
    pval = pan_als_results$pvalue[filter],
    padj = pan_als_results$padj[filter])

saveRDS(object = pan_als_all_df, file = file.path(analysis_save_dir, "03-01-02-assembled-contrasts-from-pan-als-comparison.rds"))
