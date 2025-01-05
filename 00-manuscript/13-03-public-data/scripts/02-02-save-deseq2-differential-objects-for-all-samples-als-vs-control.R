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


all_experiment_df <- do.call(rbind, lapply(seq_along(all_functions_list), function(x){
    sample_name <- names(all_functions_list)[x]
    experiment_df_list <- all_functions_list[[sample_name]](sample_df, study_specific_design = F)
    experiment_df_list$experiment
}))

all_experiment_df$condition[all_experiment_df$condition == "vcp mutant"] <- "vcp"
all_experiment_df$condition[all_experiment_df$condition == "vcp knock-in"] <- "iso"
all_experiment_df$condition[all_experiment_df$condition == "tardbp mutant"] <- "tardbp"
all_experiment_df$condition[all_experiment_df$condition == "wildtype"] <- "ctrl"
all_experiment_df$als_group <- all_experiment_df$condition
all_experiment_df$als_group[all_experiment_df$condition == "ctrl"] <- "ctrl"
all_experiment_df$als_group[!(all_experiment_df$condition %in% c("ctrl", "iso", "tdp43kd"))] <- "als"
all_experiment_df$als_group <- factor(all_experiment_df$als_group, levels = c("ctrl", "iso", "als", "tdp43kd"))


whole_cell_experiment_df <- all_experiment_df[all_experiment_df$fraction %in% c("who", "cytoplasm", "cyt") & all_experiment_df$als_group != "tdp43kd",]
whole_cell_experiment_df$study <- factor(whole_cell_experiment_df$study)
whole_cell_experiment_df$instrument <- factor(whole_cell_experiment_df$instrument)
whole_cell_experiment_df$library_layout <- factor(whole_cell_experiment_df$library_layout)
whole_cell_experiment_df <- whole_cell_experiment_df[whole_cell_experiment_df$als_group %in% c("ctrl", "als"),]
rownames(whole_cell_experiment_df) <- NULL

min_count <- 1
counts_matrix <- counts_matrix[rowMeans(counts_matrix) > min_count,]
save_path <- file.path(analysis_save_dir, paste("02-02-01", "save-deseq2-all-samples-als-vs-control-differential-object.rds", sep = "-"))
deseq_object <- DESeqDataSetFromMatrix(countData = round(counts_matrix[,whole_cell_experiment_df$samplename]), 
    colData = whole_cell_experiment_df, 
    design = ~ study + als_group)
differential_object <- DESeq(deseq_object)
saveRDS(object = differential_object, file = save_path)
