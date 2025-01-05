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
object_list <- readRDS(file.path(analysis_save_dir, "03-02-list_of_counts_stat_samplesheet_from_nygc.rds"))
isPairedEnd <- T
strand_specificity <- 2
sampleinfo_df <- read.csv(file.path(analysis_input_dir, "00_filtered_sample_gsm_list.txt"))
clinical_metadata_df <- read.csv("../0103-NYGC-consortium/input_files/samplesheet.csv")
# ==========================================================================
# Imports
# ==========================================================================
require("GenomicRanges")
require("rslurm")
require("stringr")
require("matrixStats")
# require("sva")
require("DESeq2")
# library("bladderdata")
source("../09-01-QC-compare-fractional-counts-vs-abs-counts/scripts/src/deseq_functions.R")
source("../../useful-functions/plot-utils.R")

counts_matrix <- object_list$counts

sample_df <- object_list$sample_df
sample_df <- sample_df[sample_df$sample_name %in% sampleinfo_df$Experiment,]

sample_df$mutation <- clinical_metadata_df$mutation[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$tissue <- clinical_metadata_df$sample_source[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$instrument <- sampleinfo_df$Instrument[match(sample_df$sample_name, sampleinfo_df$Experiment)]
sample_df$sex_genotype <- clinical_metadata_df$sex_genotype[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$rin <- clinical_metadata_df$rin[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$family_history_of_als_ftd <- clinical_metadata_df$family_history_of_als_ftd[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$mnd_with_ftd <- clinical_metadata_df$mnd_with_ftd[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$condition <- clinical_metadata_df$condition[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$condition[sample_df$mnd_with_ftd == "Yes"] <- "als+ftd"

no_low_rin_sample_df <- sample_df[sample_df$rin > 4 & !is.na(sample_df$rin) & sample_df$condition %in% c("als", "ctrl", "als+ftd"),]

low_counts_removed_matrix <- counts_matrix[rowMeans(counts_matrix) > 1,]

current_sampleinfo <- no_low_rin_sample_df[no_low_rin_sample_df$column_id %in% colnames(low_counts_removed_matrix),]
current_sampleinfo$type <- current_sampleinfo$mutation
current_sampleinfo$type[!(current_sampleinfo$type %in% c("sporadic", "ctrl", "none", "c9orf72", "tardbp"))] <- "mutation"

col_data <- unique(current_sampleinfo[,c("column_id", "mutation", "tissue", "type", "condition", "instrument", "sex_genotype")])
col_data$tissue[col_data$tissue %in% c("Cortex_Frontal", "Cortex_Temporal", "Cortex_Occipital")] <- "Cortex"
col_data$tissue[col_data$tissue %in% c("Cortex_Motor_Lateral", "Cortex_Motor_Medial", "Cortex_Motor_Unspecified")] <- "Cortex_Motor"
col_data$mutation <- factor(col_data$mutation, 
    levels = c("ctrl", "sporadic", "c9orf72", "tardbp", "dctn1", "ang", "setx", "optn", "anxa11", "nefh",
        "sod1", "ubqln2", "fus", "vcp", "none", "matr3", "nek1", "als2", "fig4", "grn", 
        "tbk1", "mapt"))
col_data$tissue <- factor(col_data$tissue, 
    levels = c("Spinal_Cord_Lumbar", "Spinal_Cord_Cervical", "Cerebellum",
        "Cortex_Motor", "Cortex", "Hippocampus", "Spinal_Cord_Thoracic"))
col_data$type <- factor(col_data$type, 
    levels = c("ctrl", "none", "sporadic", "mutation", 
        "tardbp", "c9orf72"))

col_data <- col_data[col_data$tissue != "Hippocampus",]
# model_1 <- model.matrix( ~ type + tissue, data = col_data)

selected_counts_matrix <- round(low_counts_removed_matrix)[,col_data$column_id]
selected_counts_matrix_with_pseudo <- selected_counts_matrix + 1

deseq_object <- DESeqDataSetFromMatrix(countData = selected_counts_matrix_with_pseudo, 
    colData = col_data, 
    design = ~ sex_genotype + instrument + tissue + condition)
differential_object <- DESeq(deseq_object)

saveRDS(object = differential_object, file = file.path(analysis_save_dir, "03-03-save-deseq2-differential-object.rds"))