data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "40-giulia-rbp-kd")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
samplesheet_df <- read.csv("/nemo/lab/patanir/home/users/ziffo/patani-collab/motor-neuron-rbp-knockdown/sample-details/samplesheet.csv")
samples_param_df <- read.csv(file.path(analysis_input_dir, "01-01-public-data-params-dataframe.txt"))
sig_result_list <- readRDS(file.path("../20-import-significant-events/rdata_files/02-significantly-differentially-expressed-events.rds"))
all_contrasts_list <- readRDS(file.path(analysis_save_dir, "01-03-list-of-contrasts-for-sfpqkd-fuskd-and-tdp43kd.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
require("DESeq2")
require("ggplot2")


feature_ranges_with_repeats <- sig_result_list$repeat_containing_ranges
sig_result_df <- sig_result_list$significant_events

repeat_containing_introns_df <- sig_result_df[sig_result_df$status == "repeat_containing" & 
    sig_result_df$region == "intron",]
repeat_containing_introns_df$hash_id <- feature_ranges$hash_id[match(repeat_containing_introns_df$feature_ids, feature_ranges$feature_id)]
rci_length <- length(unique(repeat_containing_introns_df$hash_id))

de_intron_ranges <- feature_ranges[feature_ranges$feature_id %in% sig_result_df$feature_ids & 
	feature_ranges$type == "intron"]


de_intron_df <- do.call(rbind, lapply(names(all_contrasts_list$contrasts), function(x){
    temp_df <- all_contrasts_list$contrasts[[x]]
    temp_df$direction <- "down"
    temp_df$direction[temp_df$log2FoldChange > 0] <- "up"
    de_df <- temp_df[temp_df$padj < 0.05 & !is.na(temp_df$padj) & 
    	rownames(temp_df) %in% de_intron_ranges$hash_id,]
    de_df$contrast <- x
    return(de_df)
}))
de_intron_df$contrast <- str_replace(de_intron_df$contrast, pattern = "_vs_scramble", replacement = "")

sfpqkd_de_intron_df <- de_intron_df[de_intron_df$contrast == "sfpqkd",]
sfpqkd_de_intron_df$hash_id <- rownames(sfpqkd_de_intron_df)
rownames(sfpqkd_de_intron_df) <- NULL

sfpqkd_de_intron_df$rci_status <- ifelse(sfpqkd_de_intron_df$hash_id %in% repeat_containing_introns_df$hash_id, "RCI", "non-RCI")
# non-RCI     RCI
#   21332    1272
#   

saveRDS(object = sfpqkd_de_intron_df, file = file.path(analysis_save_dir, "04-01-save-introns-differentially-expressed-in-als-vs-control-and-sfpqkd.rds"))