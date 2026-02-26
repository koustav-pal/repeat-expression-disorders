library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- file.path(project_dir, "final-figures", "figure-5")
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
feature_ranges <- readRDS(here("05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds"))
samplesheet_df <- read.csv("/nemo/lab/patanir/home/users/ziffo/patani-collab/motor-neuron-rbp-knockdown/sample-details/samplesheet.csv")
samples_param_df <- read.csv(here("40-giulia-rbp-kd/input_files/01-01-public-data-params-dataframe.txt"))
results_list <- readRDS(here("40-giulia-rbp-kd/rdata_files/01-02-list_of_counts_stat_samplesheet_from_answerals.rds"))
sig_result_list <- readRDS(here("20-import-significant-events/rdata_files/02-significantly-differentially-expressed-events.rds"))
all_contrasts_list <- readRDS(here("40-giulia-rbp-kd/rdata_files/01-03-list-of-contrasts-for-sfpqkd-fuskd-and-tdp43kd.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
require("DESeq2")
require("ggplot2")

sample_df <- results_list$sample_df
counts_matrix <- results_list$counts
min_count <- 1

rci_ranges <- sig_result_list$repeat_containing_ranges

sig_result_df <- sig_result_list$significant_events
repeat_containing_introns_df <- sig_result_df[sig_result_df$status == "repeat_containing" & 
    sig_result_df$region == "intron",]
repeat_containing_introns_df$hash_id <- feature_ranges$hash_id[match(repeat_containing_introns_df$feature_ids, feature_ranges$feature_id)]
rci_length <- length(unique(repeat_containing_introns_df$hash_id))



contrast_df <- do.call(rbind, lapply(names(all_contrasts_list$contrasts), function(x){
    temp_df <- all_contrasts_list$contrasts[[x]]
    temp_df$rci <- rownames(temp_df)
    temp_df$rbp <- x
    de_rci_df <- temp_df[temp_df$padj < 0.05 & !is.na(temp_df$padj) & rownames(temp_df) %in% repeat_containing_introns_df$hash_id,]
    return(de_rci_df)
}))

contrast_df$rbp <- str_split(string = contrast_df$rbp, pattern = "_", simplify = T)[,1]
rci_in_rbp_kd <- split(contrast_df$rci, contrast_df$rbp)

# colours <- c("#F6F4F3", "#F87575", "#FFA9A3")
require("ggvenn")
colours <- c("#E15554", "#E1BC29", "#7768AE")
pdf(file.path(analysis_plot_dir, "02-07-intersection-between-sfpq-kd-fus-kd-and-tardbp-kd.pdf"), width = 7, height = 7)
ggvenn(rci_in_rbp_kd[c("sfpqkd","fuskd","tdp43kd")], fill_color = colours)
dev.off()
