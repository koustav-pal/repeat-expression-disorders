library("here")
message("Project root:", here(), "\n")
project_dir <- here()
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
analysis_dir <- file.path(project_dir, "final-figures", "figure-5")
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
feature_ranges <- readRDS(here("05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds"))
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



lapply(names(all_contrasts_list$contrasts), function(x){
    temp_df <- all_contrasts_list$contrasts[[x]]
    temp_df$direction <- "down"
    temp_df$direction[temp_df$log2FoldChange > 0] <- "up"
    de_rci_df <- temp_df[temp_df$padj < 0.05 & !is.na(temp_df$padj),]
    data.frame(contrast = x, 
        frequency = nrow(de_rci_df),
        prop = nrow(de_rci_df) / nrow(temp_df))
})


de_intron_ranges <- feature_ranges[feature_ranges$type == "intron"]
event_frequency_df <- do.call(rbind, lapply(names(all_contrasts_list$contrasts), function(x){
    temp_df <- all_contrasts_list$contrasts[[x]]
    temp_df$direction <- "down"
    temp_df$direction[temp_df$log2FoldChange > 0] <- "up"
    de_df <- temp_df[temp_df$padj < 0.05 & !is.na(temp_df$padj) & rownames(temp_df) %in% de_intron_ranges$hash_id,]
    dir_freq <- table(de_df$direction)
    data.frame(contrast = x, 
        frequency = nrow(de_df))
}))

event_frequency_df$contrast[event_frequency_df$contrast == "sfpqkd_vs_scramble"] <- "sfpqkd"
event_frequency_df$contrast[event_frequency_df$contrast == "fuskd_vs_scramble"] <- "fuskd"
event_frequency_df$contrast[event_frequency_df$contrast == "tdp43kd_vs_scramble"] <- "tdp43kd"
event_frequency_df$contrast <- factor(event_frequency_df$contrast, levels = c("sfpqkd", "fuskd", "tdp43kd"))


current_theme <- theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))

# "AnswerALS" = "#51BBFE"
colours <- c("RCI" = "#F24236", "non-RCI" = "#F7F7F2")
# colours <- c("sfpqkd" = "#E15554" , "fuskd" = "#E1BC29", "tdp43kd" = "#7768AE")
ThePlot <- ggplot(event_frequency_df, 
    aes(x = contrast, y = prop))
ThePlot <- ThePlot + geom_col(aes(fill = rci_status), colour = "#000000")
ThePlot <- ThePlot + geom_text(aes(label = frequency), colour = "#000000")
# ThePlot <- ThePlot + stat_summary(aes(colour = condition))
ThePlot <- ThePlot + scale_x_discrete("RBP knockdown vs scramble")
ThePlot <- ThePlot + scale_y_continuous("proportion of RCIs captures upon RBP knockdown")
ThePlot <- ThePlot + scale_fill_manual("RBP", values = colours)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("Proportion of RCIs captured")
ThePlot <- ThePlot + current_theme
ggsave(width = 5, height = 10, file = file.path(analysis_plot_dir, "01-04-01-proportion-of-de-introns-captured-by-rbp-knockdown-rci-vs-non-rci.pdf"), ThePlot)



de_intron_ranges <- feature_ranges[feature_ranges$feature_id %in% sig_result_df$feature_ids & feature_ranges$type == "intron"]
event_frequency_df <- do.call(rbind, lapply(names(all_contrasts_list$contrasts), function(x){
    temp_df <- all_contrasts_list$contrasts[[x]]
    temp_df$direction <- "down"
    temp_df$direction[temp_df$log2FoldChange > 0] <- "up"
    temp_df$rci_status <- ifelse(rownames(temp_df) %in% repeat_containing_introns_df$hash_id, "RCI", "non-RCI")
    de_df <- temp_df[temp_df$padj < 0.05 & !is.na(temp_df$padj) & rownames(temp_df) %in% de_intron_ranges$hash_id,]
    dir_freq <- table(de_df$direction)
    data.frame(contrast = x, 
        rci_status = c("RCI", "non-RCI"),
        frequency = c(nrow(de_df[de_df$rci_status == "RCI",]), nrow(de_df[de_df$rci_status == "non-RCI",])),
        prop = c(nrow(de_df[de_df$rci_status == "RCI",]), nrow(de_df[de_df$rci_status == "non-RCI",])) / nrow(de_df))
}))
# & 
event_frequency_df$contrast[event_frequency_df$contrast == "sfpqkd_vs_scramble"] <- "sfpqkd"
event_frequency_df$contrast[event_frequency_df$contrast == "fuskd_vs_scramble"] <- "fuskd"
event_frequency_df$contrast[event_frequency_df$contrast == "tdp43kd_vs_scramble"] <- "tdp43kd"
event_frequency_df$contrast <- factor(event_frequency_df$contrast, levels = c("sfpqkd", "fuskd", "tdp43kd"))


current_theme <- theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))

# "AnswerALS" = "#51BBFE"
colours <- c("RCI" = "#F24236", "non-RCI" = "#F7F7F2")
# colours <- c("sfpqkd" = "#E15554" , "fuskd" = "#E1BC29", "tdp43kd" = "#7768AE")
ThePlot <- ggplot(event_frequency_df, 
    aes(x = contrast, y = prop))
ThePlot <- ThePlot + geom_col(aes(fill = rci_status), colour = "#000000")
ThePlot <- ThePlot + geom_text(aes(label = frequency), colour = "#000000")
# ThePlot <- ThePlot + stat_summary(aes(colour = condition))
ThePlot <- ThePlot + scale_x_discrete("RBP knockdown vs scramble")
ThePlot <- ThePlot + scale_y_continuous("proportion of RCIs captures upon RBP knockdown")
ThePlot <- ThePlot + scale_fill_manual("RBP", values = colours)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("Proportion of RCIs captured")
ThePlot <- ThePlot + current_theme
ggsave(width = 5, height = 10, file = file.path(analysis_plot_dir, "01-04-01-proportion-of-de-introns-captured-by-rbp-knockdown-rci-vs-non-rci.pdf"), ThePlot)



event_frequency_df <- do.call(rbind, lapply(names(all_contrasts_list$contrasts), function(x){
    temp_df <- all_contrasts_list$contrasts[[x]]
    temp_df$direction <- "down"
    temp_df$direction[temp_df$log2FoldChange > 0] <- "up"
    de_rci_df <- temp_df[temp_df$padj < 0.05 & !is.na(temp_df$padj) & rownames(temp_df) %in% repeat_containing_introns_df$hash_id,]
    dir_freq <- table(de_rci_df$direction)
    data.frame(contrast = x, 
        frequency = nrow(de_rci_df),
        prop = nrow(de_rci_df) / rci_length)
}))
event_frequency_df$contrast[event_frequency_df$contrast == "sfpqkd_vs_scramble"] <- "sfpqkd"
event_frequency_df$contrast[event_frequency_df$contrast == "fuskd_vs_scramble"] <- "fuskd"
event_frequency_df$contrast[event_frequency_df$contrast == "tdp43kd_vs_scramble"] <- "tdp43kd"
event_frequency_df$contrast <- factor(event_frequency_df$contrast, levels = c("sfpqkd", "fuskd", "tdp43kd"))


current_theme <- theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))

# "AnswerALS" = "#51BBFE"
colours <- c("up" = "#F24236", "down" = "#F7F7F2")
colours <- c("sfpqkd" = "#E15554" , "fuskd" = "#E1BC29", "tdp43kd" = "#7768AE")
ThePlot <- ggplot(event_frequency_df, 
    aes(x = contrast, y = prop))
ThePlot <- ThePlot + geom_col(aes(fill = contrast), colour = "#000000")
ThePlot <- ThePlot + geom_text(aes(label = frequency), colour = "#000000")
# ThePlot <- ThePlot + stat_summary(aes(colour = condition))
ThePlot <- ThePlot + scale_x_discrete("RBP knockdown vs scramble")
ThePlot <- ThePlot + scale_y_continuous("proportion of RCIs captures upon RBP knockdown")
ThePlot <- ThePlot + scale_fill_manual("RBP", values = colours)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("Proportion of RCIs captured")
ThePlot <- ThePlot + current_theme
ggsave(width = 5, height = 10, file = file.path(analysis_plot_dir, "01-04-01-proportion-of-rcis-captured-by-rbp-knockdown.pdf"), ThePlot)




direction_frequency_df <- do.call(rbind, lapply(names(all_contrasts_list), function(x){
    temp_df <- all_contrasts_list[[x]]
    temp_df$direction <- "down"
    temp_df$direction[temp_df$log2FoldChange > 0] <- "up"
    de_rci_df <- temp_df[temp_df$padj < 0.05 & !is.na(temp_df$padj) & rownames(temp_df) %in% repeat_containing_introns_df$hash_id,]
    dir_freq <- table(de_rci_df$direction)
    data.frame(contrast = x, 
        direction = c("up", "down"),
        frequency = c(dir_freq["up"], dir_freq["down"]),
        prop = c(dir_freq["up"], dir_freq["down"]) / sum(dir_freq))
}))
direction_frequency_df$contrast[direction_frequency_df$contrast == "sfpqkd_vs_scramble"] <- "sfpqkd"
direction_frequency_df$contrast[direction_frequency_df$contrast == "fuskd_vs_scramble"] <- "fuskd"
direction_frequency_df$contrast[direction_frequency_df$contrast == "tdp43kd_vs_scramble"] <- "tdp43kd"
direction_frequency_df$contrast <- factor(direction_frequency_df$contrast, levels = c("sfpqkd", "fuskd", "tdp43kd"))



# "AnswerALS" = "#51BBFE"
colours <- c("up" = "#F24236", "down" = "#F7F7F2")
ThePlot <- ggplot(direction_frequency_df, 
    aes(x = contrast, y = frequency))
ThePlot <- ThePlot + geom_col(aes(fill = direction), colour = "#000000", position = position_dodge(width = 0.9))
ThePlot <- ThePlot + geom_text(aes(label = frequency), colour = "#000000")
# ThePlot <- ThePlot + stat_summary(aes(colour = condition))
ThePlot <- ThePlot + scale_x_discrete("Proportion of events up/down-regulated")
ThePlot <- ThePlot + scale_y_continuous("proportion of DE RCIs up or down")
ThePlot <- ThePlot + scale_fill_manual("direction", values = colours)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("Proportion of RCIs captured")
ThePlot <- ThePlot + current_theme
ggsave(width = 5, height = 10, file = file.path(analysis_plot_dir, "01-04-02-proportion-of-rcis-up-or-downregulated.pdf"), ThePlot)
