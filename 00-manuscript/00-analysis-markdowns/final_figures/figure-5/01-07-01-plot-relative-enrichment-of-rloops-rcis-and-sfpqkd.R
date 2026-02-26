library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- file.path(project_dir, "final-figures", "figure-5")
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
feature_ranges <- readRDS(here("05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds"))
samples_param_df <- read.csv(here("40-giulia-rbp-kd/input_files/01-01-public-data-params-dataframe.txt"))
rloop_stats_df <- readRDS(here("30-sequence-analysis-of-introns/rdata_files/09-02-rloop-stats-dataframe.rds"))
sfpqkd_de_intron_df <- readRDS(here("40-giulia-rbp-kd/rdata_files/04-01-save-introns-differentially-expressed-in-als-vs-control-and-sfpqkd.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
require("DESeq2")
require("ggplot2")

# Get VCP specific differentially expressed rcis
rloop_stats_df$hash_id <- feature_ranges$hash_id[match(rloop_stats_df$feature_id, feature_ranges$feature_id)]

sfpqkd_de_intron_df$rloop_frequency <- rloop_stats_df$rloop_frequency[match(rownames(sfpqkd_de_intron_df), rloop_stats_df$hash_id)]
sfpqkd_de_intron_df$rloop_status <- ifelse(!is.na(sfpqkd_de_intron_df$rloop_frequency) & sfpqkd_de_intron_df$rloop_frequency > 0, "R-loop containing", "No R-loops")

# intersect sfpq kd associated RCIs and those with R-loops

chisq_results_1 <- chisq.test(table(sfpqkd_de_intron_df$direction, sfpqkd_de_intron_df$rloop_status))
#         Pearson's Chi-squared test with Yates' continuity correction

# data:  table(sfpqkd_de_intron_df$direction, sfpqkd_de_intron_df$rloop_status)
# X-squared = 818.67, df = 1, p-value < 2.2e-16
# chisq_results_1$observed / chisq_results_1$expected 
  #      No R-loops R-loop containing
  # down  1.2023565         0.6887626
  # up    0.8835173         1.1791578

chisq_results_2 <- chisq.test(table(sfpqkd_de_intron_df$direction, sfpqkd_de_intron_df$rci_status))
#         Pearson's Chi-squared test with Yates' continuity correction

# data:  table(sfpqkd_de_intron_df$direction, sfpqkd_de_intron_df$rci_status)
# X-squared = 39.765, df = 1, p-value = 2.865e-10
# chisq_results_2$observed / chisq_results_2$expected
  #        non-RCI       RCI
  # down 1.0135635 0.7725343
  # up   0.9921924 1.1309363

chisq_results_3 <- chisq.test(table(sfpqkd_de_intron_df$rloop_status, sfpqkd_de_intron_df$rci_status))
#         Pearson's Chi-squared test with Yates' continuity correction

# data:  table(sfpqkd_de_intron_df$rloop_status, sfpqkd_de_intron_df$rci_status)
# X-squared = 726.54, df = 1, p-value < 2.2e-16
# chisq_results_3$observed / chisq_results_3$expected
  #                     non-RCI       RCI
  # No R-loops        1.0353388 0.4073528
  # R-loop containing 0.9456467 1.9115295





frequency_table <- xtabs(~ direction + rci_status + rloop_status, data = sfpqkd_de_intron_df)

frequency_table_df <- as.data.frame(frequency_table)

full_model <- glm(Freq ~ direction * rci_status * rloop_status, family = poisson, data = frequency_table_df)
summary(full_model)

reduced_model <- glm(Freq ~ direction + rci_status + rloop_status, family = poisson, data = frequency_table_df)
summary(reduced_model)


frequency_table_df$fitted <- predict(reduced_model, type = "response")
frequency_table_df$obs_over_exp <- frequency_table_df$Freq / frequency_table_df$fitted
frequency_table_df$log_foldchange <- log2(frequency_table_df$obs_over_exp)
MinMax <- min(abs(c(min(frequency_table_df$log_foldchange), max(frequency_table_df$log_foldchange))))
frequency_table_df$label_fc <- frequency_table_df$log_foldchange
frequency_table_df$label_fc[frequency_table_df$label_fc < MinMax * -1] <- MinMax * -1 
frequency_table_df$label_fc[frequency_table_df$label_fc > MinMax] <- MinMax


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


frequency_table_df$rloop_status <- factor(frequency_table_df$rloop_status, levels = c("R-loop containing", "No R-loops"))
frequency_table_df$rci_status <- factor(frequency_table_df$rci_status, levels = c("RCI", "non-RCI"))
frequency_table_df$direction <- factor(frequency_table_df$direction, levels = c("down", "up"))
colours <- c("Same" = "#e41a1c", "Opposite" = "#377eb8")
ThePlot <- ggplot(frequency_table_df, 
    aes(x = rci_status, y = direction))
ThePlot <- ThePlot + geom_tile(aes(fill = label_fc))
ThePlot <- ThePlot + geom_text(aes(label = round(log_foldchange,2)))
ThePlot <- ThePlot + scale_x_discrete("SFPQKD direction")
ThePlot <- ThePlot + scale_y_discrete("RCI status")
# ThePlot <- ThePlot + scale_colour_manual("Condition colour", values = colours)
ThePlot <- ThePlot + scale_fill_distiller(palette = "RdBu", direction = -1)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_grid(. ~ rloop_status)
ThePlot <- ThePlot + ggtitle("Enrichment of SFPQKD and R-loop containing introns stratified by RCIs")
ThePlot <- ThePlot + current_theme
ggsave(width = 7, height = 7, file = file.path(analysis_plot_dir, "03-01-01-plot-enirchment-of-sfpqkd-rloops-and-rcis.pdf"), ThePlot)


temp_id_1 <- paste(sfpqkd_de_intron_df$direction, sfpqkd_de_intron_df$rloop_status, sfpqkd_de_intron_df$rci_status, sep = "_")
temp_id_2 <- paste(frequency_table_df$direction, frequency_table_df$rloop_status, frequency_table_df$rci_status, sep = "_")
sfpqkd_de_intron_df$log2_enrichment <- frequency_table_df$log_foldchange[match(temp_id_1, temp_id_2)]


saveRDS(object = sfpqkd_de_intron_df, file = file.path(analysis_save_dir, "04-03-save-relative-enrichment-of-different-intronic-groups.rds"))