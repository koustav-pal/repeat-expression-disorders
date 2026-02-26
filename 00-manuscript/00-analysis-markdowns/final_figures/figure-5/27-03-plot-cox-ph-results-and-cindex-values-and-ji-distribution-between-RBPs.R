library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- here("final-figures", "figure-5")
save_dir <- file.path(analysis_dir, "rdata_files")
input_dir <- file.path(analysis_dir,"input_files")
plot_dir <- file.path(analysis_dir, "plots")
output_dir <- file.path(analysis_dir, "output_files")
model_validation_list <- readRDS(file.path(save_dir, "07-02-save-all-ridge-coefficient-models-from-answerals.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
# source("../../src/plot-utils.R")
# source("../../src/basepair_coverage_functions.R")
source("src/colours.R")
require("HiCBricks")
require("GenomicRanges")
require("stringr")
require("ggplot2")
library(rlang)
library(riskRegression)
library(survival)
library(survminer)
library(transformr)
library(plotly)
library(reshape2)
library(knitr)
library(rmarkdown)
require("RColorBrewer")
require("randomForestSRC")
require("parallel")
require("matrixStats")
require("survival")
require("patchwork")
require("ggplot2")
library(broom)
library(ComplexUpset)
library(ggplot2)

# ==========================================================================
# Function
# ==========================================================================

cindex_df <- model_validation_list$cindex_list
cox_result_df <- model_validation_list$coxph_results
selected_rbp_enriched_rcis_df = model_validation_list$selected_rbp_enriched_rcis


rci_n_df <- do.call(rbind, lapply(names(selected_rbp_enriched_rcis_df), function(x){
	id_split <- str_split(x, pattern = " ", simplify = T)
	data.frame(rbp = id_split[1], repeat_class = id_split[2], n = length(unique(selected_rbp_enriched_rcis_df[[x]])))
}))
#        rbp    repeat_class   n
# 1     AATF  quadnucleotide  29
# 2    CPEB4   trinucleotide 217
# 3    CPSF6   trinucleotide 217
# 4    CSTF2 pentanucleotide   8
# 5    DDX21 pentanucleotide  11
# 6     DDX6  hexanucleotide 150
# 7   DROSHA  hexanucleotide 166
# 8    EIF3D   trinucleotide 206
# 9  FAM120A  hexanucleotide 160
# 10    GARS pentanucleotide   9
# 11  GEMIN5 pentanucleotide   6
# 12  GEMIN5   trinucleotide 124
# 13   GPKOW   trinucleotide 193
# 14  GTF2F1  hexanucleotide 145
# 15  METTL1  hexanucleotide 101
# 16  METTL1 pentanucleotide   4
# 17  METTL1   trinucleotide 143
# 18   MORC2 pentanucleotide  10
# 19   MTPAP  quadnucleotide  25
# 20   MTPAP   trinucleotide 218
# 21    NIP7  hexanucleotide 144
# 22    NIP7 pentanucleotide  10
# 23   NIPBL  quadnucleotide  30
# 24   NOL12  hexanucleotide 135
# 25   NOL12 pentanucleotide  14
# 26    NPM1  quadnucleotide  24
# 27    NPM1   trinucleotide 209
# 28   PCBP1  hexanucleotide 151
# 29   PCBP1  quadnucleotide  25
# 30   PCBP2  hexanucleotide 145
# 31   PCBP2  quadnucleotide  24
# 32    PUM1   trinucleotide 201
# 33   RBM22  quadnucleotide  24
# 34   RPS11  hexanucleotide 134
# 35   RPS11 pentanucleotide  11
# 36    RPS6  quadnucleotide  28
# 37   SF3B1  quadnucleotide  24
# 38   STAU2  quadnucleotide  20
# 39  TROVE2   trinucleotide 209
# 40   ZC3H8 pentanucleotide  10

rbp_max_lrts <- vapply(split(cox_result_df, cox_result_df$rbp), function(temp_df){
	max(temp_df$LRT_ChiSq)
},1)
rbp_max_lrts_df <- data.frame(rbp = names(rbp_max_lrts), group = "lrt", value = rbp_max_lrts)
rbp_max_lrts_df$rbp <- factor(rbp_max_lrts_df$rbp, levels = rbp_max_lrts_df$rbp[order(rbp_max_lrts_df$value, decreasing = T)])

cox_result_df$rbp <- factor(cox_result_df$rbp, levels = names(rbp_max_lrts[order(rbp_max_lrts, decreasing = T)]))
sig_cox_result_df <- cox_result_df[cox_result_df$p_LRT < 0.05 & cox_result_df$p_wald < 0.05,]
sig_cox_result_df$rbp <- droplevels(sig_cox_result_df$rbp)

rbp_max_lrts_df <- rbp_max_lrts_df[rbp_max_lrts_df$rbp %in% sig_cox_result_df$rbp,]
rbp_max_lrts_df$rbp <- droplevels(rbp_max_lrts_df$rbp)

rci_n_df <- rci_n_df[(rci_n_df$rbp %in% sig_cox_result_df$rbp) & (rci_n_df$repeat_class %in% sig_cox_result_df$repeat_class),]
rci_n_df$rbp <- factor(rci_n_df$rbp, levels = levels(rbp_max_lrts_df$rbp))

cindex_df <- cindex_df[cindex_df$rbp %in% c("baseline", levels(rbp_max_lrts_df$rbp)),]
cindex_df$rbp <- factor(cindex_df$rbp, levels = c("baseline", levels(rbp_max_lrts_df$rbp)))

max(cindex_df$C[cindex_df$rbp != "baseline"])
# percent_change_over_baseline <- (cindex_df$C[cindex_df$rbp != "baseline"] - cindex_df$C[cindex_df$rbp == "baseline"]) / cindex_df$C[cindex_df$rbp == "baseline"]
median(sig_cox_result_df$HR[sig_cox_result_df$repeat_class == "hexanucleotide"])



current_theme <- theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))

ThePlot <- ggplot(sig_cox_result_df, aes(x = rbp, y = HR))
ThePlot <- ThePlot + geom_pointrange(aes(ymin = HR_lower, ymax = HR_upper, colour = repeat_class), position = position_dodge(width = 0.9))
ThePlot <- ThePlot + geom_hline(yintercept = 1, colour = "#000000",linetype = "dashed")
ThePlot <- ThePlot + scale_x_discrete("RNA binding proteins")
ThePlot <- ThePlot + scale_y_continuous("Hazard ratio")
ThePlot <- ThePlot + scale_colour_manual("Repeat classes", values = repeat_colour_groups)
# ThePlot <- ThePlot + scale_colour_manual("Intron groups", values = Colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("Hazard ratio of ")
ThePlot <- ThePlot + current_theme

ThePlot_2 <- ggplot(rbp_max_lrts_df, aes(x = rbp, y = group))
ThePlot_2 <- ThePlot_2 + geom_tile(aes(fill = value))
ThePlot_2 <- ThePlot_2 + scale_x_discrete("RNA binding proteins")
ThePlot_2 <- ThePlot_2 + scale_y_discrete("LRT")
ThePlot_2 <- ThePlot_2 + scale_fill_distiller("LRT", palette = "Reds", direction = 1)
# ThePlot <- ThePlot + scale_colour_manual("Intron groups", values = Colours)
ThePlot_2 <- ThePlot_2 + theme_bw()
ThePlot_2 <- ThePlot_2 + ggtitle("LRTs of ANOVA comparing baseline and incorporationg RBP bound repeats")
ThePlot_2 <- ThePlot_2 + current_theme


ThePlot_3 <- ggplot(rci_n_df, aes(x = rbp, y = n))
ThePlot_3 <- ThePlot_3 + geom_col(aes(fill = repeat_class), colour = "#000000", position = position_dodge(width = 0.9))
ThePlot_3 <- ThePlot_3 + scale_x_discrete("RNA binding proteins")
ThePlot_3 <- ThePlot_3 + scale_y_continuous("n of RCIs")
ThePlot_3 <- ThePlot_3 + scale_fill_manual("Repeat classes", values = repeat_colour_groups)
# ThePlot <- ThePlot + scale_colour_manual("Intron groups", values = Colours)
ThePlot_3 <- ThePlot_3 + theme_bw()
ThePlot_3 <- ThePlot_3 + ggtitle("n of RCIs in each RBP class")
ThePlot_3 <- ThePlot_3 + current_theme


ThePlot_4 <- ggplot(cindex_df, aes(x = rbp, y = C))
ThePlot_4 <- ThePlot_4 + geom_pointrange(aes(ymin = C_lower, ymax = C_upper, colour = repeat_class), position = position_dodge(width = 0.9))
ThePlot_4 <- ThePlot_4 + scale_x_discrete("RNA binding proteins")
ThePlot_4 <- ThePlot_4 + scale_y_continuous("Bootstrap corrected C index")
ThePlot_4 <- ThePlot_4 + scale_colour_manual("Repeat classes", values = c("none" = "#2A2B2A", repeat_colour_groups))
# ThePlot <- ThePlot + scale_colour_manual("Intron groups", values = Colours)
ThePlot_4 <- ThePlot_4 + theme_bw()
ThePlot_4 <- ThePlot_4 + ggtitle("optimism corrected C index")
ThePlot_4 <- ThePlot_4 + current_theme


all_plots <- ThePlot_3 / ThePlot_2 / ThePlot / ThePlot_4

ggsave(all_plots, file = file.path(plot_dir, "09-AnswerALS-hazard-ratio-of-each-RBP-enriched-repeats-region.pdf"), width = 10, height = 10)
# ggsave(ThePlot, file = file.path(plot_dir, "01-hazard-ratio-of-each-RBP-enriched-repeats-region.pdf"), width = 10, height = 5)

# ==========================================================================
# Plot Upset of RBPs
# ==========================================================================


all_rcis <- as.vector(unique(do.call(c, selected_rbp_enriched_rcis_df)))

rci_rbp_df <- do.call(rbind, lapply(names(selected_rbp_enriched_rcis_df), function(x){
	id_split <- str_split(x, pattern = " ", simplify = T)
	data.frame(rbp = id_split[1], repeat_class = id_split[2], rci = selected_rbp_enriched_rcis_df[[x]])
}))

rci_rbp_split <- split(rci_rbp_df$rci, rci_rbp_df$rbp)


pairs_matrix <- cbind(rep(names(rci_rbp_split), each = length(rci_rbp_split)), rep(names(rci_rbp_split), times = length(rci_rbp_split)))


all_ji_df <- do.call(rbind, lapply(seq_len(nrow(pairs_matrix)), function(x){
	current_x <- pairs_matrix[x,1]
	current_y <- pairs_matrix[x,2]
	jaccard <- length(intersect(rci_rbp_split[[current_x]], rci_rbp_split[[current_y]])) / length(union(rci_rbp_split[[current_x]], rci_rbp_split[[current_y]]))
	data.frame(x = current_x, y = current_y, ji = jaccard)
}))

rbps_by_rci <- lapply(split(rci_rbp_df$rbp, rci_rbp_df$repeat_class), unique)



rbp_interconnectivity_df <- do.call(rbind, lapply(names(rbps_by_rci), function(repeat_class){
	x <- rbps_by_rci[[repeat_class]]
	current_rbp_df <- all_ji_df[all_ji_df$x %in% x & 
		all_ji_df$y %in% x & 
		all_ji_df$x != all_ji_df$y,]
	if(nrow(current_rbp_df) == 0){
		return(NULL)
	}
	median_ji <- vapply(split(current_rbp_df$ji, current_rbp_df$x), median, 1)
	current_rbp_df$repeat_class <- repeat_class
	current_rbp_df$x <- factor(current_rbp_df$x, levels = names(median_ji[order(median_ji, decreasing = T)]))
	return(current_rbp_df)
}))


ThePlot_5 <- ggplot(rbp_interconnectivity_df, aes(x = x, y = ji))
# ThePlot_5 <- ThePlot_5 + geom_violin(aes(fill = repeat_class), colour = "#000000")
ThePlot_5 <- ThePlot_5 + geom_boxplot(aes(fill = repeat_class), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
ThePlot_5 <- ThePlot_5 + scale_x_discrete("RNA binding proteins")
ThePlot_5 <- ThePlot_5 + scale_y_continuous("jaccard index")
ThePlot_5 <- ThePlot_5 + scale_fill_manual("Repeat classes", values = repeat_colour_groups)
ThePlot_5 <- ThePlot_5 + theme_bw()
ThePlot_5 <- ThePlot_5 + ggtitle("optimism corrected C index")
ThePlot_5 <- ThePlot_5 + facet_wrap(repeat_class ~ ., scales = "free", nrow = 3, ncol = 1)
ThePlot_5 <- ThePlot_5 + current_theme

ggsave(ThePlot_5, file = file.path(plot_dir, "02-jaccard-index-distribution-for-each-rbps.pdf"), width = 10, height = 10)
