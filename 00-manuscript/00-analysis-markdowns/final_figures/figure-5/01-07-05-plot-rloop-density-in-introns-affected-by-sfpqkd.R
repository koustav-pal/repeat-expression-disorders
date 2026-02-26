library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- file.path(project_dir, "final-figures", "figure-5")
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
sfpqkd_de_introns_df <- readRDS(here("40-giulia-rbp-kd/rdata_files/04-05-introns-affected-by-SFPQKD-with-rloop-density.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
require("DESeq2")
require("ggplot2")
require("BSgenome.Hsapiens.UCSC.hg38")
require("Biostrings")
require("GenomicFeatures")
require("patchwork")

require(ggbeeswarm)

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

colours <- c("up" = "#F03A47", "down" = "#276FBF")
# colours <- c("CTRL" = "#F6F4F3", "VCP" = "#F87575")

sfpqkd_de_introns_df$rci_status <- factor(sfpqkd_de_introns_df$rci_status, levels = c("RCI", "non-RCI"))
sfpqkd_de_introns_df$direction <- factor(sfpqkd_de_introns_df$direction, levels = c("up", "down"))


ThePlot <- ggplot(sfpqkd_de_introns_df[sfpqkd_de_introns_df$rloop_status == "R-loop containing",], 
    aes(x = rci_status, y = log10(rloop_width)))
ThePlot <- ThePlot + geom_violin(aes(fill = direction), colour = "#000000", shape = 21, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + geom_boxplot(aes(fill = direction), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9), alpha = 0.6)
ThePlot <- ThePlot + scale_x_discrete("RCI status")
ThePlot <- ThePlot + scale_y_continuous("log10 total R-loop width in introns")
# ThePlot <- ThePlot + scale_colour_manual("Condition colour", values = colours)
ThePlot <- ThePlot + scale_fill_manual("Condition colour", values = colours)
# ThePlot <- ThePlot + scale_fill_distiller(palette = "RdBu", direction = -1)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_wrap(. ~ rci_status, scales = "free")
ThePlot <- ThePlot + ggtitle("R-loop density in RCIs")
ThePlot <- ThePlot + current_theme

ggsave(width = 7, height = 10, file = file.path(analysis_plot_dir, "03-01-03-plot-total-rloop-width-in-introns.pdf"), ThePlot)

ThePlot <- ggplot(sfpqkd_de_introns_df[sfpqkd_de_introns_df$rloop_status == "R-loop containing",], 
    aes(x = rci_status, y = log10(longest_rloop_width)))
ThePlot <- ThePlot + geom_violin(aes(fill = direction), colour = "#000000", shape = 21, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + geom_boxplot(aes(fill = direction), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9), alpha = 0.6)
ThePlot <- ThePlot + scale_x_discrete("RCI status")
ThePlot <- ThePlot + scale_y_continuous("log10 longest R-loop width in introns")
# ThePlot <- ThePlot + scale_colour_manual("Condition colour", values = colours)
ThePlot <- ThePlot + scale_fill_manual("Condition colour", values = colours)
# ThePlot <- ThePlot + scale_fill_distiller(palette = "RdBu", direction = -1)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_wrap(. ~ rci_status, scales = "free")
ThePlot <- ThePlot + ggtitle("R-loop density in RCIs")
ThePlot <- ThePlot + current_theme

ggsave(width = 7, height = 10, file = file.path(analysis_plot_dir, "03-01-03-plot-longest-rloop-width-in-introns.pdf"), ThePlot)

# median(sfpqkd_de_introns_df[
#     sfpqkd_de_introns_df$rloop_status == "R-loop containing" & 
#     sfpqkd_de_introns_df$rci_status != "RCI" &
#     sfpqkd_de_introns_df$direction == "up","mean_rloop_width"])




ThePlot <- ggplot(sfpqkd_de_introns_df[sfpqkd_de_introns_df$rloop_status == "R-loop containing",], 
    aes(x = direction, y = log10(rloop_density + 1)))
ThePlot <- ThePlot + geom_quasirandom(aes(fill = direction), colour = "#000000", shape = 21, dodge.width = 0.9)
ThePlot <- ThePlot + geom_boxplot(aes(fill = direction), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9), alpha = 0.6)
ThePlot <- ThePlot + scale_x_discrete("Group labels")
ThePlot <- ThePlot + scale_y_continuous("log10 SFPQ iiCLIP signal per kilobase of intron")
# ThePlot <- ThePlot + scale_colour_manual("Condition colour", values = colours)
ThePlot <- ThePlot + scale_fill_manual("Condition colour", values = colours)
# ThePlot <- ThePlot + scale_fill_distiller(palette = "RdBu", direction = -1)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_wrap(. ~ rci_status, scales = "free")
ThePlot <- ThePlot + ggtitle("R-loop density in RCIs")
ThePlot <- ThePlot + current_theme

ggsave(width = 7, height = 10, file = file.path(analysis_plot_dir, "03-01-03-plot-rloop-density-in-introns.pdf"), ThePlot)



median(sfpqkd_de_introns_df$rloop_density[
        sfpqkd_de_introns_df$rloop_status == "R-loop containing" &
        sfpqkd_de_introns_df$rci_status == "RCI" & 
        sfpqkd_de_introns_df$direction == "up"])
# 25.45531

median(sfpqkd_de_introns_df$rloop_density[
        sfpqkd_de_introns_df$rloop_status == "R-loop containing" &
        sfpqkd_de_introns_df$rci_status == "RCI" & 
        sfpqkd_de_introns_df$direction == "down"])
# 11.87728



wilcox.test(sfpqkd_de_introns_df$rloop_density[
        sfpqkd_de_introns_df$rloop_status == "R-loop containing" &
        sfpqkd_de_introns_df$rci_status == "RCI" & 
        sfpqkd_de_introns_df$direction == "up"], 

        sfpqkd_de_introns_df$rloop_density[
        sfpqkd_de_introns_df$rloop_status == "R-loop containing" &
        sfpqkd_de_introns_df$rci_status == "RCI" & 
        sfpqkd_de_introns_df$direction == "down"])

#         Wilcoxon rank sum test with continuity correction

# data:  sfpqkd_de_introns_df$rloop_density[sfpqkd_de_introns_df$rloop_status == "R-loop containing" & sfpqkd_de_introns_df$rci_status == "RCI" & sfpqkd_de_introns_df$direction == "up"] and sfpqkd_de_introns_df$rloop_density[sfpqkd_de_introns_df$rloop_status == "R-loop containing" & sfpqkd_de_introns_df$rci_status == "RCI" & sfpqkd_de_introns_df$direction == "down"]
# W = 113632, p-value = 8.637e-09
# alternative hypothesis: true location shift is not equal to 0
# 

median(sfpqkd_de_introns_df$rloop_density[
        sfpqkd_de_introns_df$rloop_status == "R-loop containing" &
        sfpqkd_de_introns_df$rci_status == "non-RCI" & 
        sfpqkd_de_introns_df$direction == "up"])
# 38.51992

median(sfpqkd_de_introns_df$rloop_density[
        sfpqkd_de_introns_df$rloop_status == "R-loop containing" &
        sfpqkd_de_introns_df$rci_status == "non-RCI" & 
        sfpqkd_de_introns_df$direction == "down"])
# 42.00901




wilcox.test(sfpqkd_de_introns_df$rloop_density[
        sfpqkd_de_introns_df$rloop_status == "R-loop containing" &
        sfpqkd_de_introns_df$rci_status == "non-RCI" & 
        sfpqkd_de_introns_df$direction == "up"], 

        sfpqkd_de_introns_df$rloop_density[
        sfpqkd_de_introns_df$rloop_status == "R-loop containing" &
        sfpqkd_de_introns_df$rci_status == "non-RCI" & 
        sfpqkd_de_introns_df$direction == "down"])

#         Wilcoxon rank sum test with continuity correction

# data:  sfpqkd_de_introns_df$rloop_density[sfpqkd_de_introns_df$rloop_status == "R-loop containing" & sfpqkd_de_introns_df$rci_status == "non-RCI" & sfpqkd_de_introns_df$direction == "up"] and sfpqkd_de_introns_df$rloop_density[sfpqkd_de_introns_df$rloop_status == "R-loop containing" & sfpqkd_de_introns_df$rci_status == "non-RCI" & sfpqkd_de_introns_df$direction == "down"]
# W = 5891618, p-value = 0.9036
# alternative hypothesis: true location shift is not equal to 0

