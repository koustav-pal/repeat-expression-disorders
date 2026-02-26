library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- here("final-figures", "figure-5")
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
analysis_slurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS(here("05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds"))
sig_result_list <- readRDS(here("20-import-significant-events/rdata_files/02-significantly-differentially-expressed-events.rds"))
slurm_job_list <- readRDS(file.path(analysis_save_dir, "05-09-submit-jobs-to-compute-rbp-metaprofile.rds"))
repeat_ranges_list <- readRDS(here("final-figures/figure-1/rdata_files/00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds"))
metaprofile_results_df <- readRDS(here("30-sequence-analysis-of-introns/rdata_files/10-00-rbp-metaprofile-across-all-rcis.rds"))
sfpq_kd_rcis_df <- readRDS(here("40-giulia-rbp-kd/rdata_files/04-06-introns-affected-by-SFPQKD-with-rloop-density-and-rci-penetrance.rds"))
sfpq_kd_rcis <- rownames(sfpq_kd_rcis_df)[sfpq_kd_rcis_df$rci_status == "RCI"]
rloop_containing_sfpq_kd_rcis <- rownames(sfpq_kd_rcis_df)[sfpq_kd_rcis_df$rci_status == "RCI" & 
    sfpq_kd_rcis_df$rloop_status == "R-loop containing"]
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
source("../../src/plot-utils.R")
source("../../src/basepair_coverage_functions.R")
source("../src/colours.R")
require("cowplot")
require("ggplot2")
require("GenomicRanges")
require("HiCBricks")
require("stringr")
require("RColorBrewer")
require("rslurm")
require("dplyr")
require("rslurm")

feature_ranges_with_repeats <- repeat_ranges_list$repeat_containing_ranges
feature_ranges_with_repeats <- feature_ranges_with_repeats[feature_ranges_with_repeats$type == "intron"]
repeat_ranges <- repeat_ranges_list$repeat_ranges
feature_ranges_with_repeats$repeat_id <- names(repeat_ranges)[feature_ranges_with_repeats$which_repeat]
sfpq_kd_rci_ranges <- feature_ranges_with_repeats[feature_ranges_with_repeats$hash_id %in% sfpq_kd_rcis]
sfpq_kd_repeat_ranges <- repeat_ranges[sfpq_kd_rci_ranges$which_repeat]
sfpq_kd_rloop_rci_ranges <- feature_ranges_with_repeats[feature_ranges_with_repeats$hash_id %in% rloop_containing_sfpq_kd_rcis]
sfpq_kd_rloop_rci_repeat_ranges <- repeat_ranges[sfpq_kd_rloop_rci_ranges$which_repeat]


output_prefix = function(){
    return(format(Sys.Date(),"%Y-%m-%d"))
}

metaprofile_results_df$type <- factor(metaprofile_results_df$type, levels = c("upstream", "downstream"))
unique_rbps <- unique(metaprofile_results_df$rbp)


require(parallel)
selected_rci_metaprofile_df <- metaprofile_results_df[metaprofile_results_df$repeat_id %in% names(sfpq_kd_rloop_rci_repeat_ranges),]
# 2 bins is equivalent to 20bp on each side, i.e. 1 bin = 10bp
# That's 40bp around the repeat region
enrichment_around_repeat_df <- do.call(rbind, mclapply(unique_rbps, function(x){
    current_meta_df <- selected_rci_metaprofile_df[selected_rci_metaprofile_df$rbp == x,]
    current_meta_df_split <- split(current_meta_df, current_meta_df$type)
    meta_matrix <- reshape2::acast(data = current_meta_df, 
        formula = repeat_id ~ position, value.var = "values", 
        fun.aggregate = mean)
    selected_window <- c(-2:-1,1:2)
    repeat_classes <- current_meta_df$repeat_class[match(rownames(meta_matrix), current_meta_df$repeat_id)]
    mean_coverage_near_repeats <- rowMeans(meta_matrix[,as.character(selected_window)])
    mean_coverage_in_window <- rowMeans(meta_matrix)
    enrichment_values <- mean_coverage_near_repeats / mean_coverage_in_window
    enrichment_values_split <- split(enrichment_values, repeat_classes)
    # mean_enrichment_matrix <- 
    enrichment_df <- do.call(rbind, lapply(names(enrichment_values_split), function(y){
        current_value <- enrichment_values_split[[y]]
        current_value <- current_value[!is.na(current_value) & !is.infinite(current_value)]
        value_margin <- qt(0.975,df=length(current_value)-1) * sd(current_value)/sqrt(length(current_value))
        data.frame(rbp = x,
            repeat_class = y, 
            lower = mean(current_value) - value_margin,
            enrichment = mean(current_value),
            upper = mean(current_value) + value_margin)
    }))
    return(enrichment_df)
},mc.cores = 32))

plot_list <- lapply(split(enrichment_around_repeat_df, enrichment_around_repeat_df$repeat_class), function(class_specific_enrichment_df){
    class_specific_enrichment_df$rbp <- factor(class_specific_enrichment_df$rbp, levels = class_specific_enrichment_df$rbp[order(class_specific_enrichment_df$enrichment, decreasing = F)])
    # get top 10 and bottom 10
    all_rbps <- levels(class_specific_enrichment_df$rbp)
    selected_rbps <- all_rbps[c((length(all_rbps) - 9) : length(all_rbps))]
    selected_class_specific_enrichment_df <- class_specific_enrichment_df[class_specific_enrichment_df$rbp %in% selected_rbps,]
    selected_class_specific_enrichment_df$log_enrichment <- log2(selected_class_specific_enrichment_df$enrichment)
    selected_class_specific_enrichment_df$log_upper <- log2(selected_class_specific_enrichment_df$upper)
    selected_class_specific_enrichment_df$log_upper[is.na(selected_class_specific_enrichment_df$log_upper)] <- selected_class_specific_enrichment_df$log_enrichment[is.na(selected_class_specific_enrichment_df$log_upper)]
    selected_class_specific_enrichment_df$log_lower <- log2(selected_class_specific_enrichment_df$lower)
    selected_class_specific_enrichment_df$log_lower[is.na(selected_class_specific_enrichment_df$log_lower)] <- selected_class_specific_enrichment_df$log_enrichment[is.na(selected_class_specific_enrichment_df$log_lower)]

    ThePlot <- ggplot(selected_class_specific_enrichment_df, aes(x = rbp, y = log_enrichment))
    ThePlot <- ThePlot + geom_pointrange(aes(ymin = log_lower, ymax = log_upper, colour = repeat_class))
    ThePlot <- ThePlot + geom_hline(yintercept = 0, colour = "#000000",linetype = "dashed")
    ThePlot <- ThePlot + scale_x_discrete("RNA binding proteins")
    ThePlot <- ThePlot + scale_y_continuous("log2 enrichment of motif around repeat")
    ThePlot <- ThePlot + scale_colour_manual("Repeat classes", values = repeat_colour_groups)
    # ThePlot <- ThePlot + scale_colour_manual("Intron groups", values = Colours)
    ThePlot <- ThePlot + coord_flip() + theme_bw()
    ThePlot <- ThePlot + ggtitle("RBP motif enrichment near repeats")
    ThePlot <- ThePlot + theme(
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
    return(ThePlot)
})

require(patchwork)

all_plots <- wrap_plots(plot_list, nrow = 1, ncol = 4)

ggsave(width = 10, height = 5, 
    file = file.path(analysis_plot_dir, paste(output_prefix(), "05-12-plot-positional-preference-of-rbps-sorrounding-repeat-regions.pdf", sep = "-")), all_plots)



ind_repeat_enrichment_df <- do.call(rbind, mclapply(unique_rbps, function(x){
    current_meta_df <- selected_rci_metaprofile_df[selected_rci_metaprofile_df$rbp == x,]
    current_meta_df_split <- split(current_meta_df, current_meta_df$type)
    meta_matrix <- reshape2::acast(data = current_meta_df, 
        formula = repeat_id ~ position, value.var = "values", 
        fun.aggregate = mean)
    selected_window <- c(-2:-1,1:2)
    repeat_classes <- current_meta_df$repeat_class[match(rownames(meta_matrix), current_meta_df$repeat_id)]
    mean_coverage_near_repeats <- rowMeans(meta_matrix[,as.character(selected_window)])
    mean_coverage_in_window <- rowMeans(meta_matrix)
    enrichment_values <- mean_coverage_near_repeats / mean_coverage_in_window
    enrichment_df <- data.frame(repeat_id = names(enrichment_values), 
        repeat_class = repeat_classes,
        enrichment = log2(enrichment_values), 
        rbp = x)
},mc.cores = 32))



drosha_driven_hexanucleotide_repeat_df <- ind_repeat_enrichment_df[
    ind_repeat_enrichment_df$rbp %in% c("DROSHA") & 
    ind_repeat_enrichment_df$repeat_class == "hexanucleotide",]

selected_drosha_repeat_df <- drosha_driven_hexanucleotide_repeat_df[drosha_driven_hexanucleotide_repeat_df$enrichment > log2(1.1),]
selected_drosha_repeat_df$hash_id <- feature_ranges_with_repeats$hash_id[match(selected_drosha_repeat_df$repeat_id, feature_ranges_with_repeats$repeat_id)]

saveRDS(object = selected_drosha_repeat_df, file = file.path(analysis_save_dir, "10-01-hexanucleotide-repeats-with-increased-enrichment-for-drosha.rds"))


selected_ind_repeat_enrichment_df <- do.call(rbind, lapply(split(enrichment_around_repeat_df, enrichment_around_repeat_df$repeat_class), function(class_specific_enrichment_df){
    class_specific_enrichment_df$rbp <- factor(class_specific_enrichment_df$rbp, levels = class_specific_enrichment_df$rbp[order(class_specific_enrichment_df$enrichment, decreasing = F)])
    # get top 10 and bottom 10
    all_rbps <- levels(class_specific_enrichment_df$rbp)
    selected_rbps <- all_rbps[c((length(all_rbps) - 9) : length(all_rbps))]

    selected_ind_repeat_enrichment_df <- ind_repeat_enrichment_df[
        ind_repeat_enrichment_df$repeat_class == unique(class_specific_enrichment_df$repeat_class) &
        ind_repeat_enrichment_df$enrichment > log2(1.1) & !is.na(ind_repeat_enrichment_df$enrichment) &
        ind_repeat_enrichment_df$rbp %in% selected_rbps,]
    selected_ind_repeat_enrichment_df$hash_id <- feature_ranges_with_repeats$hash_id[match(selected_ind_repeat_enrichment_df$repeat_id, feature_ranges_with_repeats$repeat_id)]

    return(selected_ind_repeat_enrichment_df)
}))


selected_ind_repeat_enrichment_df_split <- split(selected_ind_repeat_enrichment_df, paste(selected_ind_repeat_enrichment_df$rbp, selected_ind_repeat_enrichment_df$repeat_class))
saveRDS(object = selected_ind_repeat_enrichment_df, file = file.path(analysis_save_dir, "10-02-all-repeats-with-increased-enrichment-for-top10-proteins.rds"))