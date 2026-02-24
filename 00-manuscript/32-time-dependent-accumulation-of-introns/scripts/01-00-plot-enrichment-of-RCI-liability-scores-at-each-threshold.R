library("here")
message("Project root:", here(), "\n")
project_dir <- here()
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
analysis_dir <- file.path(project_dir, "32-time-dependent-accumulation-of-introns")
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
feature_ranges <- readRDS(here("05-simulated-reads-differential-profiling", "rdata_files", "01-all-merged-exonic-and-intronic-segments-in-genes.rds"))
repeat_object_list <- readRDS(here("01-repeat-characterisation", "rdata_files", "00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds"))
sig_result_list <- readRDS(here("20-import-significant-events", "rdata_files", "02-significantly-differentially-expressed-events.rds"))
nygc_object_list <- readRDS(here("13-nygc-consortium", "rdata_files", "03-02-list_of_counts_stat_samplesheet_from_nygc.rds"))
sampleinfo_df <- read.csv(here("13-nygc-consortium", "input_files", "00_filtered_sample_gsm_list.txt"))
clinical_metadata_df <- read.csv(here("0103-NYGC-consortium", "input_files", "samplesheet.csv"))
patient_df <- readRDS(here("13-nygc-consortium", "rdata_files", "03-05-nygc-patient-sampleinfo-datasheet.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
source(here("..","useful-functions","plot-utils.R"))
source(here("..","useful-functions","basepair_coverage_functions.R"))

require("cowplot")
require("ggplot2")
require("GenomicRanges")
require("HiCBricks")
require("stringr")
require("RColorBrewer")
require("rslurm")
require("dplyr")
require("rslurm")
require("ggbeeswarm")

repeat_ranges <- repeat_object_list$repeat_ranges
feature_ranges_with_repeats <- sig_result_list$repeat_containing_ranges
feature_ranges_with_repeats$repeat_unit_size <- str_length(repeat_ranges$ssr)[feature_ranges_with_repeats$which_repeat]
feature_ranges_with_repeats$span_group <- repeat_ranges$span_group[feature_ranges_with_repeats$which_repeat]
feature_ranges_with_repeats$expansion_size <- feature_ranges_with_repeats$repeat_width / feature_ranges_with_repeats$repeat_unit_size
sig_result_df <- sig_result_list$significant_events

repeat_containing_introns_df <- sig_result_df[sig_result_df$status == "repeat_containing" & 
    sig_result_df$region == "intron" & sig_result_df$contrast %in% c("fals_vs_control", "sals_vs_control",
    "diff_expr_c9orf72_vs_ctrl", "diff_expr_c9orf72_vs_iso", "diff_expr_fus_vs_ctrl",
    "als_vs_control", "als+ftd_vs_control", "diff_expr_sod1_vs_ctrl", "diff_expr_sod1_vs_iso", 
    "diff_expr_tardbp_vs_ctrl", "diff_expr_vcp_vs_ctrl", "diff_expr_d35_mut_vs_wt", 
    "diff_expr_tardbp_vs_wt", "diff_expr_vcp_vs_wt") & 
    sig_result_df$study != "nygc_tissue_non_specific",]


repeat_containing_introns_df$cohort <- "public-data"
repeat_containing_introns_df$cohort[grepl(x = repeat_containing_introns_df$study, pattern = "nygc")] <- "nygc"
repeat_containing_introns_df$cohort[grepl(x = repeat_containing_introns_df$study, pattern = "answerals")] <- "answerals"

repeat_containing_introns_df$tissue_type <- "iPSCs"
repeat_containing_introns_df$tissue_type[grepl(x = repeat_containing_introns_df$study, pattern = "nygc")] <- "postmortem"

counts_matrix <- nygc_object_list$counts
cpm_matrix <- t(t(counts_matrix) / (colSums(counts_matrix) / 1e6))
temp_matrix <- counts_matrix / (width(feature_ranges)[match(rownames(counts_matrix), feature_ranges$hash_id)] / 1e3)
tpm_matrix <- t(t(temp_matrix) / (colSums(temp_matrix) / 1e6))



keep_hash_ids <- feature_ranges_with_repeats$hash_id[match(unique(repeat_containing_introns_df$feature_ids), feature_ranges_with_repeats$feature_id)]
keep_cpm_matrix <- cpm_matrix[keep_hash_ids,colnames(cpm_matrix) %in% patient_df$column_id]
keep_tpm_matrix <- tpm_matrix[keep_hash_ids,colnames(cpm_matrix) %in% patient_df$column_id]

repeat_lengths <- feature_ranges_with_repeats$repeat_width[match(rownames(keep_cpm_matrix), feature_ranges_with_repeats$hash_id)]


coldata_df_split <- split(patient_df, patient_df$tissue)
patient_plot_df <- do.call(rbind, lapply(coldata_df_split, function(current_patient_df){
    df_split <- split(current_patient_df, current_patient_df$subject_id)
    tpm_matrix <- do.call(cbind, lapply(df_split, function(temp_df){
        return_df <- unique(temp_df[,-c(1,2,3,4,5,6,7,8,9,11,19,21)])
        current_columns <- temp_df$column_id
        if(nrow(return_df) > 1){
            if("N/A" %in% return_df$mnd_with_ftd & "Yes" %in% return_df$mnd_with_ftd){
                return_df <- return_df[return_df$mnd_with_ftd == "Yes",]
            }else if("N/A" %in% return_df$mnd_with_ftd & "No" %in% return_df$mnd_with_ftd){
                return_df <- return_df[return_df$mnd_with_ftd == "No",]
            }else if("tardbp" %in% return_df$type & "sporadic" %in% return_df$type){
                return_df <- return_df[return_df$type == "tardbp",]
            }else if("fus" %in% return_df$mutation & "none" %in% return_df$mutation){
                return_df <- return_df[return_df$mutation == "fus",]
            }else{
                stop(message(paste(return_df$subject_id, return_df$tissue)))
            }
        }
        current_tpm_vector <- Reduce("+", lapply(current_columns, function(x){
            keep_tpm_matrix[,x]  
        })) / length(current_columns)
        return(current_tpm_vector)
    }))
    mean_feature_tpms <- rowMeans(tpm_matrix)
    tpm_thresholds <- c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)
    total_df <- do.call(rbind, lapply(seq_len(ncol(tpm_matrix)), function(x){
        current_patient <- tpm_matrix[,x]
        feature_cumsum <- vapply(tpm_thresholds, function(y){
            sum(current_patient[(current_patient) > y])
        },1)
        n_features <- vapply(tpm_thresholds, function(y){
            length(current_patient[(current_patient) > y])
        },1)
        data.frame(subject_id = colnames(tpm_matrix)[x], 
            tissue = unique(current_patient_df$tissue),
            condition = current_patient_df$condition[match(colnames(tpm_matrix)[x], current_patient_df$subject_id)],
            threshold = tpm_thresholds,
            n_features = n_features,
            prop_features_recovered = n_features / n_features[1],
            cumsum = feature_cumsum,
            prop_signal_recovered = feature_cumsum / feature_cumsum[1])
    }))
    # total_df_split <- split(total_df, total_df$condition)
    # mean_cdf_df <- do.call(rbind, lapply(total_df_split, function(temp_df){
    #     mean_signal <- rowMeans(reshape2::acast(data = temp_df, formula = threshold ~ subject_id, value.var = "cumsum"))
    #     standard_dev <- rowSds(reshape2::acast(data = temp_df, formula = threshold ~ subject_id, value.var = "cumsum"))
    #     n <- length(unique(temp_df$subject_id))
    #     margin <- qt(0.975,df=n-1)*standard_dev/sqrt(n)
    #     lower_int <- mean_signal - margin
    #     upper_int <- mean_signal + margin
    #     data.frame(
    #         condition = unique(temp_df$condition),
    #         quantiles = as.numeric(names(mean_signal)), 
    #         lower = lower_int,
    #         mean_cdf = mean_signal,
    #         upper = upper_int)
    # }))
    # mean_cdf_df$tissue <- unique(total_df$tissue)
    return(total_df)
}))
rownames(patient_plot_df) <- NULL
# patient_plot_df$tissue <- droplevels(patient_plot_df$tissue)

patient_plot_df$tissue_group <- "spinal_cord"
patient_plot_df$tissue_group[patient_plot_df$tissue %in% c("Cortex", "Cortex_Motor", "Cerebellum")] <- "brain"
patient_plot_df$tissue <- factor(patient_plot_df$tissue, levels = c("Spinal_Cord_Cervical", "Spinal_Cord_Thoracic", "Spinal_Cord_Lumbar",
    "Cortex", "Cortex_Motor", "Cerebellum"))

patient_plot_df$threshold <- factor(patient_plot_df$threshold, levels = c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1))



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

colours <- c("ctrl" = "#F6F4F3", "als" = "#F87575", "als+ftd" = "#FFA9A3")
fill_colours <- c("ctrl" = "#F6F4F3", 
    "c9orf72" = "#169873", "sporadic" = "#F49FBC",
    "tardbp" = "#FFD3BA", "mutation" = "#9EBD6E",
    "none" = "#805D93")

patient_plot_df$condition <- factor(patient_plot_df$condition, levels = c("ctrl", "als", "als+ftd"))

ThePlot <- ggplot(patient_plot_df[patient_plot_df$condition != "als+ftd",], 
    aes(x = threshold, y = prop_signal_recovered,))
ThePlot <- ThePlot + geom_violin(aes(fill = condition), colour = "#000000", 
    lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_boxplot(aes(fill = condition), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
# ThePlot <- ThePlot + geom_point(aes(fill = condition), shape = 21, colour = "#000000", position = position_jitterdodge(dodge.width = 0.9))
# ThePlot <- ThePlot + geom_smooth(method = "lm", formula = y ~ x, aes(colour = als_group))
# ThePlot <- ThePlot + stat_summary(aes(colour = condition))
ThePlot <- ThePlot + scale_x_discrete("Tissues")
ThePlot <- ThePlot + scale_y_continuous("Proportion of signal recovered")
# ThePlot <- ThePlot + scale_colour_manual("", values = colours)
ThePlot <- ThePlot + scale_fill_manual("Condition colour", values = colours)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_wrap(. ~ tissue, scales = "free")
ThePlot <- ThePlot + ggtitle("Repeat containing intron expression across tissues")
ThePlot <- ThePlot + current_theme
ggsave(width = 13, height = 6, file = file.path(analysis_plot_dir, paste(output_prefix(), "01-00-plot-RCI-TPM-cumulative-sum-at-each-threshold-for-each-tissue.pdf", sep = "-")), ThePlot)



ThePlot <- ggplot(patient_plot_df[patient_plot_df$condition != "als+ftd",], 
    aes(x = threshold, y = prop_features_recovered))
ThePlot <- ThePlot + geom_violin(aes(fill = condition), colour = "#000000", 
    lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_boxplot(aes(fill = condition), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
# ThePlot <- ThePlot + geom_point(aes(fill = condition), shape = 21, colour = "#000000", position = position_jitterdodge(dodge.width = 0.9))
# ThePlot <- ThePlot + geom_smooth(method = "lm", formula = y ~ x, aes(colour = als_group))
# ThePlot <- ThePlot + stat_summary(aes(colour = condition))
ThePlot <- ThePlot + scale_x_discrete("Tissues")
ThePlot <- ThePlot + scale_y_continuous("Number of features recovered")
# ThePlot <- ThePlot + scale_colour_manual("", values = colours)
ThePlot <- ThePlot + scale_fill_manual("Condition colour", values = colours)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_wrap(. ~ tissue, scales = "free")
ThePlot <- ThePlot + ggtitle("Repeat containing intron expression across tissues")
ThePlot <- ThePlot + current_theme
ggsave(width = 13, height = 6, file = file.path(analysis_plot_dir, paste(output_prefix(), "01-00-plot-RCI-TPM-cumulative-sum-Ns-at-each-threshold-for-each-tissue.pdf", sep = "-")), ThePlot)

expected_df <- do.call(rbind, lapply(split(patient_plot_df, paste(patient_plot_df$subject_id, patient_plot_df$tissue, patient_plot_df$condition)), function(temp_df){
    data.frame(threshold = 0,
        subject_id = unique(temp_df$subject_id),
        tissue = unique(temp_df$tissue),
        condition = unique(temp_df$condition),
        expected_signal = median(temp_df$cumsum[temp_df$threshold == 0]) / median(temp_df$n_features[temp_df$threshold == 0]))
}))
rownames(expected_df) <- NULL



enrichment_df <- do.call(rbind, lapply(split(patient_plot_df, paste(patient_plot_df$threshold, patient_plot_df$subject_id, patient_plot_df$tissue, patient_plot_df$condition)), function(temp_df){
    current_expected <- expected_df[expected_df$subject_id == unique(temp_df$subject_id) &
        expected_df$tissue == unique(temp_df$tissue) & 
        expected_df$condition == unique(temp_df$condition), "expected_signal"]
    return_df <- data.frame(threshold = unique(temp_df$threshold),
        subject_id = unique(temp_df$subject_id),
        tissue = unique(temp_df$tissue),
        condition = unique(temp_df$condition),
        observed = median(temp_df$cumsum) / median(temp_df$n_features))
    return_df$enrichment <- return_df$observed / current_expected
    return(return_df)
}))

enrichment_gain_df <- do.call(rbind, lapply(split(enrichment_df, paste(enrichment_df$subject_id, enrichment_df$tissue, enrichment_df$condition)), function(temp_df){
    reordered_df <- temp_df[order(as.numeric(temp_df$threshold), decreasing = F),]
    gain_df <- reordered_df[-1,]
    gain_df$gain_in_enrichment <- diff(reordered_df$enrichment) / reordered_df$enrichment[-nrow(reordered_df)]
    return(gain_df)
}))
rownames(enrichment_gain_df) <- NULL

ThePlot <- ggplot(enrichment_gain_df[enrichment_gain_df$condition != "als+ftd",], 
    aes(x = threshold, y = gain_in_enrichment))
ThePlot <- ThePlot + geom_violin(aes(fill = condition), colour = "#000000", 
    lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_boxplot(aes(fill = condition), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + scale_x_discrete("RCI TPM thresholds")
ThePlot <- ThePlot + scale_y_continuous("Percentage gain in signal resolution")
# ThePlot <- ThePlot + scale_colour_manual("", values = colours)
ThePlot <- ThePlot + scale_fill_manual("Condition colour", values = colours)
# ThePlot <- ThePlot + scale_fill_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_wrap(. ~ tissue, scales = "free")
ThePlot <- ThePlot + ggtitle("Repeat containing intron expression enrichment at each threshold across all tissues")
ThePlot <- ThePlot + current_theme
ggsave(width = 13, height = 6, file = file.path(analysis_plot_dir, paste(output_prefix(), "01-00-plot-RCI-TPM-enrichment-at-each-threshold-for-each-tissue.pdf", sep = "-")), ThePlot)



median_df <- do.call(rbind, lapply(split(patient_plot_df, paste(patient_plot_df$threshold, patient_plot_df$tissue, patient_plot_df$condition)), function(temp_df){
    current_expected <- expected_df[expected_df$tissue == unique(temp_df$tissue) & 
        expected_df$condition == unique(temp_df$condition),"expected_signal"]
    return_df <- data.frame(threshold = unique(temp_df$threshold),
        tissue = unique(temp_df$tissue),
        condition = unique(temp_df$condition),
        cumsum = median(temp_df$cumsum),
        n_features = median(temp_df$n_features),
        prop_signal_recovered = median(temp_df$prop_signal_recovered),
        prop_features_recovered = median(temp_df$prop_features_recovered),
        feature_to_signal_ratio = median(temp_df$prop_signal_recovered / temp_df$prop_features_recovered),
        observed = median(temp_df$cumsum) / median(temp_df$n_features))
    return_df$enrichment <- return_df$observed / current_expected
    return(return_df)
}))

lapply(split(median_df, median_df$tissue), function(temp_df){
    temp_df[which.max(temp_df$feature_to_signal_ratio),]
})
