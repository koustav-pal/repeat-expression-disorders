data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "20-import-significant-events")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
repeat_object_list <- readRDS(file.path("../01-repeat-characterisation/rdata_files", "00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds"))
answerals_als_vs_control_results_df <- readRDS(file.path(analysis_save_dir, "01-03-import-answerals-datasets-differential-expression.rds"))
nygc_als_vs_control_results_df <-  readRDS(file.path(analysis_save_dir, "01-03-import-nygc-datasets-differential-expression.rds"))
public_als_vs_control_results_df <- readRDS(file.path(analysis_save_dir, "01-02-import-public-datasets-differential-expression.rds"))
answerals_list <- readRDS("../13-01-answerals/rdata_files/01-02-list_of_counts_stat_samplesheet_from_answerals.rds")
public_data_list <- readRDS("../13-03-public-data/rdata_files/02-list_of_counts_stat_samplesheet_from_public_data.rds")
annotation_gtf <- "/nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
# ==========================================================================
# Packages and Data
# ==========================================================================
# source("../useful-functions/plot-utils.R")
source("../../useful-functions/plot-utils.R")
source("../../useful-functions/basepair_coverage_functions.R")
source("scripts/src/functions.R")
require("ggplot2")
require("GenomicRanges")
require("stringr")
require("RColorBrewer")
require("dplyr")
require("survival")
require("matrixStats")
require("rtracklayer")

annotation_ranges <- import.gff2(annotation_gtf)

repeat_ranges <- repeat_object_list$repeat_ranges
ol_object <- findOverlaps(repeat_ranges, feature_ranges)
feature_ranges_with_repeats <- feature_ranges[subjectHits(ol_object)]
feature_ranges_with_repeats$gene_name <- annotation_ranges$gene_name[match(feature_ranges_with_repeats$gene_id, annotation_ranges$gene_id)]
feature_ranges_with_repeats$which_repeat <- queryHits(ol_object)
feature_ranges_with_repeats$repeat_width <- width(repeat_ranges[queryHits(ol_object)])

all_als_vs_control_results_df <- rbind(
    public_als_vs_control_results_df[public_als_vs_control_results_df$pval < 0.05 & !is.na(public_als_vs_control_results_df$pval),], 
    nygc_als_vs_control_results_df[nygc_als_vs_control_results_df$pval < 0.05 & !is.na(nygc_als_vs_control_results_df$pval),],
    answerals_als_vs_control_results_df[answerals_als_vs_control_results_df$pval < 0.05 & !is.na(answerals_als_vs_control_results_df$pval),])

all_als_vs_control_results_df <- all_als_vs_control_results_df[all_als_vs_control_results_df$contrasts %in% 
    c("fals_vs_control", "sals_vs_control", "diff_expr_c9orf72_vs_ctrl", 
    "diff_expr_c9orf72_vs_iso", "diff_expr_fus_vs_ctrl", "als_vs_control", 
    "als+ftd_vs_control", "diff_expr_sod1_vs_ctrl", "diff_expr_sod1_vs_iso", 
    "diff_expr_tardbp_vs_ctrl", "diff_expr_vcp_vs_ctrl", "diff_expr_d35_mut_vs_wt", 
    "diff_expr_tardbp_vs_wt", "diff_expr_vcp_vs_wt") & all_als_vs_control_results_df$study != "nygc_tissue_non_specific",]

all_genotypic_backgrounds <- c("diff_expr_c9orf72_vs_ctrl" = "c9orf72", "diff_expr_fus_vs_ctrl" = "fus",
    "diff_expr_c9orf72_vs_iso" = "c9orf72", "diff_expr_sod1_vs_ctrl" = "sod1", "diff_expr_sod1_vs_iso" = "sod1",
    "diff_expr_vcp_vs_ctrl" = "vcp", "diff_expr_d35_mut_vs_wt" = "vcp", "diff_expr_tardbp_vs_wt" = "tardbp",
    "diff_expr_tardbp_vs_ctrl" = "tardbp", "diff_expr_vcp_vs_wt" = "vcp", "fals_vs_control" = "familial",
    "sals_vs_control" = "sporadic", "als_vs_control" = "als", "als+ftd_vs_control" = "als+ftd")

all_als_vs_control_results_df$genotype <- all_genotypic_backgrounds[all_als_vs_control_results_df$contrast]
all_als_vs_control_results_df$xaxis <- paste(all_als_vs_control_results_df$study, all_als_vs_control_results_df$contrast, sep = "_")
all_als_vs_control_results_df$direction <- ifelse(all_als_vs_control_results_df$logfc > 0, "up", "down")
all_als_vs_control_results_df$cohort <- "public_data"
all_als_vs_control_results_df$cohort[grepl(x = all_als_vs_control_results_df$study, pattern = "answerals")] <- "AnswerALS"
all_als_vs_control_results_df$cohort[grepl(x = all_als_vs_control_results_df$study, pattern = "nygc")] <- "NYGC"
all_als_vs_control_results_df$dir_name <- paste(all_als_vs_control_results_df$feature_id, all_als_vs_control_results_df$direction, sep = ":")
all_als_vs_control_results_df$gene_id <- str_split(all_als_vs_control_results_df$feature_ids, pattern = ":", simplify = T)[,1]
all_als_vs_control_results_df$gene_name <- annotation_ranges$gene_name[match(all_als_vs_control_results_df$gene_id, annotation_ranges$gene_id)]
all_als_vs_control_results_df$status <- ifelse(all_als_vs_control_results_df$feature_id %in% feature_ranges_with_repeats$feature_id, "repeat_containing", "non_repeat_containing")
all_als_vs_control_results_df$region <- str_split(all_als_vs_control_results_df$feature_ids, pattern = ":", simplify = T)[,6]
all_als_vs_control_results_df$tissue_type <- "iPSCs"

tissue_split <- str_split(unique(all_als_vs_control_results_df$study[grepl(x = all_als_vs_control_results_df$study, pattern = "nygc")]), pattern = "_")
names(tissue_split) <- unique(all_als_vs_control_results_df$study[grepl(x = all_als_vs_control_results_df$study, pattern = "nygc")])
tissue_types <- vapply(tissue_split, function(x){
    paste(x[x %in% c("Spinal", "Cord", "Lumbar", "Cervical", "Cerebellum", "Cortex", "Motor", "Thoracic")], collapse = " ")
},"1")
all_als_vs_control_results_df$tissue_type <- "iPSCs"
all_als_vs_control_results_df$tissue_type[grepl(x = all_als_vs_control_results_df$study, pattern = "nygc")] <- tissue_types[all_als_vs_control_results_df$study[grepl(x = all_als_vs_control_results_df$study, pattern = "nygc")]]


saveRDS(object = list(significant_events = all_als_vs_control_results_df, 
    repeat_containing_ranges = feature_ranges_with_repeats),
    file = file.path(analysis_save_dir, "02-significantly-differentially-expressed-events.rds"))