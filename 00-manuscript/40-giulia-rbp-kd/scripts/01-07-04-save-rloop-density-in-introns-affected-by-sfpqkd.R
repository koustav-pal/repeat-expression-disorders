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
analysis_tmp_dir <- file.path(analysis_dir, "tmp")
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
sfpqkd_de_introns_df <- readRDS(file.path(analysis_save_dir, "04-04-introns-affected-by-SFPQKD-with-SFPQ-binding-signal.rds"))
dripseq_list <- readRDS(file.path("../34-dripseq-analysis/rdata_files/", "02-01-list-of-total-rloop-coverage-matrix-and-library-sizes.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rslurm")
require("stringr")
require("GenomicRanges")
require("DESeq2")
require("ggplot2")
require("parallel")

reduced_rloop_ranges <- dripseq_list$reduced_rloop_ranges
reduced_rloop_ranges_split <- split(width(reduced_rloop_ranges), reduced_rloop_ranges$feature_id)
width_in_rloop <- do.call(c, mclapply(reduced_rloop_ranges_split, sum, mc.cores = 8))
mean_width_of_rloop <- do.call(c, mclapply(reduced_rloop_ranges_split, mean, mc.cores = 8))
longest_rloop_width <- do.call(c, mclapply(reduced_rloop_ranges_split, max, mc.cores = 8))


sfpqkd_de_introns_df$intronic_width <- width(feature_ranges[match(sfpqkd_de_introns_df$hash_id, feature_ranges$hash_id)])
sfpqkd_de_introns_df$rloop_width <- width_in_rloop[feature_ranges$feature_id[match(sfpqkd_de_introns_df$hash_id, feature_ranges$hash_id)]]
sfpqkd_de_introns_df$rloop_width[is.na(sfpqkd_de_introns_df$rloop_width)] <- 0
sfpqkd_de_introns_df$mean_rloop_width <- mean_width_of_rloop[feature_ranges$feature_id[match(sfpqkd_de_introns_df$hash_id, feature_ranges$hash_id)]]
sfpqkd_de_introns_df$mean_rloop_width[is.na(sfpqkd_de_introns_df$mean_rloop_width)] <- 0
sfpqkd_de_introns_df$longest_rloop_width <- longest_rloop_width[feature_ranges$feature_id[match(sfpqkd_de_introns_df$hash_id, feature_ranges$hash_id)]]
sfpqkd_de_introns_df$longest_rloop_width[is.na(sfpqkd_de_introns_df$longest_rloop_width)] <- 0

sfpqkd_de_introns_df$rloop_density <- sfpqkd_de_introns_df$rloop_width / (sfpqkd_de_introns_df$intronic_width / 1e3)

saveRDS(object = sfpqkd_de_introns_df, file = file.path(analysis_save_dir, "04-05-introns-affected-by-SFPQKD-with-rloop-density.rds"))