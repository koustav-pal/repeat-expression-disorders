# ==========================================================================
# Intial script setup
# ==========================================================================
data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "13-01-answerals")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
isPairedEnd <- T
strand_specificity <- 2
sampleinfo_df <- read.csv("/nemo/project/proj-luscombn-patani/working/answerals/sample-details/samplesheet.csv")
alignment_dir <- "/nemo/project/proj-luscombn-patani/working/"
samples <- c("answerals" = "answerals")
object_list <- readRDS(file.path(analysis_save_dir, "01-02-list_of_counts_stat_samplesheet_from_answerals.rds"))
# ==========================================================================
# Imports
# ==========================================================================
require("GenomicRanges")
require("rslurm")
require("Rsubread")
require("stringr")
require("HiCBricks")
require("DESeq2")
# ==========================================================================
# Analysis
# ==========================================================================

counts_matrix <- object_list$counts

sample_df <- object_list$sample_df
# sample_df <- sample_df[sample_df$sample_name %in% sampleinfo_df$Experiment,]

sampleinfo_df$sample_name <- str_replace_all(sampleinfo_df$participant_id, pattern = "-", replacement = "_")

sample_df$condition <- sampleinfo_df$condition[match(sample_df$sample_name, sampleinfo_df$sample_name)]
sample_df$mutation <- sampleinfo_df$mutation[match(sample_df$sample_name, sampleinfo_df$sample_name)]
sample_df$mutant <- sampleinfo_df$mutant[match(sample_df$sample_name, sampleinfo_df$sample_name)]
sample_df$mnd_group <- sampleinfo_df$mnd_group[match(sample_df$sample_name, sampleinfo_df$sample_name)]


sample_df$condition_group <- NA
sample_df$condition_group[sample_df$mutant == "yes" & sample_df$mnd_group == "als"] <- "fals"
sample_df$condition_group[sample_df$mutant == "no" & sample_df$mnd_group == "als"] <- "sals"
sample_df$condition_group[sample_df$mnd_group == "healthy"] <- "ctrl"
sample_df$condition_group[sample_df$mnd_group == "other_mnd"] <- "other_mnd"

sample_df$onset_site_detailed <- sampleinfo_df$onset_site_detailed[match(sample_df$sample_name, sampleinfo_df$sample_name)]
sample_df$onset_site_detailed[sample_df$condition_group == "ctrl"] <- "none"
sample_df$onset_site_detailed[
    sample_df$condition_group != "ctrl" & 
    sample_df$onset_site_detailed == ""] <- "unknown"


sample_df$progression <- sampleinfo_df$progression[match(sample_df$sample_name, sampleinfo_df$sample_name)]
sample_df$progression[sample_df$condition_group == "ctrl"] <- "none"
sample_df$progression[sample_df$condition_group != "ctrl" & sample_df$progression == ""] <- "unknown"

sample_df$sex <- sampleinfo_df$sex[match(sample_df$sample_name, sampleinfo_df$sample_name)]
sample_df$age_onset <- sampleinfo_df$age_onset[match(sample_df$sample_name, sampleinfo_df$sample_name)]


selected_counts_matrix <- round(counts_matrix)[,sample_df$column_id]
selected_counts_matrix <- selected_counts_matrix[rowMeans(selected_counts_matrix) > 1,]
# selected_counts_matrix_with_pseudo <- selected_counts_matrix + 1

sample_df <- sample_df[!is.na(sample_df$condition_group),]

sample_df$condition_group <- factor(sample_df$condition_group, levels = c("ctrl", "fals", "sals", "other_mnd"))
sample_df$sex <- factor(sample_df$sex, levels = c("Female", "Male"))
sample_df$progression <- factor(sample_df$progression, levels = c("none", "slow", "moderate", "fast", "unknown"))
sample_df$onset_site_detailed <- factor(sample_df$onset_site_detailed, levels = c("none", "bulbar", "limb", "axial", "mixed", "unknown"))

simple_model <- ~ condition_group
simple_model_plus_sex <- ~ sex + condition_group 
# simple_model_plus_onset_site_plus_progression <- ~ progression + onset_site_detailed + condition_group 
# simple_model_plus_onset_site_plus_progression_plus_sex <- ~ sex + progression + onset_site_detailed + condition_group 


deseq_object <- DESeqDataSetFromMatrix(countData = selected_counts_matrix[,sample_df$column_id], 
    colData = sample_df, 
    design = simple_model)
simple_differential_object <- DESeq(deseq_object)
saveRDS(object = simple_differential_object, file = file.path(analysis_save_dir, "02-01-01-save-deseq2-simple-differential-object.rds"))

design(deseq_object) <- simple_model_plus_sex
sex_differential_object <- DESeq(deseq_object)
saveRDS(object = sex_differential_object, file = file.path(analysis_save_dir, "02-01-02-save-deseq2-sex-differential-object.rds"))

# select only ALS specific 
current_sample_df <- sample_df[sample_df$condition_group != "ctrl" & sample_df$progression != "unknown",]
disease_model_plus_progression <- ~ sex + condition_group + progression + progression:condition_group
reduced_model <- ~ sex + condition_group + progression 
deseq_object <- DESeqDataSetFromMatrix(countData = selected_counts_matrix[,current_sample_df$column_id], 
    colData = current_sample_df, 
    design = disease_model_plus_progression)
progression_differential_object <- DESeq(deseq_object, test = "LRT", reduced = reduced_model)
saveRDS(object = progression_differential_object, file = file.path(analysis_save_dir, "02-01-03-save-deseq2-disease-progression-differential-object.rds"))

# fetch fALS specific onset  
current_sample_df <- sample_df[!(sample_df$condition_group %in% c("ctrl", "other_mnd")) & !(sample_df$onset_site_detailed %in% c("unknown", "axial")),]
disease_model_plus_location <- ~ sex + condition_group + onset_site_detailed + condition_group:onset_site_detailed
reduced_model <- ~ sex + condition_group + onset_site_detailed 
deseq_object <- DESeqDataSetFromMatrix(countData = selected_counts_matrix[,current_sample_df$column_id], 
    colData = current_sample_df, 
    design = disease_model_plus_location)
location_differential_object <- DESeq(deseq_object, test = "LRT", reduced = reduced_model)
saveRDS(object = progression_differential_object, file = file.path(analysis_save_dir, "02-01-04-save-deseq2-location-specific-differences-sals-vs-fals-differential-object.rds"))