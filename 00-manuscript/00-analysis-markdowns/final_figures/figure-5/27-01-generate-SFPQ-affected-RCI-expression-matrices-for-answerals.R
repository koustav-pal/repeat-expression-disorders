library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- here("final-figures", "figure-5")
save_dir <- file.path(analysis_dir, "rdata_files")
input_dir <- file.path(analysis_dir,"input_files")
plot_dir <- file.path(analysis_dir, "plots")
output_dir <- file.path(analysis_dir, "output_files")
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
# source("../../useful-functions/plot-utils.R")
# source("../../useful-functions/basepair_coverage_functions.R")
# source("../src/colours.R")
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

# ==========================================================================
# Load data
# ==========================================================================


#
# Load repeats and introns containing repeats
# --------------------------------------------------------------------------
repeat_ranges_list <- readRDS("final-figures/figure-1/rdata_files/00-repeat-ranges-introns-exons-with-mapping-to-exons-and-introns.rds")
repeat_ranges <- unique(repeat_ranges_list$repeat_ranges)
feature_ranges_with_repeats <- repeat_ranges_list$repeat_containing_ranges

#
# Load total feature ranges
# --------------------------------------------------------------------------

feature_ranges <- readRDS("05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")

#
# Load significant events and repeat containing introns
# --------------------------------------------------------------------------

sig_result_list <- readRDS("20-import-significant-events/rdata_files/02-significantly-differentially-expressed-events.rds")
sig_result_df <- sig_result_list$significant_events

repeat_containing_introns_df <- sig_result_df[sig_result_df$status == "repeat_containing" & 
    sig_result_df$region == "intron",]
repeat_containing_introns_df$cohort <- "public-data"
repeat_containing_introns_df$cohort[grepl(x = repeat_containing_introns_df$study, pattern = "nygc")] <- "nygc"
repeat_containing_introns_df$cohort[grepl(x = repeat_containing_introns_df$study, pattern = "answerals")] <- "answerals"

repeat_containing_introns_df$tissue_type <- "iPSCs"
repeat_containing_introns_df$tissue_type[grepl(x = repeat_containing_introns_df$study, pattern = "nygc")] <- "postmortem"

#
# Load Grima et al data
# --------------------------------------------------------------------------

answerals_object_list <- readRDS("39-answerals-complete/rdata_files/01-02-list_of_counts_stat_samplesheet_from_answerals.rds")
patient_df <- answerals_object_list$sample_df

#
# Load metadata
# --------------------------------------------------------------------------

sampleinfo_df <- read.csv("39-answerals-complete/input_files/metadata/aals_dataportal_datatable.csv")
cell_line_markers_csv <- read.csv("13-01-answerals/input_files/Baxi-et-al-IF-raw-data-file-for-Koustav-Pal-confirmed-data.csv")
participant_info_df <- read.csv("39-answerals-complete/input_files/metadata/aals_participants.csv")
external_controls_df <- read.csv("39-answerals-complete/input_files/metadata/external_controls.csv") 
family_history_log <- read.csv("39-answerals-complete/input_files/metadata/clinical/Family_History_Log.csv")
mutation_df <- read.csv("39-answerals-complete/input_files/metadata/clinical/ALS_Gene_Mutations.csv")
mortality_df <- read.csv("39-answerals-complete/input_files/metadata/clinical/Mortality.csv")
ANSASFD_df <- read.csv("39-answerals-complete/input_files/metadata/clinical/ANSASFD.csv")

patient_df$Participant_ID <- vapply(str_split(patient_df$sample_name, pattern = "-"), function(x){
    paste(x[c(1,2)], collapse = "-")
},"1")
patient_df$subject_id <- str_split(string = patient_df$Participant_ID, pattern = "-", simplify = T)[,2]

sample_name_matches <- match(patient_df$Participant_ID, sampleinfo_df$Participant_ID)

patient_df$di_m_ns_line_name <- sampleinfo_df$diMNs_line[sample_name_matches]



#
# Add markers
# --------------------------------------------------------------------------

colnames(cell_line_markers_csv) <- c("diMNs-line", "GUID", "condition", "SMI32", "ISL1", "NKX6.1", "TUJ1", "S100b")
cell_line_markers_csv <- cell_line_markers_csv[seq_len(nrow(cell_line_markers_csv)) < 212,]
non_na_marker_df <- cell_line_markers_csv[!vapply(seq_len(nrow(cell_line_markers_csv)), function(x){
    any(vapply(c("SMI32", "ISL1", "NKX6.1", "TUJ1", "S100b"), function(y){
            cell_line_markers_csv[x,y] == "unavailable"
        },T))
},T),]

non_na_marker_df$SMI32 <- as.numeric(non_na_marker_df$SMI32)
non_na_marker_df$ISL1 <- as.numeric(non_na_marker_df$ISL1)
non_na_marker_df$NKX6.1 <- as.numeric(non_na_marker_df$NKX6.1)
non_na_marker_df$TUJ1 <- as.numeric(non_na_marker_df$TUJ1)
non_na_marker_df$S100b <- as.numeric(non_na_marker_df$S100b)
dimn_matches <- match(patient_df$di_m_ns_line_name, non_na_marker_df[,"diMNs-line"])

#
# Add metadata
# --------------------------------------------------------------------------

patient_df$marker_data_available <- !is.na(dimn_matches)
patient_df[,"SMI32"] <- sampleinfo_df[sample_name_matches,"PCT_SMI32"]
patient_df[,"ISL1"] <- sampleinfo_df[sample_name_matches,"PCT_ISL1"]
patient_df[,"NKX6.1"] <- sampleinfo_df[sample_name_matches,"PCT_NKX61"]
patient_df[,"TUJ1"] <- sampleinfo_df[sample_name_matches,"PCT_TUJ1"]
patient_df[,"S100b"] <- sampleinfo_df[sample_name_matches,"PCT_S100b"]
patient_df$nygc_cgnd <- sampleinfo_df$NYGC_CGND_ID[sample_name_matches]
patient_df$sex <-  sampleinfo_df$SEX[sample_name_matches]
patient_df$progression <- sampleinfo_df$ALSFRS_R_PROGRESSION_SLOPE[sample_name_matches]
patient_df$age_onset <- sampleinfo_df$AGE_AT_SYMPTOM_ONSET[sample_name_matches]
patient_df$age_at_collection <- sampleinfo_df$Age_at_First_PBMC_Collection[sample_name_matches]
patient_df$age_at_collection[is.na(patient_df$age_at_collection)] <- external_controls_df[match(patient_df$Participant_ID[is.na(patient_df$age_at_collection)], external_controls_df$Participant_ID),"Age_Sample_Collection_Cedars"]
patient_df$Number_of_Visits <- sampleinfo_df$Number_of_Visits[sample_name_matches]
patient_df$family_history <- family_history_log$famrel[match(patient_df$Participant_ID, family_history_log$Participant_ID)]
patient_df$mnd_group <- participant_info_df$Cohort[match(patient_df$Participant_ID, participant_info_df$Participant_ID)]
patient_df$mnd_group[grepl(x = patient_df$Participant_ID, "CTRL")] <- "Healthy Control"
patient_df$mutation <- vapply(patient_df$Participant_ID, function(x){
    mutations <- c("ang", "c9orf72", "fus", "mutot", "progran", "setx", "sod1", "tau", "tdp43", "vapb", "vcp")
    current_mutation_df <- mutation_df[mutation_df$Participant_ID == x,mutations]
    non_na <- colnames(current_mutation_df)[!is.na(current_mutation_df) & current_mutation_df != 2]
    if(length(non_na) > 0){
        return(paste(non_na, collapse = "|"))
    }
    return("NA")
},"1")
patient_df$tested_for_c9orf72_repeat_expansion <- sampleinfo_df[match(patient_df$Participant_ID, sampleinfo_df$Participant_ID),"EH_C9orf72"] == "Yes"
patient_df$c9orf72_repeat_length <- sampleinfo_df[match(patient_df$Participant_ID, sampleinfo_df$Participant_ID),"C9orf72_repeat_length"]

patient_df$tested_for_atxn2_repeat_expansion <- sampleinfo_df[match(patient_df$Participant_ID, sampleinfo_df$Participant_ID),"EH_ATXN2"] == "Yes"
patient_df$atxn2_repeat_length <- sampleinfo_df[match(patient_df$Participant_ID, sampleinfo_df$Participant_ID),"ATXN2_repeat_length"]

patient_df$expansion_group <- "Not tested"
patient_df$expansion_group[(patient_df$tested_for_c9orf72_repeat_expansion | patient_df$tested_for_atxn2_repeat_expansion) & 
    (patient_df$c9orf72_repeat_length >= 24 | patient_df$atxn2_repeat_length >= 24) & 
    !is.na(patient_df$c9orf72_repeat_length) & !is.na(patient_df$atxn2_repeat_length)] <- "Tested and contains a C9ORF72 or ATXN2 repeat length than 24"
patient_df$expansion_group[(patient_df$tested_for_c9orf72_repeat_expansion | patient_df$tested_for_atxn2_repeat_expansion) & 
    (patient_df$c9orf72_repeat_length < 24 & patient_df$atxn2_repeat_length < 24) & 
    !is.na(patient_df$c9orf72_repeat_length) & !is.na(patient_df$atxn2_repeat_length)] <- "Tested and lacks a C9ORF72 or ATXN2 pathological repeat"


patient_df$disease_duration <- mortality_df$dieddt[match(patient_df$subject_id, mortality_df$SubjectUID)]
patient_df$subject_died <- patient_df$subject_id %in% mortality_df$SubjectUID
patient_df$last_known_alive <- ANSASFD_df$lnadt[match(patient_df$subject_id, ANSASFD_df$SubjectUID)]
patient_df$disease_duration[is.na(patient_df$disease_duration) & !is.na(patient_df$last_known_alive)] <- patient_df$last_known_alive[is.na(patient_df$disease_duration) & !is.na(patient_df$last_known_alive)]



#
# Load Nfl data
# --------------------------------------------------------------------------
nfl_df <- read.csv("39-answerals-complete/input_files/metadata/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv")
nfl_values <- vapply(split(nfl_df, nfl_df$SubjID), function(x){
    non_nas <- x$Result.Value[!is.na(x$Result.Value)]
    ret_val <- as.numeric(non_nas[which.min(non_nas)])
    if(length(ret_val) == 0){
        return(NA)
    }
    return(ret_val)
},1)

patient_df$nfl <- nfl_values[patient_df$subject_id]


#
# Load SFPQ KD penetrant RCIs
# --------------------------------------------------------------------------

sfpq_kd_rcis_df <- readRDS(file.path("40-giulia-rbp-kd/rdata_files/04-06-introns-affected-by-SFPQKD-with-rloop-density-and-rci-penetrance.rds"))
sfpq_kd_rcis <- rownames(sfpq_kd_rcis_df[sfpq_kd_rcis_df$rci_status == "RCI",])
rloop_containing_sfpq_kd_rcis <- rownames(sfpq_kd_rcis_df)[sfpq_kd_rcis_df$rci_status == "RCI" & 
    sfpq_kd_rcis_df$rloop_status == "R-loop containing"]

# ==========================================================================
# Functions
# ==========================================================================

prepare_progression_rsf_df <- function(df){
    df$time <- df$disease_duration
    df$status <- 1
    df$status[df$Deceased != 1] <- 0
    return(df)
}

# ==========================================================================
# Generate TPM matrix
# ==========================================================================

counts_matrix <- answerals_object_list$counts
temp_matrix <- counts_matrix / (width(feature_ranges)[match(rownames(counts_matrix), feature_ranges$hash_id)] / 1e3)
tpm_matrix <- t(t(temp_matrix) / (colSums(temp_matrix) / 1e6))

keep_hash_ids <- feature_ranges_with_repeats$hash_id[match(unique(sfpq_kd_rcis), feature_ranges_with_repeats$hash_id)]
keep_tpm_matrix <- tpm_matrix[keep_hash_ids,colnames(tpm_matrix) %in% patient_df$column_id]

# ==========================================================================
# Construct DESeq2 object
# ==========================================================================
require(DESeq2)

keep_patient_df <- patient_df[patient_df$mnd_group %in% c("Healthy Control", "ALS"),]
keep_patient_df$mnd_group[keep_patient_df$mnd_group == "Healthy Control"] <- "Control"
keep_patient_df$mnd_group <- factor(keep_patient_df$mnd_group, levels = c("Control", "ALS"))

deseq_object <- DESeqDataSetFromMatrix(countData = round(counts_matrix[keep_hash_ids,])[,keep_patient_df$column_id],
    colData = keep_patient_df,
    design = ~ mnd_group)

diff_object <- DESeq(deseq_object)
normalized_counts <- counts(diff_object, normalized=TRUE)

# ==========================================================================
# Generate Penetrance matrix
# ==========================================================================

patient_penetrance_matrix <- do.call(rbind, lapply(seq_len(nrow(keep_patient_df)), function(x){
    temp_df <- keep_patient_df[x,]
    current_columns <- temp_df$column_id
    current_tpm_vector <- Reduce("+", lapply(current_columns, function(x){
        normalized_counts[,x]
    })) / length(current_columns)
    current_threshold <- 0.1
    a_vector <- rep(0, length(unique(sfpq_kd_rcis)))
    names(a_vector) <- unique(sfpq_kd_rcis)
    a_vector[names(current_tpm_vector[current_tpm_vector > current_threshold])] <- 1
    return(a_vector)
}))

# ==========================================================================
# Save data
# ==========================================================================

sfpq_affected_rcis_matrix <- normalized_counts[sfpq_kd_rcis,]
keep_patient_df$tissue_group <- "iPSCs"
keep_patient_df$penetrance_row <- match(keep_patient_df$column_id, colnames(sfpq_affected_rcis_matrix))

saveRDS(object = list(patient_df = keep_patient_df, 
    penetrance_matrix = patient_penetrance_matrix,
    sfpq_affected_rcis_matrix = sfpq_affected_rcis_matrix,
    sfpq_affected_rcis = sfpq_kd_rcis,
    sfpq_affected_rloop_containing_rcis = rloop_containing_sfpq_kd_rcis),
    file = file.path(save_dir, "07-01-save-answerals-iPSCs-sfpq-affected-rci-expression-matrix.rds"))
