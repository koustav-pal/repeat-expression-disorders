# ==========================================================================
# Intial script setup
# ==========================================================================
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
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
isPairedEnd <- T
strand_specificity <- 2
samplesheet_df <- read.csv("/nemo/lab/patanir/home/users/ziffo/patani-collab/motor-neuron-rbp-knockdown/sample-details/samplesheet.csv")
alignment_dir <- "/nemo/lab/patanir/home/users/ziffo/patani-collab/motor-neuron-rbp-knockdown/nfcore/star_salmon/"
samples <- c("answerals" = "answerals")
# ==========================================================================
# Imports
# ==========================================================================
require("GenomicRanges")
require("rslurm")
require("Rsubread")
require("stringr")
require("HiCBricks")
# ==========================================================================
# Analysis
# ==========================================================================

files_in_dir <- list.files(alignment_dir)


# isPairedEnd <- vapply(samples, function(x){
# 	current_sample_samplesheet <- file.path(alignment_dir, x, "sample-details", "samplesheet.csv")
# 	sample_df <- read.csv(current_sample_samplesheet)
# 	all(sample_df$library_layout == "PAIRED")
# },T)

# strand_specificity <- vapply(samples, function(x){
# 	current_sample_samplesheet <- file.path(alignment_dir, x, "sample-details", "samplesheet.csv")
# 	sample_df <- read.csv(current_sample_samplesheet)
# 	strandedness <- unique(sample_df$strandedness)
# 	ifelse(strandedness == "reverse", 2, ifelse(strandedness == "forward", 1, 0))
# },1)
samplesheet_df$present_in_dir <- vapply(samplesheet_df$sample, function(x){
	any(grepl(x = list.files(path = alignment_dir, pattern = ".bam"), x))
},T)
keep_samplesheet_df <- samplesheet_df[samplesheet_df$present_in_dir,]
keep_samplesheet_df <- unique(keep_samplesheet_df[,c("sample", "strandedness", "fraction", "condition", "cellline", "sample_name", "sample_treatment", "group")])

all_files_list <- do.call(c,lapply(keep_samplesheet_df$sample, function(x){
	current_samples <- list.files(alignment_dir, pattern = x, recursive = T)
	return(current_samples[grepl(x = current_samples, pattern = "markdup.sorted.bam") & 
		!grepl(x = current_samples, pattern = "bai") &
		!grepl(x = current_samples, pattern = "flagstat") & 
		!grepl(x = current_samples, pattern = "idxstats") &
		!grepl(x = current_samples, pattern = "stats")])
}))
keep_samplesheet_df$bam_files <- all_files_list





feature_dataframe <- data.frame(GeneID = as.vector(feature_ranges$hash_id), 
	Chr = as.vector(seqnames(feature_ranges)),
	Start = start(feature_ranges), End = end(feature_ranges),
	Strand = as.vector(strand(feature_ranges)))

param_df <- do.call(rbind, lapply(seq_len(nrow(keep_samplesheet_df)), function(x){
	current_row <- keep_samplesheet_df[x,]
	data.frame(sample_name = stringr::str_split(string = basename(current_row$bam_files), pattern = "\\.", simplify = T)[,1], 
		alignment_file = current_row$bam_files,
		alignment_file_dir = file.path(alignment_dir),
		isPairedEnd = isPairedEnd,
		strand_specificity = strand_specificity)
}))

samples_write_out_df <- param_df
samples_write_out_df$study <- "giulia-unpublished-rbp-kd"
samples_write_out_df$sample_name <- str_split(string = basename(samples_write_out_df$alignment_file), pattern = "\\.", simplify = T)[,1]
samples_write_out_df$samplesheet_path <- "/nemo/lab/patanir/home/users/ziffo/patani-collab/motor-neuron-rbp-knockdown/sample-details/samplesheet.csv"
write.csv(samples_write_out_df, file = file.path(analysis_input_dir, "01-01-public-data-params-dataframe.txt"))

slurm_dir <- file.path(rslurm_dir, "featurecounts_of_exons_and_introns")
if(!dir.exists(slurm_dir)){
	dir.create(slurm_dir, recursive = T)
}
setwd(slurm_dir)
# file.remove(file.path("_rslurm_exon_intron_counts", list.files("./_rslurm_exon_intron_counts/")))
slurm_opts <- list(
	"time" = "72:00:00", 
	"partition" = "ncpu",
	"mem-per-cpu" = "64GB")


process_bam_file <- function(alignment_file, alignment_file_dir, sample_name, 
	feature_df, isPairedEnd, strand_specificity){
	process_file <- file.path(alignment_file_dir, alignment_file)
	count_file <- file.path(alignment_file_dir, paste(alignment_file, 
		"featureCounts", sep = "."))
	an_object <- featureCounts(files = process_file,
		annot.ext = feature_df, 
		useMetaFeatures = T, 
		minMQS = 255,
		isPairedEnd = isPairedEnd,
		strandSpecific = strand_specificity,
		allowMultiOverlap = TRUE,
		requireBothEndsMapped = isPairedEnd,
		reportReads = NULL,
		fraction = TRUE,
		reportReadsPath = alignment_file_dir)
}

slurm_job_list <- slurm_apply(f = process_bam_file, 
		params = param_df, jobname = "answerals_fractional_counts",
		global_objects = c("process_bam_file"),
		feature_df = feature_dataframe,
		pkgs = c("Rsubread"),
		submit = TRUE, cpus_per_node = 1, nodes = 50, 
		slurm_options = slurm_opts)
saveRDS(object = slurm_job_list, file = file.path(analysis_save_dir, "01-01-submit-jobs-for-exonic-intronic-fractional-counts.rds"))