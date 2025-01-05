# ==========================================================================
# Intial script setup
# ==========================================================================
data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "13-nygc-consortium")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
isPairedEnd <- T
strand_specificity <- 2
sampleinfo_df <- read.csv(file.path(analysis_input_dir, "00_filtered_sample_gsm_list.txt"))
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
feature_dataframe <- data.frame(GeneID = as.vector(feature_ranges$hash_id), 
	Chr = as.vector(seqnames(feature_ranges)),
	Start = start(feature_ranges), End = end(feature_ranges),
	Strand = as.vector(strand(feature_ranges)))

param_df <- data.frame(sample_name = sampleinfo_df[,14], 
	alignment_file = file.path(sampleinfo_df$bam_file),
	alignment_file_dir = file.path(sampleinfo_df$alignment_dir),
	isPairedEnd = T,
	strand_specificity = 2)


param_df <- unique(param_df)
samples_write_out_df <- param_df
samples_write_out_df$study <- "NYGC"
samples_write_out_df$samplesheet_path <- file.path(analysis_input_dir, "00_filtered_sample_gsm_list.txt")
write.csv(samples_write_out_df, file = file.path(analysis_input_dir, "01-01-nygc-params-dataframe.txt"))



slurm_dir <- file.path(rslurm_dir, "featurecounts_of_exons_and_introns")
if(!dir.exists(slurm_dir)){
	dir.create(slurm_dir, recursive = T)
}
setwd(slurm_dir)
# file.remove(file.path("_rslurm_exon_intron_counts", list.files("./_rslurm_exon_intron_counts/")))
slurm_opts <- list(
	"time" = "72:00:00", 
	"partition" = "cpu",
	"mem-per-cpu" = "8GB")


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

# boo <- process_bam_file(alignment_file = "GSM4659746.Aligned.sortedByCoord.out.bam",
# 	alignment_file_dir = "/nemo/project/proj-patanir-seq/palk/0103-NYGC-consortium/output-files/GSM4659746",
# 	sample_name = "SRX8681220",
# 	feature_df = feature_dataframe,
# 	isPairedEnd = T,
# 	strand_specificity = 2)

slurm_job_list <- slurm_apply(f = process_bam_file, 
		params = param_df, jobname = "exon_intron_fractional_counts",
		global_objects = c("process_bam_file"),
		feature_df = feature_dataframe,
		pkgs = c("Rsubread"),
		submit = TRUE, cpus_per_node = 8, nodes = 250, 
		slurm_options = slurm_opts)
saveRDS(object = slurm_job_list, file = file.path(analysis_save_dir, "03-submit-jobs-for-exonic-intronic-fractional-counts.rds"))