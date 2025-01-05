# ==========================================================================
# Intial script setup
# ==========================================================================
data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "13-04-tyzack-et-al-data-processing")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
alignment_dir <- file.path(project_dir, "analysis", "Alignments")
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
samples <- c("vcp-tyzack-2021" = "motor-neuron-vcp-ipsc-tyzack-2021")
isPairedEnd <- c("vcp-tyzack-2021" = TRUE)
strand_specificity <- c("vcp-tyzack-2021" = 2)
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

param_df <- do.call(rbind, lapply(names(samples), function(sample_name){
	data.frame(sample_name = sample_name, 
		alignment_file = str_replace_all(list.files(file.path(alignment_dir, samples[sample_name], "star_rsem"), pattern = ".bam.bai"), pattern = ".bai", replacement = ""),
		alignment_file_dir = file.path(alignment_dir, samples[sample_name], "star_rsem"),
		isPairedEnd = isPairedEnd[sample_name],
		strand_specificity = strand_specificity[sample_name])
}))



samples_write_out_df <- param_df
samples_write_out_df$samplesheet_path <- NA
write.csv(samples_write_out_df, file = file.path(analysis_input_dir, "01-01-public-data-params-dataframe.txt"))



slurm_dir <- file.path(rslurm_dir, "featurecounts_of_exons_and_introns")
if(!dir.exists(slurm_dir)){
	dir.create(slurm_dir, recursive = T)
}
setwd(slurm_dir)
file.remove(file.path("_rslurm_exon_intron_counts", list.files("./_rslurm_exon_intron_counts/")))
slurm_opts <- list(
	"time" = "72:00:00", 
	"partition" = "cpu",
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
		params = param_df, jobname = "tyzacketal_fractional_counts",
		global_objects = c("process_bam_file"),
		feature_df = feature_dataframe,
		pkgs = c("Rsubread"),
		submit = TRUE, cpus_per_node = 1, nodes = 96, 
		slurm_options = slurm_opts)
saveRDS(object = slurm_job_list, file = file.path(analysis_save_dir, "01-01-submit-jobs-for-exonic-intronic-fractional-counts.rds"))