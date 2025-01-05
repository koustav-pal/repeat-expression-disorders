# ==========================================================================
# Intial script setup
# ==========================================================================
data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "13-03-public-data")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
alignment_dir <- "/nemo/project/proj-luscombn-patani/working/public-data"
feature_ranges <- readRDS("../05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
samples <- c(
	"c9orf72-dafinca-2020" = "ipsc-mn-c9orf72-dafinca-2020",
	"c9orf72-catanese-2021" = "ipsc-mn-c9orf72-fus-catanese-2021",
	"c9orf72-sareen-2013" = "ipsc-mn-c9orf72-sareen-2013",
	"c9orf72-sommer-2022" = "ipsc-mn-c9orf72-sommer-2022",
	"c9orf72-sterneckert-2020" = "ipsc-mn-c9orf72-sterneckert-2020",
	"c9orf72-casteli-2021" = "motor-neuron-c9orf72-ipsc-casteli-2021",
	"fus-desantis-2017" = "ipsc-mn-fus-desantis-2017",
	"fus-hawkins-2022" = "ipsc-mn-fus-hawkins-2022",
	"fus-kapeli-2016" = "ipsc-mn-fus-kapeli-2016",
	"hnrnpa2b1-markmiller-2021" = "ipsc-mn-hnrnpa2b1-markmiller-2021",
	"sod1-bhinge-2017" = "ipsc-mn-sod1-bhinge-2017",
	"sod1-moccia-2014" = "ipsc-mn-sod1-moccia-2014",
	"sod1-tardbp-dash-2022" = "ipsc-mn-sod1-tardbp-dash-2022",
	"sod1-wang-2017" = "ipsc-mn-sod1-wang-2017",
	"tardbp-dafinca-2020" = "ipsc-mn-tardbp-dafinca-2020",
	"tardbp-kd-brown-2022" = "ipsc-mn-tardbp-kd-brown-2022",
	"tardbp-kd-klim-2022" = "ipsc-mn-tardbp-kd-klim-2019",
	"tardbp-melamed-2019" = "ipsc-mn-tardbp-melamed-2019",
	"tardbp-smith-2021" = "ipsc-mn-tardbp-smith-2021",
	"tardbp-kd-kapeli-2016" = "ipsc-mn-tdp43-kd-kapeli-2016",
	"vcp-luisier-2018" = "ipsc-mn-vcp-luisier-2017",
	"vcp-ziff-2023" = "motor-neuron-vcpinhibitor-ipsc-ziff-2023")
samples_bam_dir <- c(
	"c9orf72-dafinca-2020" = "star_salmon",
	"c9orf72-catanese-2021" = "star_salmon",
	"c9orf72-sareen-2013" = "star_salmon",
	"c9orf72-sommer-2022" = "star_salmon",
	"c9orf72-sterneckert-2020" = "star_salmon",
	"c9orf72-casteli-2021" = "star_salmon",
	"fus-desantis-2017" = "star_salmon",
	"fus-hawkins-2022" = "star_salmon",
	"fus-kapeli-2016" = "star_salmon",
	"hnrnpa2b1-markmiller-2021" = "star_salmon",
	"sod1-bhinge-2017" = "star_salmon",
	"sod1-moccia-2014" = "star_salmon",
	"sod1-tardbp-dash-2022" = "star_salmon",
	"sod1-wang-2017" = "star_salmon",
	"tardbp-dafinca-2020" = "star_salmon",
	"tardbp-kd-brown-2022" = "star_salmon",
	"tardbp-kd-klim-2022" = "star_salmon",
	"tardbp-melamed-2019" = "star_salmon",
	"tardbp-smith-2021" = "star_salmon",
	"tardbp-kd-kapeli-2016" = "star_salmon",
	"vcp-luisier-2018" = "star_rsem",
	"vcp-ziff-2023" = "star_rsem")
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
isPairedEnd <- vapply(samples, function(x){
	current_sample_samplesheet <- file.path(alignment_dir, x, "sample-details", "samplesheet.csv")
	sample_df <- read.csv(current_sample_samplesheet)
	if(!is.null(sample_df$library_layout)){
		return(all(sample_df$library_layout == "PAIRED"))
	}else{
		nfcore_df <- read.csv(file.path(alignment_dir, x, "sample-details", "nfcore_sampleinfo.csv"))
		return(all(!is.na(nfcore_df$fastq_1) & !is.na(nfcore_df$fastq_2)))
	}
	return(FALSE)
},T)

strand_specificity <- vapply(samples, function(x){
	message(x)
	current_sample_samplesheet <- file.path(alignment_dir, x, "sample-details", "samplesheet.csv")
	sample_df <- read.csv(current_sample_samplesheet)
	if("strandedness" %in% colnames(sample_df)){
	}else{
		sample_df <- read.csv(file.path(alignment_dir, x, "sample-details", "nfcore_sampleinfo.csv"))
	}
	strandedness <- unique(sample_df$strandedness)
	ifelse(strandedness == "reverse", 2, ifelse(strandedness == "forward", 1, 0))
},1)


feature_dataframe <- data.frame(GeneID = as.vector(feature_ranges$hash_id), 
	Chr = as.vector(seqnames(feature_ranges)),
	Start = start(feature_ranges), End = end(feature_ranges),
	Strand = as.vector(strand(feature_ranges)))

param_df <- do.call(rbind, lapply(names(samples), function(sample_name){
	message(sample_name)
	data.frame(sample_name = sample_name, 
		alignment_file = str_replace_all(list.files(file.path(alignment_dir, samples[sample_name], "nfcore", samples_bam_dir[sample_name]), pattern = ".bam.bai"), pattern = ".bai", replacement = ""),
		alignment_file_dir = file.path(alignment_dir, samples[sample_name], "nfcore", samples_bam_dir[sample_name]),
		isPairedEnd = isPairedEnd[sample_name],
		strand_specificity = strand_specificity[sample_name])
}))

samples_write_out_df <- param_df
samples_write_out_df$samplesheet_path <- file.path(alignment_dir, samples[param_df$sample_name], "sample-details", "samplesheet.csv")
write.csv(samples_write_out_df, file = file.path(analysis_input_dir, "01-01-public-data-params-dataframe.txt"))




slurm_dir <- file.path(rslurm_dir, "featurecounts_of_exons_and_introns")
if(!dir.exists(slurm_dir)){
	dir.create(slurm_dir, recursive = T)
}
setwd(slurm_dir)
# file.remove(file.path("_rslurm_exon_intron_counts", list.files("./_rslurm_exon_intron_counts/")))
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
		params = param_df, jobname = "exon_intron_counts",
		global_objects = c("process_bam_file"),
		feature_df = feature_dataframe,
		pkgs = c("Rsubread"),
		submit = FALSE, cpus_per_node = 1, nodes = nrow(param_df), 
		slurm_options = slurm_opts)
saveRDS(object = slurm_job_list, file = file.path(analysis_save_dir, "01-submit-jobs-for-exonic-intronic-counts.rds"))