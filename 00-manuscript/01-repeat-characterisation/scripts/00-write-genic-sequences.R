data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "01-repeat-characterisation")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
alignment_dir <- file.path(project_dir, "analysis", "Alignments")
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rtracklayer")
require("GenomicRanges")
require("BSgenome.Hsapiens.UCSC.hg38")
require("Biostrings")
require("GenomeInfoDb")

annotation_gtf <- "/camp/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
annotation_ranges <- import.gff2(annotation_gtf)

genes_ranges <- annotation_ranges[annotation_ranges$type == "gene"]

current_seqlevels <- seqlevels(genes_ranges)
map_seqlevels <- mapSeqlevels(current_seqlevels, "ucsc")

genes_ranges <- dropSeqlevels(genes_ranges, 
	value = current_seqlevels[current_seqlevels %in% names(map_seqlevels)[is.na(map_seqlevels)]],
	pruning.mode = "coarse")
seqlevels(genes_ranges) <- map_seqlevels[!is.na(map_seqlevels)]

sequence_list <- getSeq(BSgenome.Hsapiens.UCSC.hg38, genes_ranges)

names(sequence_list) <- paste(genes_ranges$gene_id, as.vector(seqnames(genes_ranges)), 
		start(genes_ranges), end(genes_ranges), strand(genes_ranges), sep = ":")

writeXStringSet(x = sequence_list, filepath = file.path(analysis_input_dir, "00-all_genic_sequences_chr1-22.fa"))