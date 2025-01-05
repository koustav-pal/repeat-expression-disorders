data_dir <- normalizePath("/camp/lab/luscomben/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "05-simulated-reads-differential-profiling")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
rslurm_dir <- file.path(analysis_dir, "rslurm")
alignment_dir <- file.path(project_dir, "analysis", "Alignments")
# ==========================================================================
# Packages and Data
# ==========================================================================
require("rtracklayer")
require("GenomicFeatures")
require("GenomeInfDb")
require("GenomicRanges")
require("parallel")
require("Rsamtools")
require("QuasR")
require("digest")


genome_file <- "/camp/home/palk/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fasta_index <- indexFa(genome_file)
fa_index_read <- scanFaIndex(file = genome_file)

chom_sizes_df <- data.frame(chrom = as.character(seqnames(fa_index_read)),
length = width(fa_index_read),
is_circular = F)

chom_sizes_df <- chom_sizes_df[chom_sizes_df$chrom %in% c(1:22, "X", "Y", "MT"),]
chom_sizes_df$chr_starts <- cumsum(as.numeric(chom_sizes_df$length)) - chom_sizes_df$length

annotation_gtf <- "/camp/home/palk/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
annotation_ranges <- import.gff2(annotation_gtf)
annotation_ranges <- dropSeqlevels(annotation_ranges, seqlevels(annotation_ranges)[!(seqlevels(annotation_ranges) %in% chom_sizes_df$chrom)], pruning.mode = "coarse")
txdb <- makeTxDbFromGRanges(annotation_ranges)

intron_ranges <- intronicParts(txdb, linked.to.single.gene.only=TRUE)
exon_ranges <- exonicParts(txdb, linked.to.single.gene.only=TRUE)

intron_ranges$type <- "intron"
exon_ranges$type <- "exon"


feature_ranges <- c(exon_ranges, intron_ranges)

feature_starts <- start(feature_ranges) + chom_sizes_df$chr_starts[match(as.vector(seqnames(feature_ranges)), chom_sizes_df$chrom)]
feature_ends <- end(feature_ranges) + chom_sizes_df$chr_starts[match(as.vector(seqnames(feature_ranges)), chom_sizes_df$chrom)]

sorted_feature_ranges <- feature_ranges[order(feature_starts)]

sorted_feature_ranges$feature_id <- paste(sorted_feature_ranges$gene_id, 
as.vector(seqnames(sorted_feature_ranges)), start(sorted_feature_ranges), 
end(sorted_feature_ranges), as.vector(strand(sorted_feature_ranges)), sorted_feature_ranges$type, sep = ":")

sorted_feature_ranges$hash_id <- vapply(seq_along(sorted_feature_ranges), digest, "1", algo = "crc32")


if(any(duplicated(sorted_feature_ranges$hash_id))){
    stop("Duplicated hash_ids found!")
}

saveRDS(object = sorted_feature_ranges, file = file.path(analysis_save_dir, "01-all-merged-exonic-and-intronic-segments-in-genes.rds"))