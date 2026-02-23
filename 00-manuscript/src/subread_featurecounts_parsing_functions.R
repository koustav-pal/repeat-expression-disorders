.check_featurecounts_out_list <- function(job_result_list){
	if(!all(c("counts", "targets", "annotation", "stat") %in% names(job_result_list))){
		stop("Object is not originating from Rsubread::featurecounts function")
	}
	return(T)
}
.return_sample_sheet_with_cols <- function(job_result_list, sample_name, samplesheet_df){
	.check_featurecounts_out_list(job_result_list)
	current_colnames <- job_result_list$targets
	which_rows_in_sheet <- vapply(current_colnames, function(alignment_file){
		if(alignment_file %in% samplesheet_df$alignment_file){
			return(which(samplesheet_df$alignment_file == alignment_file & samplesheet_df$sample_name == sample_name))	
		}else{
			return(which(grepl(x = samplesheet_df$alignment_file, pattern = alignment_file) & samplesheet_df$sample_name == sample_name))	
		}
	},1)
	column_ids <- paste(samplesheet_df$sample_name[which_rows_in_sheet], 
		samplesheet_df$alignment_file[which_rows_in_sheet], sep = "_")
	samplesheet_subset_df <- samplesheet_df[which_rows_in_sheet,]
	samplesheet_subset_df$column_id <- column_ids
	return(samplesheet_subset_df)
}
.return_counts_matrix <- function(job_result_list, sample_name, samplesheet_df){
	.check_featurecounts_out_list(job_result_list)
	counts_matrix <- job_result_list$counts
	new_samplesheet_df <- .return_sample_sheet_with_cols(job_result_list, sample_name, samplesheet_df)
	colnames(counts_matrix) <- new_samplesheet_df$column_id
	return(counts_matrix)
}

.return_stat_matrix <- function(job_result_list, sample_name, samplesheet_df){
	.check_featurecounts_out_list(job_result_list)
	stat_matrix <- job_result_list$stat
	new_samplesheet_df <- .return_sample_sheet_with_cols(job_result_list, sample_name, samplesheet_df)
	colnames(stat_matrix) <- new_samplesheet_df$column_id
	return(stat_matrix)
}

.deconvolute_stats_and_counts <- function(results_list){
	count_matrix <- do.call(cbind, lapply(results_list, function(y){
		y$counts
	}))
	stat_matrix <- do.call(cbind, lapply(results_list, function(y){
		y$stat
	}))
	return_sample_df <- do.call(rbind, lapply(results_list, function(y){
		y$sample_df
	}))
	return_list <- list(sample_df = return_sample_df, 
		counts = count_matrix,
		stat = stat_matrix)
	return(return_list)
}
assemble_all_counts_from_jobs <- function(slurm_job_list, slurm_dir, samplesheet_df){
	all_iteration_list <- lapply(seq(from = 0, to = slurm_job_list$nodes - 1), function(x){
		result_file <- file.path(slurm_dir, 
			paste("", "rslurm", slurm_job_list$jobname, sep = "_"), 
			paste(paste("results", x, sep = "_"), "RDS", sep = "."))
		if(!file.exists(result_file)){
			stop(paste("Did not find relevant results for job. Error in node ", x))
		}
		results_out_list <- readRDS(result_file)
		all_results_list <- lapply(names(results_out_list)[names(results_out_list) %in% samplesheet_df$sample_name], function(sample_name){
			current_out_list <- results_out_list[[sample_name]]
			return_sample_df <- .return_sample_sheet_with_cols(current_out_list, sample_name, samplesheet_df)
			return_counts_matrix <- .return_counts_matrix(current_out_list, sample_name, samplesheet_df)
			return_stats_matrix <- .return_stat_matrix(current_out_list, sample_name, samplesheet_df)
			return_list <- list(sample_df = return_sample_df, 
				counts = return_counts_matrix,
				stat = return_stats_matrix)
			return(return_list)
		})
		return_list <- .deconvolute_stats_and_counts(all_results_list)
		return(return_list)
	})
	return_list <- .deconvolute_stats_and_counts(all_iteration_list)
	return(return_list)
}