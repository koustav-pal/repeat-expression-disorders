return_cache_filepath <- function(cache_object, query_name, meta_name){
    the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
    if(nrow(the_query) > 0){
        if(!file.exists(the_query$rpath)){
            new_rpath <- file.path(normalizePath(bfccache(cache_object)), basename(the_query$rpath))
            if(!file.exists(new_rpath)){
                file.create(new_rpath)
            }
            file_reset_status <- reset_cache_filepath(cache_object, query_name)
            if(file_reset_status){
                the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
                return(the_query$rpath)
            }
        }else{
            file_status <- reset_cache_filepath(cache_object, query_name = query_name)
            if(file_status){
                the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
                return(the_query$rpath)
            }else{
                stop(paste(the_query$rpath, "does not exist!"))
            }
        }
    }else{
        new_path <- bfcnew(cache_object, rname = query_name, rtype = "local", ext = ".txt")
        meta_df <- data.frame(list(rid = names(new_path), seed = NA, iteration = NA, complete = F))
        if(meta_name %in% bfcmetalist(cache_object)){
            bfcmeta(cache_object, name = meta_name, append = T) <- meta_df
        }else{
            bfcmeta(cache_object, name = meta_name, overwrite = T) <- meta_df
        }
        return(new_path)
    }
}

reset_cache_filepath <- function(cache_object, query_name){
    the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
    rpath <- the_query$rpath
    # fpath <- the_query$fpath
    if(!file.exists(rpath)){
        new_rpath <- file.path(normalizePath(bfccache(cache_object)), basename(rpath))
        # new_fpath <- file.path(normalizePath(bfccache(cache_object)), basename(fpath))
        # cache_object[[the_query$rid]] <- new_rpath
        bfcupdate(cache_object, rids = the_query$rid, rname = query_name, rpath = new_rpath)
    }
    the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
    return(file.exists(the_query$rpath))
}

return_cache_seed <- function(cache_object, query_name, meta_name){
    the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
    return_df <- bfcmeta(cache_object, name = meta_name)
    current_df <- return_df[return_df$rid == the_query$rid,]
    if(is.na(current_df$seed)){
        current_seed <- sample(x = seq(1, to = .Machine$integer.max), size = 1)
        bfcmetaremove(x = cache_object, name = meta_name)
        return_df$seed[return_df$rid == the_query$rid] <- current_seed
        bfcmeta(cache_object, name = meta_name, overwrite = T) <- return_df
    }else{
        current_seed <- current_df$seed
    }
    return(current_seed)
}

update_bfcmeta_iteration <- function(cache_object, query_name, iteration, meta_name){
    the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
    return_df <- bfcmeta(cache_object, name = meta_name)
    bfcmetaremove(x = cache_object, name = meta_name)
    return_df$iteration[return_df$rid == the_query$rid] <- iteration
    bfcmeta(cache_object, name = meta_name, overwrite = T) <- return_df
    return(TRUE)
}

mark_bfcmeta_complete <- function(cache_object, query_name, meta_name){
    the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
    return_df <- bfcmeta(cache_object, name = meta_name)
    bfcmetaremove(x = cache_object, name = meta_name)
    return_df$complete[return_df$rid == the_query$rid] <- T
    bfcmeta(cache_object, name = meta_name, overwrite = T) <- return_df
    return(TRUE)
}

check_bfcmeta_complete <- function(cache_object, query_name, meta_name){
    the_query <- bfcquery(x = cache_object, query = query_name, field = "rname", exact = T)
    return_df <- bfcmeta(cache_object, name = meta_name)
    current_df <- return_df[return_df$rid == the_query$rid,]
    if(as.logical(current_df$complete)){
        return(T)
    }else{
        return(F)
    }
}

binomial_generated_counts <- function(a_matrix){
    col_sums <- colSums(a_matrix)
    prop_matrix <- t(t(a_matrix) / col_sums)
    sizes_vector <- as.vector(rep(col_sums, each = nrow(prop_matrix)))
    col_means_vector <- rep(colMeans(a_matrix), each = nrow(a_matrix))

    seed_string <- attributes(a_matrix)$seed
    set.seed(as.integer(seed_string))

    simulated_counts_matrix <- matrix(rbinom(n = length(prop_matrix), 
        size = as.integer(sizes_vector), p = as.vector(prop_matrix)),
        nrow = nrow(a_matrix), ncol = ncol(a_matrix))
    
    colnames(simulated_counts_matrix) <- colnames(a_matrix)
    rownames(simulated_counts_matrix) <- rownames(a_matrix)

    attr(simulated_counts_matrix, "seed") <- seed_string
    return(simulated_counts_matrix)
}

perform_deseq_de_profiling <- function(a_matrix, ...){
    args <- list(...)
    coldata <- args$coldata
    contrast_df <- args$contrasts
    design <- as.formula(paste("~", paste(args$design_vars, collapse = "*"), sep = " ")) 
    pvalue_threshold <- ifelse("pvalue_threshold" %in% names(args), args$pvalue_threshold, 0.05)
    feature_list <- args$feature_list
    if("prior_filter" %in% names(args)){
        prior_filters <- args$prior_filter
    }else{
        prior_filters <- rep(T, nrow(a_matrix))
    }

    # message(length(prior_filters))
    if(length(prior_filters) != nrow(a_matrix)){
        stop("Prior filter length must be equivalent to matrix nrows")
    }
    keep_informative_rows <- rowMeans(a_matrix) > 1
    message("Analysing ", length(keep_informative_rows[keep_informative_rows]), " features...")
    sig_diff_matrix_list <- lapply(feature_list, function(features){
        keep_features <- rownames(a_matrix) %in% features
        a_deseq_object <- DESeqDataSetFromMatrix(countData = a_matrix[keep_informative_rows & prior_filters & keep_features,], 
            colData = coldata, design = design)
        de_results <- DESeq(a_deseq_object)
        contrast_list <- build_contrasts(de_results, contrast_df)
        return_list <- lapply(names(contrast_list), function(a_contrast_name){
            differential_results = results(de_results, contrast = contrast_list[[a_contrast_name]], 
                independentFiltering = T, pAdjustMethod = "BH")
            feature_diff_matrix <- differential_results[
                differential_results$pvalue < pvalue_threshold & 
                !is.na(differential_results$pvalue),]
            feature_diff_matrix$contrast <- a_contrast_name
            return(feature_diff_matrix)
        })
        return_df <- do.call(rbind, return_list)
    })
    sig_diff_matrix <- do.call(rbind, sig_diff_matrix_list)
    return(sig_diff_matrix)
}

build_contrasts <- function(deseq_object, contrast_df, 
    contrast_group = "contrast_name", within_contrast_group = "group"){

    col_data <- colData(deseq_object)
    model_matrix <- model.matrix(design(deseq_object), col_data)

    contrast_df_split <- split(contrast_df, contrast_df[,contrast_group])
    return_contrast_list <- lapply(contrast_df_split, function(a_contrast_df){
        group_split_df <- split(a_contrast_df, a_contrast_df[,within_contrast_group])
        all_comparison_list <- lapply(group_split_df, function(group_df){
            from_filter_list <- lapply(seq_along(group_df[,"from_var"]), function(x){
                col_data[,group_df[x,"from_var"]] == group_df[x,"from_value"]
            })
            from_filter <- Reduce("*", from_filter_list) == 1
            to_filter_list <- lapply(seq_along(group_df[,"to_var"]), function(x){
                col_data[,group_df[x,"to_var"]] == group_df[x,"to_value"]
            })
            to_filter <- Reduce("*", to_filter_list) == 1
            return(colMeans(model_matrix[from_filter,]) - colMeans(model_matrix[to_filter,]))
        })
        return_contrast <- Reduce("-", all_comparison_list)
        return(return_contrast)
    })
    return(return_contrast_list)
}

assemble_deseq_de_results <- function(a_deseq_result_matrix){
    return_df <- data.frame(feature_ids = rownames(a_deseq_result_matrix),
        contrast = a_deseq_result_matrix$contrast,
        lfc = a_deseq_result_matrix$log2FoldChange,
        padj = a_deseq_result_matrix$padj, 
        pvalue = a_deseq_result_matrix$pvalue)
    return(return_df)
}

generate_iterated_deseq_de_result <- function(num_iterations, original_matrix, ...){
    args <- list(...)
    seed_present <- F
    cache_dir <- args$biocfilecache
    sample_name <- args$sample_name
    if("prior_seeds" %in% names(attributes(original_matrix))){
        prior_seeds <- attributes(original_matrix)$prior_seeds
        if(num_iterations != length(prior_seeds)){
            stop("Length of attached seed list does not match that of iterations! Seed list length is ",length(prior_seeds))
        }
        seed_present <- T
    }
    cache_object <- BiocFileCache(file.path(cache_dir, sample_name), ask = FALSE)
    a_result_list <- lapply(seq(from = 1, to = num_iterations), function(x){
        cache_outpath <- return_cache_filepath(cache_object, 
            query_name = paste(sample_name, x, sep = "_"),
            meta_name = "ResourceData")
        if(check_bfcmeta_complete(cache_object, 
            query_name = paste(sample_name, x, sep = "_"),
            meta_name = "ResourceData")){
            return(cache_outpath)
        }
        current_cache_seed <- return_cache_seed(cache_object, 
            query_name = paste(sample_name, x, sep = "_"),
            meta_name = "ResourceData")
        cache_seed_present <- T
        if(seed_present){
            current_seed <- prior_seeds[x]
            attr(original_matrix, "seed") <- current_seed
        }
        if(cache_seed_present){
            attr(original_matrix, "seed") <- current_cache_seed
        }
        count_matrix <- binomial_generated_counts(original_matrix)
        deseq_de_results <- perform_deseq_de_profiling(count_matrix, ...)
        deseq_de_df <- assemble_deseq_de_results(deseq_de_results)
        deseq_de_df$iterator <- x
        deseq_de_df$seed <- attributes(count_matrix)$seed
        fwrite(x = deseq_de_df, quote = F, col.names = T, file = cache_outpath,
            sep = "\t", append = F, row.names = F)
        update_bfcmeta_iteration(cache_object, 
            query_name = paste(sample_name, x, sep = "_"),
            iteration = x,
            meta_name = "ResourceData")
        mark_bfcmeta_complete(cache_object, 
            query_name = paste(sample_name, x, sep = "_"),
            meta_name = "ResourceData")
        return(cache_outpath)
    })
    deseq_results_path <- do.call(c, a_result_list)
    return(deseq_results_path)
}

fetch_all_files <- function(cache_dir, sample_name, num_iterations){
    cache_object <- BiocFileCache(cache_dir, ask = FALSE)
    the_query <- bfcquery(x = cache_object, query = sample_name, field = "rname", exact = F)
    all_files_df <- do.call(rbind, lapply(seq(from = 1, to = num_iterations), function(x){
        data.frame(sample_name = sample_name, 
            iteration = x, filepath = the_query$rpath[the_query$rname == paste(sample_name, x, sep = "_")],
            complete = the_query$complete[the_query$rname == paste(sample_name, x, sep = "_")])
    }))
    return(all_files_df)
}

fetch_results_matrix <- function(original_matrix, cache_dir, 
    sample_name, num_iterations){
    message("fetching database...")
    current_files_df <- fetch_all_files(cache_dir, sample_name, num_iterations)
    if(any(current_files_df$complete != 1) | nrow(current_files_df) != num_iterations){
        stop("Some iterations are not complete")
    }
    
    a_vector <- rep(NA, nrow(original_matrix))
    names(a_vector) <- rownames(original_matrix)

    message("fetching values...")
    
    all_iter_values_list <- lapply(seq_len(nrow(current_files_df)), function(x){
        current_df <- fread(current_files_df[x,"filepath"], data.table = F)
        current_df <- current_df[current_df$feature_ids %in% rownames(original_matrix),]
        lfcs_vector <- a_vector
        lfcs_vector[current_df$feature_ids] <- current_df$lfc
        padj_vector <- a_vector
        padj_vector[current_df$feature_ids] <- current_df$padj
        pval_vector <- a_vector
        pval_vector[current_df$feature_ids] <- current_df$pval
        list(lfcs = lfcs_vector,
            padj = padj_vector,
            pval = pval_vector)
    })

    message("assembling values...")

    lfcs_matrix_list <- lapply(all_iter_values_list, function(a_list){
        a_list[["lfcs"]]
    })
    lfcs_matrix <- do.call(cbind, lfcs_matrix_list)

    pval_matrix_list <- lapply(all_iter_values_list, function(a_list){
        a_list[["pval"]]
    })
    pval_matrix <- do.call(cbind, pval_matrix_list)
    
    padj_matrix_list <- lapply(all_iter_values_list, function(a_list){
        a_list[["padj"]]
    })
    padj_matrix <- do.call(cbind, padj_matrix_list)
    
    which_keep_lfcs_matrix <- rowSums(!is.na(lfcs_matrix)) > 0
    keep_lfcs_matrix <- lfcs_matrix[which_keep_lfcs_matrix,]
    keep_pval_matrix <- pval_matrix[which_keep_lfcs_matrix,]
    keep_padj_matrix <- padj_matrix[which_keep_lfcs_matrix,]
    return(list(lfc = keep_lfcs_matrix, 
        pval = keep_pval_matrix, 
        padj = keep_padj_matrix))
}

estimate_standard_errors <- function(original_matrix, cache_dir, 
    sample_name, num_iterations){

    results_list <- fetch_results_matrix(original_matrix, 
        cache_dir, sample_name, num_iterations)
    lfc_matrix <- results_list[["lfc"]]
    pval_matrix <- results_list[["pval"]]
    padj_matrix <- results_list[["padj"]]

    n <- rowSums(!is.na(lfc_matrix))
    lfc_se <- rowSds(lfc_matrix, na.rm = T) / sqrt(rowSums(!is.na(lfc_matrix)))
    pval_se <- rowSds(pval_matrix, na.rm = T) / sqrt(rowSums(!is.na(pval_matrix)))
    padj_se <- rowSds(padj_matrix, na.rm = T) / sqrt(rowSums(!is.na(padj_matrix)))
    return_df <- data.frame(feature_ids = rownames(lfc_matrix),
        sample_name = sample_name, n = n, lfc_se = lfc_se, 
        pval_se = pval_se, padj_se = padj_se)
    return(return_df)
}

estimate_distances <- function(standard_error_df, return_distances = F) {
    se_summary_df <- standard_error_df %>% 
    group_by(n) %>%
    mutate(mean_lfc_se = mean(lfc_se[!is.na(lfc_se)]),
        mean_pval_se = mean(pval_se[!is.na(pval_se)]),
        mean_padj_se = mean(padj_se[!is.na(padj_se)])) %>%
    select(sample_name, n, mean_lfc_se, mean_pval_se, mean_padj_se) %>%
    unique() %>% ungroup() %>% as.data.frame()

    se_cumsum_df <- se_summary_df %>%
    filter(!is.na(se_summary_df$mean_lfc_se) &
        !is.na(se_summary_df$mean_pval_se) &
        !is.na(se_summary_df$mean_padj_se)) %>%
    arrange(n) %>%
    mutate(lfc_se_cumsum = cumsum(mean_lfc_se),
        pval_se_cumsum = cumsum(mean_pval_se),
        padj_se_cumsum = cumsum(mean_padj_se)) %>%
    select(sample_name, n, lfc_se_cumsum, pval_se_cumsum, padj_se_cumsum) %>% 
    unique() %>% as.data.frame()

    plot_df <- se_cumsum_df[,c(1,2)]
    se <- c(se_cumsum_df$lfc_se_cumsum, se_cumsum_df$pval_se_cumsum, se_cumsum_df$padj_se_cumsum)
    variables <- rep(c("lfc", "pval", "padj"), each = nrow(plot_df))
    return_df <- do.call(rbind, lapply(unique(variables), function(a_var){
        var_temp_df <- plot_df
        var_temp_df$se <- (se[variables == a_var])
        var_temp_df$se_scaled <- var_temp_df$se / max(var_temp_df$se)
        var_temp_df$var <- variables[variables == a_var]
        return(var_temp_df)
    }))
    if(!return_distances){
        return(return_df)
    }

    distance_df <- do.call(rbind, lapply(unique(return_df$var), function(x){
        var_df <- return_df[return_df$var == x,]

        x_coords <- seq(from = min(var_df$n), to = max(var_df$n))
        y_coords <- seq(from = min(var_df$se_scaled), 
            to = max(var_df$se_scaled), length.out = length(x_coords))

        var_df$distances <- vapply(seq_len(nrow(var_df)), function(z){
            a_row <- var_df[z,]
            distances <- sqrt((a_row$n - x_coords)^2 + (a_row$se_scaled - y_coords)^2)
            which_min <- which.min(distances)
            distances[which_min]
        },1)
        n_xcoords <- vapply(seq_len(nrow(var_df)), function(z){
            a_row <- var_df[z,]
            distances <- sqrt((a_row$n - x_coords)^2 + (a_row$se_scaled - y_coords)^2)
            which_min <- which.min(distances)
            x_coords[which_min]
        },1)
        se_ycoords <- vapply(seq_len(nrow(var_df)), function(z){
            a_row <- var_df[z,]
            distances <- sqrt((a_row$n - x_coords)^2 + (a_row$se_scaled - y_coords)^2)
            which_min <- which.min(distances)
            y_coords[which_min]
        },1)
        dist_max <- which.max(var_df$distances)
        var_df$group <- "data"
        temp_df <- var_df
        temp_df$n <- n_xcoords
        temp_df$se_scaled <- se_ycoords
        temp_df$group <- "model"
        return(rbind(var_df[dist_max,], temp_df[dist_max,]))
    }))
    return(distance_df)
}

report_ci <- function(a_matrix, ci = 0.95, var){
    feature_ids <- rownames(a_matrix)
    all_sds <- rowSds(a_matrix, na.rm = T) 
    all_n <- rowSums(!is.na(a_matrix))
    all_means <- rowMeans(a_matrix, na.rm = T)

    error <- (qnorm(1 - ((1 - ci)/2))*all_sds)/sqrt(all_n)
    return(data.frame(group = var, feature_id = feature_ids, 
            mean = all_means, sd = all_sds, se = error, n = all_n))
}

get_final_results <- function(original_matrix, cache_dir, 
    sample_name, num_iterations, ci = 0.95){

    results_list <- fetch_results_matrix(original_matrix, 
        cache_dir, sample_name, num_iterations)
    lfc_matrix <- results_list[["lfc"]]
    pval_matrix <- results_list[["pval"]]
    padj_matrix <- results_list[["padj"]]

    lfc_stats <- report_ci(lfc_matrix, ci = ci, var = "lfc")
    pval_stats <- report_ci(pval_matrix, ci = ci, var = "pval")
    padj_stats <- report_ci(padj_matrix, ci = ci, var = "padj")

    assembled_results_df <- rbind(lfc_stats, pval_stats, padj_stats)
    return(assembled_results_df)
}

filter_results <- function(results, distances){
    threshold_df <- distances[distances$group == "data",]
    n_thresholds <- threshold_df$n
    names(n_thresholds) <- threshold_df$var

    n_thresholds_rep <- n_thresholds[match(results$group, names(n_thresholds))]

    results_split <- split(results$n > n_thresholds_rep, results$feature_id)

    filter <- vapply(results_split, function(x){ all(x) },T)

    de_features <- names(results_split)[filter]

    filtered_results_df <- results[results$feature_id %in% de_features,]
    rownames(filtered_results_df) <- NULL

    return(filtered_results_df)
}