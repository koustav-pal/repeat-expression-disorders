# .load_bigwig <- function(path_to_file){
#     a_ranges_object <- import.bw(path_to_file)
#     # a_ranges_object$mids <- start(a_ranges_object) + 
#     #     ceiling(width(a_ranges_object)/2)
#     # a_new_ranges_object <- MakeGRangesObject(
#     #     Chrom = seqnames(a_ranges_object), 
#     #     Start = a_ranges_object$mids,
#     #     End =  a_ranges_object$mids)
#     # a_new_ranges_object$score <- a_ranges_object$score/width(a_ranges_object)
#     return(a_ranges_object)
# }

.return_overlap_weighted_score_between_ranges <- function(query_ranges, bin_to, subject_ranges){
    query_start <- start(query_ranges)
    query_stop <- end(query_ranges)
    Main_seq <- seq(from = query_start, to = query_stop)
    Main_pos <- Main_seq - query_start
    end_dist_from_start <- end(subject_ranges) - query_start
    start_dist_to_end <- query_stop - start(subject_ranges)

    subject_ranges_in_feature <- subject_ranges[(end_dist_from_start >= 0 & start_dist_to_end >= 0)]
    
    Step_size <- width(query_ranges)/bin_to
    average_score <- rep(0, bin_to)
    if(length(subject_ranges_in_feature) == 0){
        return(average_score)
    }
    score_per_basepair_list <- lapply(seq_along(subject_ranges_in_feature), function(x){
        a_temp_ranges <- subject_ranges_in_feature[x]
        a_start <- start(a_temp_ranges)
        an_end <- end(a_temp_ranges)
        Seq <- seq(from = a_start, to = an_end)
        position_offset_by_start <- Seq - query_start
        if(any(position_offset_by_start < 0)){
            position_offset_by_start <- position_offset_by_start[position_offset_by_start >= 0]
        }
        if(any(position_offset_by_start > max(Main_pos))){
            position_offset_by_start <- position_offset_by_start[position_offset_by_start <= max(Main_pos)]
        }
        data.frame(pos = position_offset_by_start, 
            score = 1, 
            frag = x)
    })
    score_per_basepair_df <- do.call(rbind, score_per_basepair_list)
    if(any(!(Main_pos %in% score_per_basepair_df$pos))){
        temp_df <- rbind(score_per_basepair_df,
            data.frame(pos = Main_pos[!(Main_pos %in% score_per_basepair_df$pos)],
                score = 0, frag = length(subject_ranges_in_feature)+1))
        score_per_basepair_df <- rbind(score_per_basepair_df, temp_df)
        score_per_basepair_df <- score_per_basepair_df[order(score_per_basepair_df$pos),]
    }
    score_per_basepair_df$pos[score_per_basepair_df$pos == 0] <- score_per_basepair_df$pos[score_per_basepair_df$pos == 0] + 1
    score_per_basepair_df$bins <- ceiling(score_per_basepair_df$pos/Step_size)
    Score_split <- split(score_per_basepair_df$score, score_per_basepair_df$bins)
    frag_split <-  split(score_per_basepair_df$frag, score_per_basepair_df$bins)
    average_score <- vapply(seq_along(Score_split), function(x){
        frag_table <- as.vector(table(frag_split[[x]]))
        scores <- vapply(split(Score_split[[x]], frag_split[[x]]), unique, 1)
        weighted.mean(scores, frag_table)
    }, 1)
    return(average_score)
}

.return_bigwig_signal <- function(query_ranges, bin_to = 1, score_ranges){
    query_start <- start(query_ranges)
    query_stop <- end(query_ranges)
    Main_seq <- seq(from = query_start, to = query_stop)
    Main_pos <- Main_seq - query_start
    end_dist_from_start <- end(score_ranges) - query_start
    start_dist_to_end <- query_stop - start(score_ranges)
    ranges_in_feature <- score_ranges[(end_dist_from_start >= 0 & start_dist_to_end >= 0)]
    Step_size <- width(query_ranges)/bin_to
    a_score <- 0
    if(length(ranges_in_feature) == 0){
        return(a_score)
    }
    score_per_basepair_list <- lapply(seq_along(ranges_in_feature), function(x){
        a_temp_ranges <- ranges_in_feature[x]
        a_start <- start(a_temp_ranges)
        an_end <- end(a_temp_ranges)
        Seq <- seq(from = a_start, to = an_end)
        position_offset_by_start <- Seq - query_start
        if(any(position_offset_by_start < 0)){
            position_offset_by_start <- position_offset_by_start[position_offset_by_start >= 0]
        }
        if(any(position_offset_by_start > max(Main_pos))){
            position_offset_by_start <- position_offset_by_start[position_offset_by_start <= max(Main_pos)]
        }
        data.frame(pos = position_offset_by_start, 
            score = a_temp_ranges$score, 
            frag = x)
    })
    score_per_basepair_df <- do.call(rbind, score_per_basepair_list)
    
    if(any(!(Main_pos %in% score_per_basepair_df$pos))){
        temp_df <- rbind(score_per_basepair_df,
            data.frame(pos = Main_pos[!(Main_pos %in% score_per_basepair_df$pos)],
                score = 0, frag = length(ranges_in_feature)+1))
        score_per_basepair_df <- rbind(score_per_basepair_df, temp_df)
        score_per_basepair_df <- score_per_basepair_df[order(score_per_basepair_df$pos),]
    }
    score_per_basepair_df$pos[score_per_basepair_df$pos == 0] <- score_per_basepair_df$pos[score_per_basepair_df$pos == 0] + 1
    score_per_basepair_df$bins <- ceiling(score_per_basepair_df$pos/Step_size)
    Score_split <- split(score_per_basepair_df$score, score_per_basepair_df$bins)
    frag_split <-  split(score_per_basepair_df$frag, score_per_basepair_df$bins)
    a_score <- vapply(seq_along(Score_split), function(x){
        frag_table <- as.vector(table(frag_split[[x]]))
        scores <- vapply(split(Score_split[[x]], frag_split[[x]]), unique, 1)
        weighted.mean(scores, frag_table)
    }, 1)
    return(a_score)
}

fetch_signal_from_bigwig <- function(query_ranges, bigwig_file, minus_strand = F){
    if(is.null(names(query_ranges))){
        names(query_ranges) <- paste("Seq_id", seq_along(query_ranges), sep = "_")
    }
    scored_ranges <- import.bw(bigwig_file, which = query_ranges)
    if(minus_strand){
        scored_ranges$score <- scored_ranges$score * -1
    }
    message("Select only those SSRs covered by atleast 1 corresponding coverage ranges mid...")
    ol_object <- findOverlaps(query_ranges, scored_ranges)
    scored_ranges_split <- split(scored_ranges[subjectHits(ol_object)], 
        names(query_ranges[queryHits(ol_object)]))
    query_ranges$score <- 0
    values_list <- lapply(names(scored_ranges_split), function(a_ranges_name){
        a_ranges <- query_ranges[a_ranges_name]
        Strand_info <- unique(as.vector(strand(a_ranges)))
        current_score_ranges <- scored_ranges_split[[a_ranges_name]]
        values <- .return_bigwig_signal(query_ranges = a_ranges,
            bin_to = 1, score_ranges = current_score_ranges)
        return(values)
    })
    query_ranges[names(scored_ranges_split)]$score <- do.call(c, values_list)
    ordered_score <- query_ranges$score
    return(ordered_score)
}