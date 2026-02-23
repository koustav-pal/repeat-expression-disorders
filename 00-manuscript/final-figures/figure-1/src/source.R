# Author: Koustav Pal 
# Affiliation: The Francis Crick Institute
# Contact: koustav.pal@crick.ac.uk
# ==========================================================================
# Description
# ==========================================================================
# Source functions for computing overlap weighted basepair coverage from
# provided bw files. 
# The main function is the .return_bigwig_signal function and 
# fetch_bigwig_signal is just a wrapper around this function

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

fetch_bigwig_signal <- function(repeat_ranges, index_start, index_end, bw_signal){
    current_ranges <- repeat_ranges[index_start:index_end]
    repeat_matrix <- do.call(rbind, lapply(seq_along(current_ranges), function(x){
        message(paste("working on repeats", x))
        .return_bigwig_signal(query_ranges = current_ranges[x], bin_to = 10, score_ranges = bw_signal)
    }))

    repeat_start_ranges <- current_ranges
    end(repeat_start_ranges) <- start(repeat_start_ranges)
    start(repeat_start_ranges) <- start(repeat_start_ranges) - 100
    end(repeat_start_ranges) <- end(repeat_start_ranges) - 1

    repeat_start_matrix <- do.call(rbind, lapply(seq_along(repeat_start_ranges), function(x){
        message(paste("working on upstream 100bp", x))
        .return_bigwig_signal(query_ranges = repeat_start_ranges[x], bin_to = 10, score_ranges = bw_signal)
    }))

    repeat_end_ranges <- current_ranges
    start(repeat_end_ranges) <- end(repeat_end_ranges)
    end(repeat_end_ranges) <- end(repeat_end_ranges) + 100
    start(repeat_end_ranges) <- start(repeat_end_ranges) + 1
    repeat_end_matrix <- do.call(rbind, lapply(seq_along(repeat_end_ranges), function(x){
        message(paste("working on downstream 100bp", x))
        .return_bigwig_signal(query_ranges = repeat_end_ranges[x], bin_to = 10, score_ranges = bw_signal)
    }))
    return_matrix <- cbind(repeat_start_matrix, repeat_matrix, repeat_end_matrix)
    return(return_matrix)
}

summarise_for_barplot <- function(events_df, total_length){
    require("dplyr")
    events_df_split <- split(events_df, paste(events_df$study, events_df$contrasts, sep = "_"))
    frequency_df <- do.call(rbind, lapply(events_df_split, function(temp_df){
        events_up <- length(unique(temp_df$feature_ids[temp_df$logfc > 0]))
        events_down <- length(unique(temp_df$feature_ids[temp_df$logfc < 0]))
        return_df <- data.frame(study = unique(temp_df$study),
            contrast = unique(temp_df$contrasts),
            direction = c("up", "down"),
            frequency = c(events_up, events_down))
        return_df$proportion <- return_df$frequency / total_length
        return_df <- return_df %>% 
            arrange(desc(direction)) %>% 
            mutate(prop_lab_y = cumsum(proportion) - (0.5 * proportion),
                freq_lab_y = cumsum(frequency) - (0.5 * frequency)) %>%
            as.data.frame()
        return(return_df)
    }))
    rownames(frequency_df) <- NULL
    return(frequency_df)
}

summarise_for_plot_by_background <- function(events_df, total_length, split_groups){
    require("dplyr")
    groups <- do.call(paste, lapply(split_groups, function(x){
            events_df[,x]  
        }))
    events_df_split <- split(events_df, paste(events_df$background, events_df$study, sep = "_"))
    frequency_df <- do.call(rbind, lapply(events_df_split, function(temp_df){
        num_events <- length(unique(temp_df$feature_ids))
        return_df <- data.frame(study = unique(temp_df$study),
            
            contrast = unique(temp_df$contrasts),
            frequency = num_events)
        return_df$proportion <- return_df$frequency / total_length
        return_df <- return_df %>% 
            arrange(desc(direction)) %>% 
            mutate(prop_lab_y = cumsum(proportion) - (0.5 * proportion),
                freq_lab_y = cumsum(frequency) - (0.5 * frequency)) %>%
            as.data.frame()
        return(return_df)
    }))
    rownames(frequency_df) <- NULL
    return(frequency_df)
}
