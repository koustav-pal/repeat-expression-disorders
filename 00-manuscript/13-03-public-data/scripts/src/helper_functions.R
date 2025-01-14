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