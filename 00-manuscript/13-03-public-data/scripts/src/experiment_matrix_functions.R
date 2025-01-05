perform_diff_expr <- function(counts_matrix, sample_df, FUN, sample_name, min_count = 1, 
    study_specific_design = T, save_path){
    experiment_df_list <- FUN(sample_df, study_specific_design = study_specific_design)
    experiment_matrix <- experiment_df_list$experiment
    full_model <- experiment_df_list$design
    counts_matrix <- counts_matrix[rowMeans(counts_matrix) > min_count,]
    deseq_object <- DESeqDataSetFromMatrix(countData = round(counts_matrix[,experiment_matrix$samplename]), 
        colData = experiment_matrix, 
        design = full_model)
    differential_object <- DESeq(deseq_object)
    saveRDS(object = differential_object, file = save_path)
    return(data.frame(sample_name = sample_name,
            save_path = save_path))
}

return_function_list <- function(){
    c9orf72_casteli_2021 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["c9orf72-casteli-2021"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "c9orf72-casteli-2021",
            condition = sample_split[,2], 
            knockdown = sample_split[,3], 
            fraction = sample_split[,4],
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "c9orf72"))
        experiment_df$knockdown <- factor(experiment_df$knockdown, levels = c("scramble", "sfrsf1"))
        experiment_df$fraction <- factor(experiment_df$fraction, levels = c("who", "cyt"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, 
            experiment_df$knockdown, experiment_df$fraction, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition * knockdown * fraction    
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[experiment_df$knockdown == "scramble" & 
                experiment_df$fraction == "scramble",
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + condition + fraction
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    c9orf72_catanese_2021 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["c9orf72-catanese-2021"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "c9orf72-catanese-2021",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "iso", "c9orf72", "fus"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    c9orf72_dafinca_2021 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["c9orf72-dafinca-2020"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "c9orf72-dafinca-2020",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "iso", "c9orf72"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    c9orf72_sareen_2013 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["c9orf72-sareen-2013"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "c9orf72-sareen-2013",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "c9orf72"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    c9orf72_sommer_2022 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["c9orf72-sommer-2022"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "c9orf72-sommer-2022",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "c9orf72"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    c9orf72_sterneckert_2020 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["c9orf72-sterneckert-2020"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "c9orf72-sterneckert-2020",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "c9orf72"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    fus_desantis_2017 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["fus-desantis-2017"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "fus-desantis-2017",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "fus"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    fus_hawkins_2022 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["fus-hawkins-2022"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "fus-hawkins-2022",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "fus"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    fus_kapeli_2016 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["fus-kapeli-2016"]]
        current_df <- current_df[1:6,]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)[1:6,]
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "fus-kapeli-2016",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "fus"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    hnrnpa2b1_markmiller_2021 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["hnrnpa2b1-markmiller-2021"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "hnrnpa2b1-markmiller-2021",
            condition = sampleinfo_df$condition[match(sample_split[,2], sampleinfo_df$cellline)], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "als"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    sod1_bhinge_2017 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["sod1-bhinge-2017"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = "sod1-bhinge-2017",
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("iso", "sod1"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    sod1_moccia_2014 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["sod1-moccia-2014"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("iso", "sod1"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    sod1_tardbp_dash_2022 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["sod1-tardbp-dash-2022"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "sod1", "tardbp"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    sod1_wang_2017 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["sod1-wang-2017"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("iso", "sod1"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(
                experiment = experiment_df,
                design = simple_model
            ))
        }
    }
    tardbp_dafinca_2020 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["tardbp-dafinca-2020"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "tardbp"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    tardbp_kd_brown_2022 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["tardbp-kd-brown-2022"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "tdp43kd"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    tardbp_kd_kapeli_2016 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["tardbp-kd-kapeli-2016"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "tdp43kd"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model
            ))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    tardbp_kd_klim_2022 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["tardbp-kd-klim-2022"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "tdp43kd"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    tardbp_melamed_2019 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["tardbp-melamed-2019"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "tardbp"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    tardbp_smith_2021 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["tardbp-smith-2021"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,2], 
            fraction = "who", 
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "tardbp"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition
            return(list(
                experiment = experiment_df,
                design = full_model))
        }else{
            experiment_df <- experiment_df[,
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    vcp_luisier_2017 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["vcp-luisier-2018"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,5], 
            fraction = "who", 
            timepoint = sample_split[,4],
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded")))
        experiment_df$condition[experiment_df$condition == "Control"] <- "ctrl"
        experiment_df$condition[experiment_df$condition == "VCP"] <- "vcp"
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "vcp"))
        experiment_df$timepoint <- factor(experiment_df$timepoint, levels = c("iPSC", "d7", "d14", "d21", "d35"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, experiment_df$timepoint, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition * timepoint
            return(list(
                experiment = experiment_df,
                design = full_model))
        }else{
            experiment_df <- experiment_df[experiment_df$timepoint == "d35",
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    vcp_ziff_2023 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["vcp-ziff-2023"]]
        sampleinfo_df <- read.csv(unique(current_df$samplesheet_path))
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        sample_split_4_split <- str_split(sample_split[,4], pattern = "\\.", simplify = T)
        conditions <- sampleinfo_df[match(sample_split_4_split[,1], sampleinfo_df[,"cell.line.ch1"]) , "genotype.ch1"]
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            genotype = sample_split_4_split[,1],
            condition = conditions, 
            fraction = sample_split[,2], 
            treatment = sample_split[,3],
            instrument = unique(sampleinfo_df$instrument),
            library_layout = ifelse(current_df$isPairedEnd, "paired-end", "single-end"),
            directionality = ifelse(current_df$strand_specificity == 2, "reverse", 
                ifelse(current_df$strand_specificity == 1, "forward", "unstranded"))) 
        experiment_df$condition <- factor(experiment_df$condition, 
            levels = c("wildtype", "vcp knock-in", "vcp mutant", "tardbp mutant"))
        experiment_df$fraction <- factor(experiment_df$fraction, levels = c("nucleus", "cytoplasm"))
        experiment_df$genotype <- factor(experiment_df$genotype)
        experiment_df$treatment <- factor(experiment_df$treatment, levels = c("untreated", "ml240"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, experiment_df$fraction, experiment_df$treatment, sep = "_")
        experiment_df$unique_sample_name <- str_replace_all(string = experiment_df$unique_sample_name, pattern = " ", replacement = "_")
        if(study_specific_design){
            full_model <- ~ condition * fraction * treatment
            return(list(
                experiment = experiment_df,
                design = full_model))
        }else{
            experiment_df <- experiment_df[experiment_df$treatment == "untreated",
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + instrument + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    vcp_tyzack_2021 <- function(sample_df, study_specific_design = T){
        sample_df_split <- split(sample_df, sample_df$sample_name)
        current_df <- sample_df_split[["vcp-tyzack-2021"]]
        sample_split <- str_split(current_df$column_id, pattern = "_", simplify = T)
        experiment_df <- data.frame(samplename = current_df$column_id,
            study = sample_split[,1],
            condition = sample_split[,5], 
            fraction = sample_split[,3], 
            timepoint = sample_split[,4],
            instrument = "Illumina NovaSeq 6000",
            library_layout = "paired-end",
            directionality = "reverse")
        experiment_df$condition <- factor(experiment_df$condition, levels = c("ctrl", "vcp"))
        experiment_df$fraction <- factor(experiment_df$fraction, levels = c("nuc", "cyt"))
        experiment_df$timepoint <- factor(experiment_df$timepoint, levels = c("d0", "d3", "d7", "d14", "d22", "d35"))
        experiment_df$library_layout <- factor(experiment_df$library_layout)
        experiment_df$directionality <- factor(experiment_df$directionality)
        experiment_df$unique_sample_name <- paste(experiment_df$study, experiment_df$condition, experiment_df$timepoint, sep = "_")
        if(study_specific_design){
            full_model <- ~ condition * fraction * timepoint
            return(list(
                experiment = experiment_df,
                design = full_model))
        }else{
            experiment_df <- experiment_df[experiment_df$timepoint == "d35",
                c("samplename", "study", "instrument", "library_layout", "condition", "fraction")]
            simple_model <- ~ study + library_layout + fraction + condition
            return(list(experiment = experiment_df,
                design = simple_model))
        }
    }
    all_functions <- list(
        "c9orf72-dafinca-2020" = c9orf72_dafinca_2021,
        "c9orf72-catanese-2021" = c9orf72_catanese_2021,
        "c9orf72-sareen-2013" = c9orf72_sareen_2013,
        "c9orf72-sommer-2022" = c9orf72_sommer_2022,
        "c9orf72-sterneckert-2020" = c9orf72_sterneckert_2020,
        "c9orf72-casteli-2021" = c9orf72_casteli_2021,
        "fus-desantis-2017" = fus_desantis_2017,
        "fus-hawkins-2022" = fus_hawkins_2022,
        "fus-kapeli-2016" = fus_kapeli_2016,
        "hnrnpa2b1-markmiller-2021" = hnrnpa2b1_markmiller_2021,
        "sod1-bhinge-2017" = sod1_bhinge_2017,
        "sod1-moccia-2014" = sod1_moccia_2014,
        "sod1-tardbp-dash-2022" = sod1_tardbp_dash_2022,
        "sod1-wang-2017" = sod1_wang_2017,
        "tardbp-dafinca-2020" = tardbp_dafinca_2020,
        "tardbp-kd-brown-2022" = tardbp_kd_brown_2022,
        "tardbp-kd-klim-2022" = tardbp_kd_klim_2022,
        "tardbp-melamed-2019" = tardbp_melamed_2019,
        "tardbp-smith-2021" = tardbp_smith_2021,
        "tardbp-kd-kapeli-2016" = tardbp_kd_kapeli_2016,
        "vcp-luisier-2018" = vcp_luisier_2017,
        "vcp-ziff-2023" = vcp_ziff_2023,
        "vcp-tyzack-2021" = vcp_tyzack_2021)
    return(all_functions)
}
