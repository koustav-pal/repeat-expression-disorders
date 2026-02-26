library("here")
message("Project root:", here(), "\n")
project_dir <- here()
analysis_dir <- here("final-figures", "figure-5")
save_dir <- file.path(analysis_dir, "rdata_files")
input_dir <- file.path(analysis_dir,"input_files")
plot_dir <- file.path(analysis_dir, "plots")
output_dir <- file.path(analysis_dir, "output_files")
model_validation_list <- readRDS(file.path(save_dir, "06-03-save-all-ridge-coefficient-models-from-grima-et-al.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
# source("../../src/plot-utils.R")
# source("../../src/basepair_coverage_functions.R")
# source("src/colours.R")
source("48-test-SFPQ-affected-RCI-expression-in-survival-models/scripts/src/sources.R")
require("HiCBricks")
require("GenomicRanges")
require("stringr")
require("ggplot2")
library(rlang)
library(riskRegression)
library(survival)
library(survminer)
library(transformr)
library(plotly)
library(reshape2)
library(knitr)
library(rmarkdown)
require("RColorBrewer")
require("randomForestSRC")
require("parallel")
require("matrixStats")
require("survival")
require("patchwork")
require("ggplot2")
library(broom)
library(ComplexUpset)
library(ggplot2)
library(ggrepel)
require(rms)
require(patchwork)
# ==========================================================================
# Function
# ==========================================================================
library(survival)
library(broom)
library(survminer)
library(ggplot2)

# Helper: given a survdiff object, compute "relative hazard" like you were using before
relative_hazard_from_survdiff <- function(sd) {
  # sd$obs and sd$exp are in the order of factor levels
  # With levels c("Low","High"), this returns (High/Low) hazard proxy
  if (length(sd$obs) != 2) stop("relative_hazard_from_survdiff expects exactly 2 groups.")
  (sd$obs[2] / sd$exp[2]) / (sd$obs[1] / sd$exp[1])
}

# Main: build KM plot + extract HR, log-rank stat, and relative hazard
km_and_stats <- function(surv_df, score_var, time_var = "time", status_var = "status",
                         palette = c("#FFCCC9", "#FF5666"),
                         conf_level = 0.95) {

	current_theme <- theme(
	    legend.position = "bottom",
	    text = element_text(size = 12),
	    axis.text.x = element_text(size = 10),
	    axis.text.y = element_text(size = 12),
	    legend.text = element_text(size = 10),
	    legend.title = element_text(size = 10),
	    plot.background = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.grid.major = element_line(colour = "azure3",
	    linetype = "dashed",
	    size = 0.2))

  # Make a local copy (avoid overwriting outside)
  df <- surv_df

  # risk group by median split
  med <- median(df[[score_var]], na.rm = TRUE)
  df$risk_group <- ifelse(df[[score_var]] >= med, "High", "Low")
  df$risk_group <- factor(df$risk_group, levels = c("Low", "High"))

  # KM fit
  fit_km <- survfit(Surv(time,status) ~ risk_group, data = df)
 
  # Log-rank test
  sd <- survdiff(Surv(time,status) ~ risk_group, data = df)

  logrank_chisq <- unname(sd$chisq)
  logrank_df <- length(sd$n) - 1
  logrank_p <- pchisq(logrank_chisq, df = logrank_df, lower.tail = FALSE)

  # "Relative hazard" proxy from survdiff (High vs Low)
  rel_haz <- relative_hazard_from_survdiff(sd)

  # Cox PH for HR (High vs Low)
  fit_cox <- coxph(Surv(time,status) ~ risk_group, data = df)
  s_cox <- summary(fit_cox)

  hr <- unname(s_cox$conf.int["risk_groupHigh", "exp(coef)"])
  hr_low <- unname(s_cox$conf.int["risk_groupHigh", "lower .95"])
  hr_high <- unname(s_cox$conf.int["risk_groupHigh", "upper .95"])
  p_wald <- unname(s_cox$coefficients["risk_groupHigh", "Pr(>|z|)"])

  # Optional: concordance of this 1-variable Cox model
  c_idx <- unname(s_cox$concordance[1])
  c_se  <- unname(s_cox$concordance[2])
  zcrit <- qnorm(1 - (1 - conf_level)/2)
  c_low <- max(0, min(1, c_idx - zcrit * c_se))
  c_high <- max(0, min(1, c_idx + zcrit * c_se))

  # Plot (no table returned, just the ggplot like your code)
  g <- ggsurvplot(
    fit_km, data = df,
    palette = palette,
    pval = FALSE,           # we'll annotate ourselves so it's consistent with extracted stats
    conf.int = TRUE,
    risk.table = TRUE, risk.table.height = 0.22,
    xlab = "Time", ylab = "Survival probability",
    legend.title = "Composite score",
    ggtheme = theme_bw() + current_theme
  )$plot +
    ggtitle(score_var) +
    annotate(
      "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
      label = sprintf(
        "HR (High vs Low) = %.2f [%.2f, %.2f]\nLog-rank χ² = %.2f (p = %.3g)\nRel. hazard (O/E) = %.2f\nC = %.3f [%.3f, %.3f]",
        hr, hr_low, hr_high,
        logrank_chisq, logrank_p,
        rel_haz,
        c_idx, c_low, c_high
      ),
      size = 3.3
    )

  stats <- data.frame(
    score = score_var,
    n = nrow(df),
    events = sum(df[[status_var]] == 1, na.rm = TRUE),

    HR_high_vs_low = hr,
    HR_low95 = hr_low,
    HR_high95 = hr_high,
    p_wald = p_wald,

    logrank_chisq = logrank_chisq,
    logrank_df = logrank_df,
    logrank_p = logrank_p,

    relative_hazard_OE_high_vs_low = rel_haz,

    cindex = c_idx,
    cindex_low = c_low,
    cindex_high = c_high,

    stringsAsFactors = FALSE
  )

  list(plot = g, stats = stats, fit_km = fit_km, fit_cox = fit_cox, survdiff = sd)
}


surv_df <- model_validation_list$surv_df
scores <- model_validation_list$all_rbp_scores
selected_rbp_enriched_rcis_df = model_validation_list$selected_rbp_enriched_rcis
cox_result_df <- model_validation_list$coxph_results


surv_df$drosha_score <- scores[["DROSHA hexanucleotide"]]
surv_df$cpeb4_score <- scores[["CPEB4 trinucleotide"]]
surv_df$cpsf6_score <- scores[["CPSF6 trinucleotide"]]

surv_df$drosha_score_scaled <- scale(surv_df$drosha_score)
surv_df$cpeb4_score_scaled <- scale(surv_df$cpeb4_score)
surv_df$cpsf6_score_scaled <- scale(surv_df$cpsf6_score)


current_theme <- theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))


scores <- c("drosha_score_scaled", "cpeb4_score_scaled", "cpsf6_score_scaled")

km_results <- lapply(scores, function(x) {
	km_and_stats(surv_df, x, time_var = "time")
})
names(km_results) <- scores

# Extract plots
km_plots <- lapply(km_results, `[[`, "plot")

# Extract a single reporting table
km_stats_table <- do.call(rbind, lapply(km_results, `[[`, "stats"))
km_stats_table

all_plots <- wrap_plots(km_plots, nrow = 1, ncol = 3)

ggsave(all_plots, file = file.path(plot_dir, "04-DROSHA-CPSF6-CPEB4-ridge-regression-based-survival-curves.pdf"), width = 10, height = 5)

forest_plots <- lapply(scores, function(x){
		surv_df$score <- surv_df[,x]
		model <- coxph(Surv(time, status) ~ age_at_onset + Sex + score, data = surv_df, x=TRUE)
    ggforest(model, data = surv_df) + ggtitle(x)
})

all_forest_plots <- wrap_plots(forest_plots, nrow = 1, ncol = 3)
ggsave(all_forest_plots, file = file.path(plot_dir, "05-DROSHA-CPSF6-CPEB4-ridge-regression-based-forest-plots.pdf"), width = 10, height = 5)

saveRDS(object = surv_df, file = file.path(save_dir, "06-05-drosha-cpfs6-cpeb4-scores-in-grima-et-al.rds"))
