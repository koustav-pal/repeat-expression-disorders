data_dir <- normalizePath("/scratch/projects/CFP03/CFP03-SF-111/koustav.pal/")
project_dir <- file.path(data_dir, "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript")
setwd(analysis_dir)
save_dir <- file.path(analysis_dir, "48-test-SFPQ-affected-RCI-expression-in-survival-models", "rdata_files")
input_dir <- file.path(analysis_dir,"48-test-SFPQ-affected-RCI-expression-in-survival-models", "input_files")
plot_dir <- file.path(analysis_dir, "48-test-SFPQ-affected-RCI-expression-in-survival-models", "plots")
output_dir <- file.path(analysis_dir, "48-test-SFPQ-affected-RCI-expression-in-survival-models", "output_files")
gtac_predicted_scores_df <- readRDS(file.path(save_dir, "06-04-save-gtac-predicted-score.rds"))
grima_surv_df <- readRDS(file.path(save_dir, "06-05-drosha-cpfs6-cpeb4-scores-in-grima-et-al.rds"))
# ==========================================================================
# Packages and Data
# ==========================================================================
# Change to normal R conda environment
# source("../../src/plot-utils.R")
# source("../../src/basepair_coverage_functions.R")
source("src/colours.R")
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
make_frozen_cuts <- function(dfA, score_vars,
                             method = c("median", "quantile"),
                             probs = c(0.5),   # e.g. c(0.5) median; c(1/3,2/3) tertiles; c(0.25,0.5,0.75) quartiles
                             na.rm = TRUE) {

  method <- match.arg(method)

  if (!all(score_vars %in% names(dfA))) {
    missing <- score_vars[!score_vars %in% names(dfA)]
    stop("These score_vars are not in dfA: ", paste(missing, collapse = ", "))
  }

  # compute cut(s) per score in cohort A
  cuts <- setNames(vector("list", length(score_vars)), score_vars)

  for (v in score_vars) {
    x <- dfA[[v]]
    if (method == "median") {
      cuts[[v]] <- stats::median(x, na.rm = na.rm)
    } else {
      qs <- stats::quantile(x, probs = probs, na.rm = na.rm, names = FALSE, type = 7)
      qs <- as.numeric(qs)
      if (any(!is.finite(qs))) stop("Non-finite cutpoints for ", v)
      cuts[[v]] <- qs
    }
  }

  # function that applies the learned cut(s) to a new cohort
  apply_fun <- function(dfB, suffix = "_group", include_score_z = FALSE,
                        zscore_using_A = FALSE) {
    out <- dfB

    for (v in score_vars) {
      if (!(v %in% names(out))) stop("Variable ", v, " not found in dfB")

      xB <- out[[v]]

      # Optional: frozen z-scoring using cohort A mean/sd (recommended if scores are on comparable scale)
      if (zscore_using_A) {
        muA <- mean(dfA[[v]], na.rm = TRUE)
        sdA <- stats::sd(dfA[[v]], na.rm = TRUE)
        if (!is.finite(sdA) || sdA == 0) stop("Cannot z-score: sd in A is 0/NA for ", v)
        xB_use <- (xB - muA) / sdA
        if (include_score_z) out[[paste0(v, "_z")]] <- xB_use
      } else {
        xB_use <- xB
      }

      cA <- cuts[[v]]

      # median split (single cut)
      if (length(cA) == 1) {
        out[[paste0(v, suffix)]] <- ifelse(xB_use >= cA, "High", "Low")
        out[[paste0(v, suffix)]] <- factor(out[[paste0(v, suffix)]], levels = c("Low", "High"))
      } else {
        # multiple cutpoints -> quantile bins
        # Labels: Q1..Q(k+1)
        breaks <- c(-Inf, cA, Inf)
        labs <- paste0("Q", seq_len(length(breaks) - 1))
        out[[paste0(v, suffix)]] <- cut(xB_use, breaks = breaks, labels = labs, include.lowest = TRUE)
      }
    }

    out
  }

  list(
    cuts = cuts,
    method = method,
    probs = probs,
    score_vars = score_vars,
    apply = apply_fun
  )
}
# Helper: given a survdiff object, compute "relative hazard" like you were using before
relative_hazard_from_survdiff <- function(sd) {
  # sd$obs and sd$exp are in the order of factor levels
  # With levels c("Low","High"), this returns (High/Low) hazard proxy
  if (length(sd$obs) != 2) stop("relative_hazard_from_survdiff expects exactly 2 groups.")
  (sd$obs[2] / sd$exp[2]) / (sd$obs[1] / sd$exp[1])
}
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
  # med <- median(df[[score_var]], na.rm = TRUE)
  # df$risk_group <- ifelse(df[[score_var]] >= med, "High", "Low")
  # df$risk_group <- factor(df$risk_group, levels = c("Low", "High"))

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
    stringsAsFactors = FALSE)

  list(plot = g, stats = stats, fit_km = fit_km, fit_cox = fit_cox, survdiff = sd)
}

frozen <- make_frozen_cuts(grima_surv_df, method = "median",
  probs =  c(0.5),
  score_vars = c("drosha_score", "cpeb4_score", "cpsf6_score"))

stratified_gtac_df <- frozen$apply(gtac_predicted_scores_df)



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


scores <- c("drosha_score_group", "cpeb4_score_group", "cpsf6_score_group")
km_results <- lapply(scores[1], function(x) {
  stratified_gtac_df$risk_group <- stratified_gtac_df[,x]
  if(length(unique(stratified_gtac_df$risk_group)) == 1){
    return(NULL)
  }
	km_and_stats(stratified_gtac_df, score_var = x)
})
names(km_results) <- scores[1]

table(stratified_gtac_df$drosha_score_group)


# Extract plots
km_plots <- lapply(km_results, `[[`, "plot")

# Extract a single reporting table
km_stats_table <- do.call(rbind, lapply(km_results, `[[`, "stats"))
km_stats_table

all_plots <- wrap_plots(km_plots[c(1)], nrow = 1, ncol = 1)

ggsave(all_plots, file = file.path(plot_dir, "06-external-validation-DROSHA-CPSF6-CPEB4-ridge-regression-based-survival-curves.pdf"), width = 5, height = 5)

forest_plots <- lapply(scores[1], function(x){
		stratified_gtac_df$score <- stratified_gtac_df[,x]
		model <- coxph(Surv(time, status) ~ age_at_onset + sex + score, data = stratified_gtac_df, x=TRUE)
    ggforest(model, data = stratified_gtac_df) + ggtitle(x)
})

all_forest_plots <- wrap_plots(forest_plots, nrow = 1, ncol = 1)
ggsave(all_forest_plots, file = file.path(plot_dir, "07-external-validation-DROSHA-CPSF6-CPEB4-ridge-regression-based-forest-plots.pdf"), width = 5, height = 5)
