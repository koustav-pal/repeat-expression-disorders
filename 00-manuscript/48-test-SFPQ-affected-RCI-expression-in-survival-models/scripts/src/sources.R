

extractC <- function(val) {

  # Extract optimism-corrected Dxy
  dxy <- val["Dxy","index.corrected"]

  # Convert to C-index
  cindex <- (1 + dxy) / 2

  # Extract bootstrap CI bounds for Dxy
  dxy_lower <- val["Dxy","Lower"]
  dxy_upper <- val["Dxy","Upper"]

  # Convert CI bounds to C-index scale
  cindex_lower <- (1 + dxy_lower) / 2
  cindex_upper <- (1 + dxy_upper) / 2

  # Other useful metrics
  slope <- val["Slope","index.corrected"]
  optimism <- val["Dxy","optimism"]

  # Return summary table
  data.frame(
    C = cindex,
    C_lower = cindex_lower,
    C_upper = cindex_upper,
    Dxy = dxy,
    Slope = slope,
    Optimism = optimism
  )
}

extract_cox_report <- function(base_fit,
                               full_fit,
                               term,
                               label = NULL,
                               conf_level = 0.95) {

  stopifnot(inherits(base_fit, "coxph"), inherits(full_fit, "coxph"))
  if (is.null(label)) label <- term

  # ---- Cox summary (term-level) ----
  s_full <- summary(full_fit)

  if (!(term %in% rownames(s_full$coefficients))) {
    stop(sprintf("term '%s' not found in full_fit coefficients.", term))
  }

  # coef table
  coef_tab <- s_full$coefficients
  ci_tab   <- s_full$conf.int

  hr       <- unname(ci_tab[term, "exp(coef)"])
  hr_low   <- unname(ci_tab[term, "lower .95"])
  hr_high  <- unname(ci_tab[term, "upper .95"])
  p_wald   <- unname(coef_tab[term, "Pr(>|z|)"])

  # ---- Model-level concordance (C-index) + approx CI ----
  c_idx <- unname(s_full$concordance[1])
  c_se  <- unname(s_full$concordance[2])

  alpha <- 1 - conf_level
  zcrit <- qnorm(1 - alpha/2)

  c_low <- c_idx - zcrit * c_se
  c_high <- c_idx + zcrit * c_se

  # clamp to [0,1]
  c_low <- max(0, min(1, c_low))
  c_high <- max(0, min(1, c_high))

  # ---- LRT vs baseline (anova) ----
  a <- anova(base_fit, full_fit, test = "LRT")
  chisq <- unname(a$Chisq[2])
  df    <- unname(a$Df[2])
  p_lrt <- unname(a$`Pr(>|Chi|)`[2])

  # ---- Score (logrank) test for full model ----
  # In Cox PH output, the "Score (logrank) test" is available via summary()$sctest
  # sctest = c(test_statistic, df, pvalue)
  if (!is.null(s_full$sctest) && length(s_full$sctest) >= 3) {
    logrank_chisq <- unname(s_full$sctest[1])
    logrank_df    <- unname(s_full$sctest[2])
    logrank_p     <- unname(s_full$sctest[3])
  } else {
    # Fallback: compute it directly from the fitted model
    # (rarely needed, but safer across versions)
    logrank_chisq <- unname(full_fit$score)
    logrank_df    <- unname(length(full_fit$coefficients))
    logrank_p     <- unname(pchisq(logrank_chisq, df = logrank_df, lower.tail = FALSE))
  }

  # ---- Sample size + events ----
  n <- full_fit$n
  events <- full_fit$nevent

  # ---- Log-likelihoods ----
  ll_base <- as.numeric(logLik(base_fit))
  ll_full <- as.numeric(logLik(full_fit))

  data.frame(
    model = label,
    term = term,

    # term effect
    HR = hr,
    HR_lower = hr_low,
    HR_upper = hr_high,
    p_wald = p_wald,

    # model comparison
    logLik_base = ll_base,
    logLik_full = ll_full,
    LRT_ChiSq = chisq,
    LRT_df = df,
    p_LRT = p_lrt,

    # full-model Score (logrank) test
    Score_LogRank_ChiSq = logrank_chisq,
    Score_LogRank_df    = logrank_df,
    Score_LogRank_p     = logrank_p,

    # discrimination
    C_index = c_idx,
    C_index_lower = c_low,
    C_index_upper = c_high,
    C_index_se = c_se,

    # counts
    n = n,
    events = events,

    row.names = NULL,
    check.names = FALSE
  )
}


extract_pairwise_cox_report <- function(base_fit,
                               full_fit,
                               term1,
                               term2,
                               label = NULL,
                               conf_level = 0.95) {

  stopifnot(inherits(base_fit, "coxph"), inherits(full_fit, "coxph"))
  if (is.null(label)) label <- paste(term1, term2, sep = "+")

  # ---- Cox summary (term-level) ----
  s_full <- summary(full_fit)

  if (!(term1 %in% rownames(s_full$coefficients)) | !(term2 %in% rownames(s_full$coefficients))) {
    stop(sprintf("terms '%s' or '%s' not found in full_fit coefficients.", term1, term2))
  }

  # coef table
  coef_tab <- s_full$coefficients
  ci_tab   <- s_full$conf.int

  term1_hr       <- unname(ci_tab[term1, "exp(coef)"])
  term1_hr_low   <- unname(ci_tab[term1, "lower .95"])
  term1_hr_high  <- unname(ci_tab[term1, "upper .95"])
  term1_p_wald   <- unname(coef_tab[term1, "Pr(>|z|)"])

  term2_hr       <- unname(ci_tab[term2, "exp(coef)"])
  term2_hr_low   <- unname(ci_tab[term2, "lower .95"])
  term2_hr_high  <- unname(ci_tab[term2, "upper .95"])
  term2_p_wald   <- unname(coef_tab[term2, "Pr(>|z|)"])

  # ---- Model-level concordance (C-index) + approx CI ----
  c_idx <- unname(s_full$concordance[1])
  c_se  <- unname(s_full$concordance[2])

  alpha <- 1 - conf_level
  zcrit <- qnorm(1 - alpha/2)

  c_low <- c_idx - zcrit * c_se
  c_high <- c_idx + zcrit * c_se

  # clamp to [0,1]
  c_low <- max(0, min(1, c_low))
  c_high <- max(0, min(1, c_high))

  # ---- LRT vs baseline (anova) ----
  a <- anova(base_fit, full_fit, test = "LRT")
  chisq <- unname(a$Chisq[2])
  df    <- unname(a$Df[2])
  p_lrt <- unname(a$`Pr(>|Chi|)`[2])

  # ---- Score (logrank) test for full model ----
  # In Cox PH output, the "Score (logrank) test" is available via summary()$sctest
  # sctest = c(test_statistic, df, pvalue)
  if (!is.null(s_full$sctest) && length(s_full$sctest) >= 3) {
    logrank_chisq <- unname(s_full$sctest[1])
    logrank_df    <- unname(s_full$sctest[2])
    logrank_p     <- unname(s_full$sctest[3])
  } else {
    # Fallback: compute it directly from the fitted model
    # (rarely needed, but safer across versions)
    logrank_chisq <- unname(full_fit$score)
    logrank_df    <- unname(length(full_fit$coefficients))
    logrank_p     <- unname(pchisq(logrank_chisq, df = logrank_df, lower.tail = FALSE))
  }

  # ---- Sample size + events ----
  n <- full_fit$n
  events <- full_fit$nevent

  # ---- Log-likelihoods ----
  ll_base <- as.numeric(logLik(base_fit))
  ll_full <- as.numeric(logLik(full_fit))

  data.frame(
    model = label,
    term1 = term1,
    term2 = term2,

    # term1 effect
    term1_HR = term1_hr,
    term1_HR_lower = term1_hr_low,
    term1_HR_upper = term1_hr_high,
    term1_p_wald = term1_p_wald,

    # term2 effect
    term2_HR = term2_hr,
    term2_HR_lower = term2_hr_low,
    term2_HR_upper = term2_hr_high,
    term2_p_wald = term2_p_wald,

    # model comparison
    logLik_base = ll_base,
    logLik_full = ll_full,
    LRT_ChiSq = chisq,
    LRT_df = df,
    p_LRT = p_lrt,

    # full-model Score (logrank) test
    Score_LogRank_ChiSq = logrank_chisq,
    Score_LogRank_df    = logrank_df,
    Score_LogRank_p     = logrank_p,

    # discrimination
    C_index = c_idx,
    C_index_lower = c_low,
    C_index_upper = c_high,
    C_index_se = c_se,

    # counts
    n = n,
    events = events,

    row.names = NULL,
    check.names = FALSE
  )
}

