# plot_diagnostic_recommendations.R
# Failure-mode-specific interpretation for the diagnostic-report interpretation
# page (see .probit_interpretation_text() in plot_probit_diagnostics.R).
#
# The problem this file fixes: a single generic line -- "consider increasing
# adapt_delta or simplifying the model" -- fired for BOTH divergence failures
# and energy failures alike. adapt_delta specifically resolves step-size-
# induced divergences; recommending it when n_divergent == 0 (an energy or
# Omega-geometry failure instead) is a category error, not just imprecise
# wording. This file replaces that one-size-fits-all line with a ranked
# primary/secondary/not-indicated classification driven by which canonical
# diagnostic fields actually failed, plus a compact diagnostic table.
#
# Nothing here recomputes a diagnostic value -- every number it reads
# (n_divergent, ebfmi_min, max_rhat_structural, fit$diagnostics_detail$grouped,
# Omega degeneracy stats) is already computed once elsewhere in the package.

#' Central thresholds for diagnostic-report interpretation
#'
#' Single place these thresholds are defined, so the interpretation logic and
#' the report table stay in sync instead of drifting independently.
#'
#' @param ebfmi_min Warn below this minimum chain E-BFMI. Default \code{0.3}
#'   (Stan's own guidance threshold).
#' @param rhat_max Warn above this Rhat. Default \code{1.01}.
#' @param rhat_severe Structural Rhat above this is treated as severe rather
#'   than borderline. Default \code{1.05}.
#' @param treedepth_max Treedepth-saturation count above this triggers the
#'   treedepth failure mode. Default \code{0L} (any saturation).
#' @param divergent_max Divergent-transition count above this triggers the
#'   divergence failure mode. Default \code{0L} (any divergence).
#' @param omega_small_eigen_threshold Per-draw smallest-eigenvalue-of-Omega
#'   value below which a draw is counted as near-singular. Default \code{0.05}.
#' @param omega_small_eigen_fraction_warning Fraction of posterior draws
#'   flagged near-singular (via \code{omega_small_eigen_threshold}) above
#'   which Omega geometry is reported as a primary concern. Default \code{0.50}.
#' @param omega_condition_number_warning Median Omega condition number above
#'   which it is called out explicitly. Default \code{50}.
#' @return Named list of thresholds.
#' @keywords internal
.probit_diagnostic_thresholds <- function(ebfmi_min = 0.3,
                                          rhat_max = 1.01,
                                          rhat_severe = 1.05,
                                          treedepth_max = 0L,
                                          divergent_max = 0L,
                                          omega_small_eigen_threshold = 0.05,
                                          omega_small_eigen_fraction_warning = 0.50,
                                          omega_condition_number_warning = 50) {
  list(
    ebfmi_min = ebfmi_min,
    rhat_max = rhat_max,
    rhat_severe = rhat_severe,
    treedepth_max = treedepth_max,
    divergent_max = divergent_max,
    omega_small_eigen_threshold = omega_small_eigen_threshold,
    omega_small_eigen_fraction_warning = omega_small_eigen_fraction_warning,
    omega_condition_number_warning = omega_condition_number_warning
  )
}

.probit_num1 <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  v <- suppressWarnings(as.numeric(x[[1L]]))
  if (length(v) == 0L) NA_real_ else v
}

.probit_group_row <- function(grouped_diag, group_name) {
  if (is.null(grouped_diag) || nrow(grouped_diag) == 0L || !"parameter_group" %in% names(grouped_diag)) {
    return(NULL)
  }
  row <- grouped_diag[grouped_diag$parameter_group == group_name, , drop = FALSE]
  if (nrow(row) == 0L) NULL else row
}

#' Classify why a fit's sampler diagnostics look the way they do
#'
#' Produces a ranked primary/secondary/not-indicated recommendation set
#' instead of one generic line, so a reviewer is told the SPECIFIC failure
#' mode (divergence vs. energy vs. Omega geometry vs. fixed-effect sparsity vs.
#' treedepth) rather than a boilerplate "increase adapt_delta or simplify the
#' model" that is only actually appropriate for divergences.
#'
#' @param diag \code{fit$diagnostics} (canonical scalar diagnostic fields).
#' @param grouped_diag \code{fit$diagnostics_detail$grouped} -- canonical
#'   per-parameter_group Rhat/ESS table (rows: beta, Omega, lp__, z_free, any
#'   declared random-effect blocks). May be \code{NULL}.
#' @param degeneracy_stats Result of \code{.omega_degeneracy_stats()}, or
#'   \code{NULL} for identity-residual fits.
#' @param residual_structure \code{"identity"} or \code{"correlated"}.
#' @param beta_summary Per-covariate-family beta Rhat/ESS breakdown from
#'   \code{plot_probit_beta_family_diagnostics()$summary}, or \code{NULL}.
#' @param thresholds List from \code{.probit_diagnostic_thresholds()}.
#' @return Named list: \code{primary}, \code{secondary}, \code{not_indicated}
#'   (character vectors, possibly empty) and \code{bucket} (single string
#'   identifying which top-level failure mode was selected as primary, for
#'   testing).
#' @keywords internal
.probit_classify_fitting_issues <- function(diag, grouped_diag, degeneracy_stats,
                                            residual_structure, beta_summary = NULL,
                                            thresholds = .probit_diagnostic_thresholds()) {
  n_divergent <- .probit_num1(diag$n_divergent)
  n_treedepth_sat <- .probit_num1(diag$n_treedepth_sat)
  ebfmi_min <- .probit_num1(diag$ebfmi_min)

  has_divergence <- !is.na(n_divergent) && n_divergent > thresholds$divergent_max
  has_treedepth  <- !is.na(n_treedepth_sat) && n_treedepth_sat > thresholds$treedepth_max
  has_energy     <- !is.na(ebfmi_min) && ebfmi_min < thresholds$ebfmi_min

  omega_row <- .probit_group_row(grouped_diag, "Omega")
  beta_row  <- .probit_group_row(grouped_diag, "beta")

  omega_rhat_pct <- if (!is.null(omega_row)) .probit_num1(omega_row$pct_rhat_above_1_01) else NA_real_
  beta_rhat_pct  <- if (!is.null(beta_row))  .probit_num1(beta_row$pct_rhat_above_1_01)  else NA_real_

  omega_rhat_bad <- !is.na(omega_rhat_pct) && omega_rhat_pct > 0
  beta_rhat_bad  <- !is.na(beta_rhat_pct)  && beta_rhat_pct  > 0

  omega_near_singular <- identical(residual_structure, "correlated") && !is.null(degeneracy_stats) &&
    !is.na(degeneracy_stats$pct_near_degenerate) &&
    (degeneracy_stats$pct_near_degenerate / 100) > thresholds$omega_small_eigen_fraction_warning

  primary <- character(0)
  secondary <- character(0)
  not_indicated <- character(0)
  bucket <- "clean"

  # -- adapt_delta / max_treedepth: state explicitly when NOT indicated, so a
  # reader never has to infer their absence from silence. --
  if (!has_divergence) {
    not_indicated <- c(not_indicated,
      "Increasing adapt_delta -- not indicated (zero divergent transitions; adapt_delta specifically targets step-size-induced divergences).")
  }
  if (!has_treedepth) {
    not_indicated <- c(not_indicated,
      "Increasing max_treedepth -- not indicated (no treedepth saturations).")
  }

  # -- A. Divergence failure: the ONLY case where adapt_delta is a primary remedy --
  if (has_divergence) {
    bucket <- "divergence"
    primary <- c(primary, sprintf(
      "Divergent transitions were detected (n = %.0f). Increasing adapt_delta may help, but persistent divergences after that indicate the model geometry or parameterization should be revised, not just re-tuned.",
      n_divergent
    ))
  }

  # -- B/E/F. Energy failure and/or Omega near-singularity/Rhat failure: these
  # are reported as ONE combined geometry concern when both are present,
  # since near-singular Omega is a plausible direct cause of poor energy
  # exploration in a correlated-residual fit, not an independent finding. --
  geometry_issue <- has_energy || omega_near_singular || (identical(residual_structure, "correlated") && omega_rhat_bad)
  if (geometry_issue) {
    parts <- character(0)
    if (has_energy) {
      parts <- c(parts, sprintf("poor posterior energy exploration (min E-BFMI = %.3f, below the %.2f threshold)", ebfmi_min, thresholds$ebfmi_min))
    }
    if (identical(residual_structure, "correlated")) {
      omega_parts <- character(0)
      if (omega_near_singular) {
        omega_parts <- c(omega_parts, sprintf(
          "%.1f%% of posterior draws have Omega's smallest eigenvalue below %.2f (near-singular)",
          degeneracy_stats$pct_near_degenerate, thresholds$omega_small_eigen_threshold
        ))
      }
      if (!is.null(degeneracy_stats) && !is.na(degeneracy_stats$med_condition_number) &&
          degeneracy_stats$med_condition_number > thresholds$omega_condition_number_warning) {
        omega_parts <- c(omega_parts, sprintf("median condition number %.1f", degeneracy_stats$med_condition_number))
      }
      if (omega_rhat_bad) {
        omega_parts <- c(omega_parts, sprintf("%.1f%% of Omega/L_Omega parameters exceed Rhat %.2f", omega_rhat_pct, thresholds$rhat_max))
      }
      if (length(omega_parts) > 0L) {
        parts <- c(parts, paste0("near-singular residual correlation geometry (", paste(omega_parts, collapse = "; "), ")"))
      }
    }
    if (length(parts) > 0L) {
      combined <- sprintf("Primary fitting concern: %s.", paste(parts, collapse = " and "))
      if (identical(residual_structure, "correlated") && !is.na(beta_rhat_pct)) {
        combined <- paste(combined, sprintf(
          "%s zero divergences and no treedepth saturation, increasing adapt_delta is unlikely to address this.",
          if (has_divergence) "Despite" else "With"
        ))
      }
      if (identical(bucket, "clean")) bucket <- "energy_or_omega_geometry"
      target <- if (identical(bucket, "divergence")) "secondary" else "primary"
      if (identical(target, "primary")) primary <- c(primary, combined) else secondary <- c(secondary, combined)
    }
  }

  # -- C. Treedepth, if not already the headline divergence/geometry issue --
  if (has_treedepth) {
    line <- sprintf(
      "Treedepth was saturated on %.0f iteration(s), indicating slow posterior exploration. max_treedepth may be raised for diagnosis, but frequent saturation points to the parameterization rather than just the treedepth cap.",
      n_treedepth_sat
    )
    if (identical(bucket, "clean")) {
      bucket <- "treedepth"
      primary <- c(primary, line)
    } else {
      secondary <- c(secondary, line)
    }
  }

  # -- D. Fixed-effect Rhat failures, named by family when possible --
  if (beta_rhat_bad) {
    family_txt <- if (!is.null(beta_summary) && nrow(beta_summary) > 0L) {
      worst_family <- beta_summary[which.max(beta_summary$pct_rhat_gt_1_01), , drop = FALSE]
      if (!is.na(worst_family$pct_rhat_gt_1_01[[1L]]) && worst_family$pct_rhat_gt_1_01[[1L]] > 0) {
        sprintf("Convergence problems are concentrated in %s fixed effects (%.1f%% of that family exceeds Rhat %.2f, worst Rhat = %.4f).",
               worst_family$family[[1L]], worst_family$pct_rhat_gt_1_01[[1L]], thresholds$rhat_max, worst_family$worst_rhat[[1L]])
      } else NA_character_
    } else NA_character_
    line <- if (!is.na(family_txt)) {
      family_txt
    } else {
      sprintf("%.1f%% of fixed-effect (beta) parameters exceed Rhat %.2f.", beta_rhat_pct, thresholds$rhat_max)
    }
    if (identical(bucket, "clean")) {
      bucket <- "beta_rhat"
      primary <- c(primary, line)
    } else {
      secondary <- c(secondary, line)
    }
  }

  if (identical(bucket, "clean") && length(primary) == 0L) {
    primary <- c(primary, "No divergence, energy, treedepth, or Rhat concerns were found for this fit.")
  }

  list(primary = primary, secondary = secondary, not_indicated = not_indicated, bucket = bucket)
}

#' Compact diagnostic-summary table for the interpretation page
#'
#' @inheritParams .probit_classify_fitting_issues
#' @return Character vector of pre-formatted monospace table lines (header,
#'   separator, one row per diagnostic).
#' @keywords internal
.probit_diagnostic_table_lines <- function(diag, degeneracy_stats, residual_structure,
                                           thresholds = .probit_diagnostic_thresholds()) {
  n_divergent <- .probit_num1(diag$n_divergent)
  n_treedepth_sat <- .probit_num1(diag$n_treedepth_sat)
  ebfmi_min <- .probit_num1(diag$ebfmi_min)
  max_rhat <- .probit_num1(diag$max_rhat_structural)

  .verdict <- function(bad) if (isTRUE(bad)) "Fail" else if (isFALSE(bad)) "Pass" else "NA"

  rows <- list(
    list("Divergences", sprintf("%.0f", n_divergent), sprintf(">%d", thresholds$divergent_max),
        .verdict(!is.na(n_divergent) && n_divergent > thresholds$divergent_max)),
    list("Treedepth hits", sprintf("%.0f", n_treedepth_sat), sprintf(">%d", thresholds$treedepth_max),
        .verdict(!is.na(n_treedepth_sat) && n_treedepth_sat > thresholds$treedepth_max)),
    list("Min E-BFMI", sprintf("%.3f", ebfmi_min), sprintf("<%.2f", thresholds$ebfmi_min),
        .verdict(!is.na(ebfmi_min) && ebfmi_min < thresholds$ebfmi_min)),
    list("Max structural Rhat", sprintf("%.4f", max_rhat), sprintf(">%.2f", thresholds$rhat_max),
        .verdict(!is.na(max_rhat) && max_rhat > thresholds$rhat_max))
  )
  if (identical(residual_structure, "correlated") && !is.null(degeneracy_stats)) {
    pct <- degeneracy_stats$pct_near_degenerate
    rows <- c(rows, list(list(
      "Omega near-singular draws", sprintf("%.1f%%", pct),
      sprintf(">%.0f%%", 100 * thresholds$omega_small_eigen_fraction_warning),
      .verdict(!is.na(pct) && (pct / 100) > thresholds$omega_small_eigen_fraction_warning)
    )))
  }

  header <- sprintf("%-26s %10s %14s %14s", "Diagnostic", "Value", "Threshold", "Interpretation")
  body <- vapply(rows, function(r) sprintf("%-26s %10s %14s %14s", r[[1]], r[[2]], r[[3]], r[[4]]), character(1L))
  c(header, strrep("-", nchar(header)), body)
}
