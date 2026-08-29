# plot_beta_family_diagnostics.R
# Fixed-effect (beta) convergence broken out by SOURCE OF VARIATION (hospital,
# age, sex, specimen, ...) rather than lumped into one "beta" parameter group.
# fit$diagnostics_detail$grouped already reports Rhat/ESS quantiles per
# top-level parameter_group (beta/Omega/z_free/lp__/RE blocks -- see
# .probit_parameter_group() in daly_resistance_profiles.R), but every
# fixed-effect coefficient -- one intercept, one age slope, one sex contrast,
# and one dummy per hospital -- collapses into the single "beta" row there,
# which hides whether e.g. hospital contrasts specifically are the problem.
# This file adds that finer split for beta only (the other parameter_group
# rows are already appropriately granular); it is a display/grouping
# breakdown of quantities already computed at fit time, nothing here
# recomputes a diagnostic value.

.covariate_family <- function(covariate) {
  fam <- rep("Other", length(covariate))
  fam[grepl("^\\(Intercept\\)$", covariate)] <- "Intercept"
  fam[grepl("^center_name", covariate)] <- "Hospital"
  fam[grepl("^age", covariate, ignore.case = TRUE)] <- "Age"
  fam[grepl("^(gender|sex)", covariate, ignore.case = TRUE)] <- "Sex"
  fam[grepl("^(specimen|sample_type|source)", covariate, ignore.case = TRUE)] <- "Specimen"
  fam[grepl("^(icu|ward|admission_type)", covariate, ignore.case = TRUE)] <- "Ward/ICU"
  fam
}

.beta_parameter_table <- function(fit) {
  X <- fit$X_event %||% fit$X_design
  monitored <- fit$diagnostics_detail$monitored_parameters
  if (is.null(X) || is.null(colnames(X)) || is.null(monitored) || nrow(monitored) == 0L) return(NULL)
  beta_rows <- monitored[monitored$parameter_group == "beta", , drop = FALSE]
  if (nrow(beta_rows) == 0L) return(NULL)
  covariate_names <- colnames(X)
  class_cols <- fit$class_cols
  K <- length(covariate_names)

  m <- regmatches(beta_rows$variable, regexec("^beta\\[(\\d+),(\\d+)\\]$", beta_rows$variable))
  k <- suppressWarnings(as.integer(vapply(m, `[`, character(1L), 2L)))
  d <- suppressWarnings(as.integer(vapply(m, `[`, character(1L), 3L)))
  ok <- !is.na(k) & !is.na(d) & k <= K & d <= length(class_cols)
  beta_rows <- beta_rows[ok, , drop = FALSE]
  k <- k[ok]; d <- d[ok]
  if (nrow(beta_rows) == 0L) return(NULL)

  beta_rows$covariate <- covariate_names[k]
  beta_rows$family <- .covariate_family(beta_rows$covariate)
  beta_rows$class_display <- class_display_label(class_cols[d])
  beta_rows$rhat <- suppressWarnings(as.numeric(beta_rows$rhat))
  beta_rows$ess_bulk <- suppressWarnings(as.numeric(beta_rows$ess_bulk))
  beta_rows$ess_tail <- suppressWarnings(as.numeric(beta_rows$ess_tail))
  beta_rows
}

#' Fixed-effect convergence broken out by source of variation
#'
#' Splits the fixed-effect ("beta") parameter group by which covariate
#' family each coefficient belongs to (Intercept / Hospital / Age / Sex /
#' Specimen / Ward-ICU / Other, by covariate-name pattern), so a reviewer can
#' see directly whether e.g. hospital contrasts specifically are destabilising
#' the fit, rather than only a single aggregate Rhat/ESS for all of "beta".
#'
#' @param fit List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param title_base Character. Prefixed to plot titles.
#' @return Named list with \code{table} (a formatted summary page) and
#'   \code{bar} (a percent-Rhat-above-1.01-by-family bar chart), or
#'   \code{NULL} if fixed-effect parameter diagnostics are unavailable.
#' @export
plot_probit_beta_family_diagnostics <- function(fit, title_base = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  tbl <- .beta_parameter_table(fit)
  if (is.null(tbl)) return(NULL)

  summary_tbl <- do.call(rbind, lapply(split(tbl, tbl$family), function(g) {
    data.frame(
      family = g$family[[1L]],
      n_parameters = nrow(g),
      pct_rhat_gt_1_01 = round(100 * mean(g$rhat > 1.01, na.rm = TRUE), 1),
      pct_rhat_gt_1_05 = round(100 * mean(g$rhat > 1.05, na.rm = TRUE), 1),
      worst_rhat = round(max(g$rhat, na.rm = TRUE), 4),
      min_ess_bulk = round(min(g$ess_bulk, na.rm = TRUE), 0),
      min_ess_tail = round(min(g$ess_tail, na.rm = TRUE), 0),
      stringsAsFactors = FALSE
    )
  }))
  summary_tbl <- summary_tbl[order(-summary_tbl$pct_rhat_gt_1_01, -summary_tbl$worst_rhat), , drop = FALSE]

  header <- sprintf("%-12s %6s %10s %10s %10s %10s %10s",
                    "Family", "n", "%Rhat>.01", "%Rhat>.05", "worstRhat", "minESSblk", "minESStl")
  body <- sprintf("%-12s %6d %10.1f %10.1f %10.4f %10.0f %10.0f",
                  summary_tbl$family, summary_tbl$n_parameters,
                  summary_tbl$pct_rhat_gt_1_01, summary_tbl$pct_rhat_gt_1_05,
                  summary_tbl$worst_rhat, summary_tbl$min_ess_bulk, summary_tbl$min_ess_tail)
  table_lines <- c(header, strrep("-", nchar(header)), body)

  table_page <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = rev(seq_along(table_lines)), label = table_lines,
                      hjust = 0, size = 3.4, family = "mono") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, length(table_lines) + 1) +
    ggplot2::labs(
      title = paste(title_base, "-- Fixed-Effect Convergence by Source of Variation"),
      subtitle = "Every beta coefficient split by the covariate family it belongs to (not lumped into one 'beta' aggregate) -- answers which part of the fixed-effect structure is destabilising the fit."
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11), plot.subtitle = ggplot2::element_text(size = 8))

  summary_tbl$family <- factor(summary_tbl$family, levels = rev(summary_tbl$family))
  bar_page <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = .data$family, y = .data$pct_rhat_gt_1_01)) +
    ggplot2::geom_col(fill = "#f28e2b") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(title_base, "-- % of Fixed-Effect Coefficients with Rhat > 1.01, by Family"),
      x = NULL, y = "% of parameters in family"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

  list(table = table_page, bar = bar_page, summary = summary_tbl)
}

.human_readable_parameter_label <- function(variable, parameter_group, fit) {
  out <- variable
  is_beta <- parameter_group == "beta" & grepl("^beta\\[\\d+,\\d+\\]$", variable)
  if (any(is_beta)) {
    X <- fit$X_event %||% fit$X_design
    covariate_names <- if (!is.null(X)) colnames(X) else NULL
    class_cols <- fit$class_cols
    if (!is.null(covariate_names) && !is.null(class_cols)) {
      m <- regmatches(variable[is_beta], regexec("^beta\\[(\\d+),(\\d+)\\]$", variable[is_beta]))
      k <- suppressWarnings(as.integer(vapply(m, `[`, character(1L), 2L)))
      d <- suppressWarnings(as.integer(vapply(m, `[`, character(1L), 3L)))
      ok <- !is.na(k) & !is.na(d) & k <= length(covariate_names) & d <= length(class_cols)
      idx <- which(is_beta)[ok]
      covariate_disp <- covariate_names[k[ok]]
      is_hosp <- grepl("^center_name", covariate_disp)
      covariate_disp[is_hosp] <- hospital_display_label(sub("^center_name", "", covariate_disp[is_hosp]))
      out[idx] <- sprintf("%s x %s", covariate_disp, class_display_label(class_cols[d[ok]]))
    }
  }
  is_omega <- parameter_group == "Omega" & grepl("^Omega\\[\\d+,\\d+\\]$", variable)
  if (any(is_omega)) {
    class_cols <- fit$class_cols
    if (!is.null(class_cols)) {
      m <- regmatches(variable[is_omega], regexec("^Omega\\[(\\d+),(\\d+)\\]$", variable[is_omega]))
      i <- suppressWarnings(as.integer(vapply(m, `[`, character(1L), 2L)))
      j <- suppressWarnings(as.integer(vapply(m, `[`, character(1L), 3L)))
      ok <- !is.na(i) & !is.na(j) & i <= length(class_cols) & j <= length(class_cols) & i != j
      idx <- which(is_omega)[ok]
      out[idx] <- sprintf("Correlation: %s x %s", class_display_label(class_cols[i[ok]]), class_display_label(class_cols[j[ok]]))
    }
  }
  out[parameter_group == "lp__"] <- "Log posterior density"
  out
}

#' Human-readable worst-converging structural parameters
#'
#' Replaces raw Stan parameter names (e.g. \code{"beta[12,2]"}) with the
#' substantive covariate/class (or correlation pair) they refer to, and shows
#' only parameters that actually exceed the Rhat threshold -- not an
#' arbitrary fixed top-N -- alongside their bulk and tail ESS: a high Rhat
#' with healthy ESS is a different diagnostic situation than a high Rhat with
#' very low ESS, and showing only one of the two obscures that distinction.
#'
#' @param fit List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param title_base Character. Prefixed to plot titles.
#' @param rhat_threshold Numeric. Only parameters with Rhat above this are shown. Default \code{1.01}.
#' @param max_rows Integer cap on rows shown, in case very many parameters fail. Default \code{40}.
#' @return A ggplot object (table page), or \code{NULL} if no structural
#'   parameter exceeds the threshold, or diagnostics are unavailable.
#' @export
plot_probit_worst_parameters <- function(fit, title_base = "", rhat_threshold = 1.01, max_rows = 40L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  monitored <- fit$diagnostics_detail$monitored_parameters
  if (is.null(monitored) || nrow(monitored) == 0L) return(NULL)
  monitored$rhat <- suppressWarnings(as.numeric(monitored$rhat))
  failing <- monitored[!is.na(monitored$rhat) & monitored$rhat > rhat_threshold, , drop = FALSE]
  if (nrow(failing) == 0L) return(NULL)
  failing <- failing[order(-failing$rhat), , drop = FALSE]
  n_total_failing <- nrow(failing)
  failing <- head(failing, max_rows)

  failing$label <- .human_readable_parameter_label(failing$variable, failing$parameter_group, fit)
  failing$ess_bulk <- suppressWarnings(as.numeric(failing$ess_bulk))
  failing$ess_tail <- suppressWarnings(as.numeric(failing$ess_tail))

  header <- sprintf("%-45s %10s %10s %10s", "Parameter", "Rhat", "ESSbulk", "ESStail")
  body <- sprintf("%-45s %10.4f %10.0f %10.0f", failing$label, failing$rhat, failing$ess_bulk, failing$ess_tail)
  table_lines <- c(header, strrep("-", nchar(header)), body)
  subtitle <- sprintf("%d structural parameter(s) exceed Rhat > %.2f%s",
                      n_total_failing, rhat_threshold,
                      if (n_total_failing > max_rows) sprintf(" (showing worst %d)", max_rows) else "")

  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = rev(seq_along(table_lines)), label = table_lines,
                      hjust = 0, size = 3.0, family = "mono") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, length(table_lines) + 1) +
    ggplot2::labs(
      title = paste(title_base, "-- Parameters with the Poorest Chain Convergence"),
      subtitle = subtitle
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11), plot.subtitle = ggplot2::element_text(size = 8))
}
