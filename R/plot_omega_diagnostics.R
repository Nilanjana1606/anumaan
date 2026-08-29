# plot_omega_diagnostics.R
# Latent cross-class correlation (Omega) diagnostic figures for correlated-
# residual multivariate probit fits. Deliberately kept as TWO separate
# plots rather than one: a visually attractive correlation heatmap is
# meaningless if the sampler never agreed on those values, so convergence
# is shown as its own figure a reader cannot skip past.

#' Mirror an upper-triangle class_1/class_2 correlation table into a full
#' symmetric matrix (both (i,j) and (j,i), plus a diagonal of 1s/NA), for a
#' complete-looking heatmap instead of a half-populated triangle.
#' @keywords internal
.omega_full_matrix <- function(corr_summary, classes, diagonal_value) {
  mirrored <- corr_summary
  mirrored$class_1 <- corr_summary$class_2
  mirrored$class_2 <- corr_summary$class_1
  diag_tbl <- data.frame(class_1 = classes, class_2 = classes, stringsAsFactors = FALSE)
  for (col in setdiff(names(corr_summary), c("class_1", "class_2"))) diag_tbl[[col]] <- NA_real_
  diag_tbl[[diagonal_value]] <- 1
  full <- dplyr::bind_rows(corr_summary, mirrored, diag_tbl)
  full$class_1 <- factor(class_display_label(full$class_1), levels = class_display_label(classes))
  full$class_2 <- factor(class_display_label(full$class_2), levels = rev(class_display_label(classes)))
  full
}

#' Posterior-median latent cross-class correlation heatmap (Omega)
#'
#' Rows/columns are antimicrobial classes; each cell is the posterior
#' median of the corresponding \code{Omega[i,j]} (or \code{R_block[r,i,j]})
#' entry. Only meaningful for correlated-residual fits -- pass \code{NULL}
#' (e.g. an identity-residual fit's \code{summarize_fit_correlation_matrix()}
#' result) to get \code{NULL} back rather than an empty/misleading plot.
#'
#' @param corr_summary Tibble from \code{summarize_fit_correlation_matrix()}
#'   (with \code{class_1, class_2, correlation_median, rhat} columns), or
#'   \code{NULL}.
#' @param class_cols Character vector of all class names, in the fit's
#'   canonical order (for consistent row/column ordering).
#' @param title_base Character. Prefixed to the plot title.
#' @param rhat_flag_threshold Numeric. Cells with \code{rhat} above this are
#'   marked with an asterisk, cross-referencing the companion convergence
#'   heatmap rather than repeating full diagnostics on this plot.
#' @return A ggplot object, or \code{NULL} if \code{corr_summary} is
#'   \code{NULL}/empty.
#' @export
plot_omega_correlation_heatmap <- function(corr_summary, class_cols, title_base = "",
                                           rhat_flag_threshold = 1.01) {
  if (is.null(corr_summary) || nrow(corr_summary) == 0L) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

  full <- .omega_full_matrix(corr_summary, class_cols, "correlation_median")
  full$label <- ifelse(is.na(full$correlation_median), "1.00", sprintf("%.2f", full$correlation_median))
  full$flag <- !is.na(full$rhat) & full$rhat > rhat_flag_threshold
  full$label <- ifelse(full$flag, paste0(full$label, "*"), full$label)

  ggplot2::ggplot(full, ggplot2::aes(x = .data$class_1, y = .data$class_2, fill = .data$correlation_median)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-1, 1), na.value = "grey90",
                                  name = "Posterior\nmedian rho") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste(title_base, "-- Latent Within-Event Cross-Class Correlation (Omega)"),
      subtitle = "Posterior median latent correlations -- NOT empirical correlations of observed resistance.",
      caption = paste(
        "Z is a latent, unobserved probit-scale variable; Omega describes correlation among Z, not among the observed AST results.",
        if (any(full$flag, na.rm = TRUE))
          sprintf("\n* = Rhat > %.2f -- see the companion convergence heatmap; do not treat that cell as a settled estimate.", rhat_flag_threshold)
        else NULL
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8, face = "italic"),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)
    )
}

#' Sampling-diagnostic (Rhat) heatmap for Omega
#'
#' Companion to \code{plot_omega_correlation_heatmap()} -- shows Rhat for
#' each \code{Omega[i,j]} entry, so a reader cannot look at an attractive
#' correlation heatmap without also seeing whether those specific
#' parameters actually converged. A uniformly high correlation heatmap
#' alongside a uniformly poor Rhat heatmap here is the signature of a
#' near-degenerate, weakly-identified correlation matrix, not a settled
#' finding.
#'
#' @inheritParams plot_omega_correlation_heatmap
#' @return A ggplot object, or \code{NULL} if \code{corr_summary} is
#'   \code{NULL}/empty or has no \code{rhat} column (e.g. computed by an
#'   older \code{summarize_fit_correlation_matrix()} call).
#' @export
plot_omega_convergence_heatmap <- function(corr_summary, class_cols, title_base = "") {
  if (is.null(corr_summary) || nrow(corr_summary) == 0L || !"rhat" %in% names(corr_summary)) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

  full <- .omega_full_matrix(corr_summary, class_cols, "rhat")
  full$label <- ifelse(is.na(full$rhat), "--", sprintf("%.3f", full$rhat))

  ggplot2::ggplot(full, ggplot2::aes(x = .data$class_1, y = .data$class_2, fill = .data$rhat)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", limits = c(1, NA),
                                 na.value = "grey90", name = "Rhat") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste(title_base, "-- Sampling Diagnostics for Omega"),
      subtitle = "Rhat per Omega[i,j] entry.",
      caption = "Values above ~1.01 mean the chains have not agreed on that correlation -- do not interpret the corresponding cell in the correlation heatmap as a settled estimate.",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8, face = "italic"),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)
    )
}

#' Full per-pair Omega summary table (rho, 95\% CrI, Rhat, ESS bulk/tail)
#'
#' The two heatmaps show the correlation estimate and its Rhat as separate
#' pictures; this page puts every number for every class pair -- median rho,
#' credible interval, Rhat, bulk and tail ESS -- in one table, for a reader
#' who wants the exact figures rather than reading them off a colour scale.
#'
#' @inheritParams plot_omega_correlation_heatmap
#' @return A ggplot object (table page), or \code{NULL} if \code{corr_summary}
#'   is \code{NULL}/empty.
#' @export
plot_omega_summary_table <- function(corr_summary, class_cols, title_base = "") {
  if (is.null(corr_summary) || nrow(corr_summary) == 0L) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

  tbl <- corr_summary
  tbl$pair <- sprintf("%s x %s", class_short_label(tbl$class_1), class_short_label(tbl$class_2))
  tbl <- tbl[order(-tbl$rhat), , drop = FALSE]
  ess_tail <- if ("ess_tail" %in% names(tbl)) tbl$ess_tail else rep(NA_real_, nrow(tbl))

  header <- sprintf("%-16s %8s %18s %10s %10s %10s", "Pair", "rho", "95% CrI", "Rhat", "ESSbulk", "ESStail")
  body <- sprintf("%-16s %8.3f %8.3f,%8.3f %10.4f %10.0f %10.0f",
                  tbl$pair, tbl$correlation_median, tbl$correlation_lower, tbl$correlation_upper,
                  tbl$rhat, tbl$ess_bulk, ess_tail)
  table_lines <- c(header, strrep("-", nchar(header)), body)

  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = rev(seq_along(table_lines)), label = table_lines,
                      hjust = 0, size = 3.2, family = "mono") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, length(table_lines) + 1) +
    ggplot2::labs(
      title = paste(title_base, "-- Omega Pairwise Correlation Summary"),
      subtitle = "Posterior median latent correlation, 95% credible interval, and per-pair convergence, for every class pair."
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11), plot.subtitle = ggplot2::element_text(size = 8))
}

#' Smallest-eigenvalue/condition-number stats for the posterior Omega matrix
#' @keywords internal
.omega_degeneracy_stats <- function(fit, class_cols, degenerate_threshold = 0.05) {
  if (!requireNamespace("posterior", quietly = TRUE)) return(NULL)
  if (!identical(fit$residual_structure, "correlated") || is.null(fit$draws)) return(NULL)
  D <- length(class_cols)
  if (D < 2L) return(NULL)
  omega_vars <- sprintf("Omega[%d,%d]", rep(seq_len(D), D), rep(seq_len(D), each = D))
  omega_vars <- intersect(omega_vars, posterior::variables(fit$draws))
  if (length(omega_vars) < D * D) return(NULL)

  mat <- as.matrix(posterior::as_draws_matrix(posterior::subset_draws(fit$draws, variable = omega_vars)))
  mat <- mat[, omega_vars, drop = FALSE]
  n_draws <- nrow(mat)
  eig <- vapply(seq_len(n_draws), function(i) {
    M <- matrix(mat[i, ], nrow = D, ncol = D)
    ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    c(lambda_min = min(ev), lambda_max = max(ev))
  }, numeric(2))
  lambda_min <- eig["lambda_min", ]
  lambda_max <- eig["lambda_max", ]
  condition_number <- lambda_max / pmax(lambda_min, .Machine$double.eps)
  list(
    lambda_min = lambda_min,
    pct_near_degenerate = round(100 * mean(lambda_min < degenerate_threshold), 1),
    med_condition_number = stats::median(condition_number),
    degenerate_threshold = degenerate_threshold
  )
}

#' Near-degeneracy diagnostic for the posterior Omega correlation matrix
#'
#' A multivariate-probit correlated-residual fit can push Omega toward a
#' near-singular boundary (very high pairwise correlations across the board)
#' without that being visible from any single Omega\[i,j\] entry. For each
#' posterior draw this computes the smallest eigenvalue of the full D x D
#' Omega matrix -- a value approaching 0 means that specific draw's
#' correlation matrix is nearly singular -- and reports the distribution
#' across draws, plus the condition number (largest / smallest eigenvalue)
#' as a supplementary summary.
#'
#' @param fit List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param class_cols Character vector of all class names, in the fit's
#'   canonical order.
#' @param title_base Character. Prefixed to the plot title.
#' @param degenerate_threshold Numeric. Draws with smallest eigenvalue below
#'   this are counted as "near-degenerate" in the subtitle. Default \code{0.05}.
#' @return A ggplot object, or \code{NULL} if this is not a correlated-residual
#'   fit, or Omega draws are unavailable.
#' @export
plot_omega_degeneracy_diagnostic <- function(fit, class_cols, title_base = "", degenerate_threshold = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  stats_out <- .omega_degeneracy_stats(fit, class_cols, degenerate_threshold)
  if (is.null(stats_out)) return(NULL)
  lambda_min <- stats_out$lambda_min
  pct_near_degenerate <- stats_out$pct_near_degenerate
  med_cond <- stats_out$med_condition_number

  ggplot2::ggplot(data.frame(lambda_min = lambda_min), ggplot2::aes(x = .data$lambda_min)) +
    ggplot2::geom_histogram(bins = 40, fill = "#B2182B", colour = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = degenerate_threshold, colour = "black", linetype = "dashed") +
    ggplot2::labs(
      title = paste(title_base, "-- Smallest Eigenvalue of the Posterior Omega Matrix"),
      subtitle = sprintf(
        "%.1f%% of draws have smallest eigenvalue < %.2f (near-singular Omega) | median condition number = %.1f",
        pct_near_degenerate, degenerate_threshold, med_cond
      ),
      caption = "Values approaching 0 mean that draw's D x D correlation matrix is nearly singular -- the sampler is repeatedly visiting a near-degenerate boundary of the correlation-matrix space, not a settled interior estimate.",
      x = "Smallest eigenvalue of Omega (per draw)", y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8, face = "italic")
    )
}
