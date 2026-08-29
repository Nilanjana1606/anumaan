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
