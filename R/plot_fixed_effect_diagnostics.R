# plot_fixed_effect_diagnostics.R
# Fixed-effect (beta) coefficient diagnostics, separated from
# plot_probit_diagnostics()'s existing beta-density page (which only shows
# a max_params-truncated SAMPLE of coefficients, typically 10 of however
# many exist -- fine as a quick sampler-health glance, but never intended
# to be a complete view). With hospital as a fixed effect, K can easily
# exceed 20 (one dummy per hospital); plotting all of them in one page,
# unfaceted, the way a naive forest plot would, is unreadable, and mixing
# hospital contrasts with substantive covariates (age, gender, ...) in one
# list obscures which effects are actually of scientific interest.

#' Fixed-effect coefficient diagnostics, split by contrast type and faceted
#' by antimicrobial class
#'
#' Returns a named list of up to two ggplot objects: \code{other} (every
#' non-hospital-dummy covariate -- intercept, age, gender, etc.) and
#' \code{hospital} (hospital contrasts only, if the fit used hospital as a
#' fixed effect rather than a random effect). Both facet by antimicrobial
#' class rather than putting every class's coefficients on one axis.
#'
#' @param fit List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param title_base Character. Prefixed to plot titles.
#' @param ci_level Numeric. Credible interval level. Default \code{0.95}.
#' @return Named list of ggplot objects (\code{other}, \code{hospital}),
#'   either of which may be absent if that group has no covariates in this
#'   fit. \code{NULL} if \code{fit$draws} or \code{fit$X_event} is missing.
#' @export
plot_probit_fixed_effect_diagnostics <- function(fit, title_base = "", ci_level = 0.95) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("posterior", quietly = TRUE)) {
    return(NULL)
  }
  X <- fit$X_event %||% fit$X_design
  if (is.null(fit$draws) || is.null(X) || is.null(colnames(X))) return(NULL)

  covariate_names <- colnames(X)
  class_cols <- fit$class_cols
  K <- length(covariate_names)
  D <- length(class_cols)

  beta_vars <- sprintf("beta[%d,%d]", rep(seq_len(K), D), rep(seq_len(D), each = K))
  beta_vars <- intersect(beta_vars, posterior::variables(fit$draws))
  if (length(beta_vars) == 0L) return(NULL)

  summ <- posterior::summarise_draws(
    posterior::subset_draws(fit$draws, variable = beta_vars),
    "median", ~ stats::quantile(.x, probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2))
  )
  names(summ)[names(summ) == paste0(round(100 * (1 - ci_level) / 2, 1), "%")] <- "lower"
  names(summ)[names(summ) == paste0(round(100 * (1 - (1 - ci_level) / 2), 1), "%")] <- "upper"

  m <- regmatches(summ$variable, regexec("^beta\\[(\\d+),(\\d+)\\]$", summ$variable))
  summ$k <- as.integer(vapply(m, `[`, character(1L), 2L))
  summ$d <- as.integer(vapply(m, `[`, character(1L), 3L))
  summ$covariate <- covariate_names[summ$k]
  summ$class_display <- class_display_label(class_cols[summ$d])
  summ$is_hospital <- grepl("^center_name", summ$covariate)

  .plot_group <- function(df, subtitle_extra, y_lab) {
    if (nrow(df) == 0L) return(NULL)
    df$covariate_display <- factor(df$covariate_display,
                                   levels = sort(unique(df$covariate_display), decreasing = TRUE))
    ggplot2::ggplot(df, ggplot2::aes(y = .data$covariate_display)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                             height = 0, colour = "#4292C6", alpha = 0.7) +
      ggplot2::geom_point(ggplot2::aes(x = .data$median), colour = "#2166AC", size = 1.6) +
      ggplot2::facet_wrap(~class_display) +
      ggplot2::labs(
        title = paste(title_base, "-- Fixed-Effect Coefficients", subtitle_extra),
        subtitle = "Posterior median and 95% credible interval, probit scale. Dashed line = no effect.",
        x = "Probit-scale coefficient", y = y_lab
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11), strip.text = ggplot2::element_text(size = 8))
  }

  plots <- list()

  other_df <- summ[!summ$is_hospital, , drop = FALSE]
  if (nrow(other_df) > 0L) {
    other_df$covariate_display <- other_df$covariate
    plots$other <- .plot_group(other_df, "(Non-Hospital Covariates)", "Covariate / contrast")
  }

  hosp_df <- summ[summ$is_hospital, , drop = FALSE]
  if (nrow(hosp_df) > 0L) {
    hosp_df$covariate_display <- hospital_display_label(sub("^center_name", "", hosp_df$covariate))
    plots$hospital <- .plot_group(hosp_df, "(Hospital Contrasts)", "Hospital")
  }

  plots
}
