# plot_worst_offender_diagnostics.R
# Rhat/ESS say chains disagree, but not WHY -- a high Rhat can mean slow
# mixing, chain separation, a sticky region, or genuine multimodality, and
# those call for different fixes. Trace/rank plots distinguish them visually.
# Producing them for every parameter is unreadable (100+ beta coefficients);
# this auto-selects only the worst-converging structural parameters -- top-N
# by Rhat and top-N by lowest bulk ESS, deduplicated -- and relabels them
# with the same human-readable names used on the worst-parameters table
# page, so the trace/rank facet strips read as substance, not "beta[12,2]".

.worst_offender_variables <- function(monitored, n_worst) {
  if (is.null(monitored) || nrow(monitored) == 0L) return(character(0))
  monitored$rhat <- suppressWarnings(as.numeric(monitored$rhat))
  monitored$ess_bulk <- suppressWarnings(as.numeric(monitored$ess_bulk))
  by_rhat <- monitored[order(-monitored$rhat), , drop = FALSE]
  by_ess <- monitored[order(monitored$ess_bulk), , drop = FALSE]
  unique(c(
    utils::head(by_rhat$variable[!is.na(by_rhat$rhat)], n_worst),
    utils::head(by_ess$variable[!is.na(by_ess$ess_bulk)], n_worst)
  ))
}

#' Trace and rank plots for the worst-converging structural parameters
#'
#' Auto-selects the union of the top \code{n_worst} parameters by Rhat and
#' the top \code{n_worst} by lowest bulk ESS (from
#' \code{fit$diagnostics_detail$monitored_parameters}), and shows their
#' chain trace and rank-overlay plots with human-readable facet labels
#' (covariate x class, or correlation pair) instead of raw Stan parameter
#' names -- enough to tell slow-mixing, chain separation, a sticky region,
#' and multimodality apart, which a single Rhat number cannot.
#'
#' @param fit List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param draws A \code{posterior::draws_array} (typically \code{fit$draws}).
#' @param title_base Character. Prefixed to plot titles.
#' @param n_worst Integer. How many parameters to select by each criterion
#'   (Rhat, bulk ESS) before deduplicating. Default \code{6}.
#' @return Named list with \code{trace} and \code{rank} ggplot objects, or
#'   \code{NULL} if no structural parameter diagnostics are available.
#' @export
plot_probit_worst_offender_diagnostics <- function(fit, draws, title_base = "", n_worst = 6L) {
  if (!requireNamespace("bayesplot", quietly = TRUE) || !requireNamespace("posterior", quietly = TRUE)) return(NULL)
  monitored <- fit$diagnostics_detail$monitored_parameters
  worst_vars <- .worst_offender_variables(monitored, n_worst)
  worst_vars <- intersect(worst_vars, posterior::variables(draws))
  if (length(worst_vars) == 0L) return(NULL)

  sub_draws <- posterior::subset_draws(draws, variable = worst_vars)
  parameter_group <- monitored$parameter_group[match(worst_vars, monitored$variable)]
  human_labels <- .human_readable_parameter_label(worst_vars, parameter_group, fit)
  human_labels <- make.unique(human_labels, sep = " #")
  rename_map <- stats::setNames(worst_vars, human_labels)
  sub_draws <- posterior::rename_variables(sub_draws, !!!rename_map)

  list(
    trace = bayesplot::mcmc_trace(sub_draws) +
      ggplot2::labs(
        title = paste(title_base, "-- Chain Traces for the Worst-Converging Parameters"),
        subtitle = sprintf(
          "Top %d by Rhat, top %d by lowest bulk ESS (deduplicated, %d shown) -- distinguishes slow mixing, chain separation, sticky regions, and multimodality, which a single Rhat number cannot.",
          n_worst, n_worst, length(worst_vars)
        )
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10), plot.subtitle = ggplot2::element_text(size = 7.5)),
    rank = bayesplot::mcmc_rank_overlay(sub_draws) +
      ggplot2::labs(
        title = paste(title_base, "-- Rank Plots for the Worst-Converging Parameters"),
        subtitle = "Uniform distribution across ranks indicates good chain mixing; a chain's rank histogram skewed high or low indicates that chain is systematically above/below the others."
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10), plot.subtitle = ggplot2::element_text(size = 7.5))
  )
}
