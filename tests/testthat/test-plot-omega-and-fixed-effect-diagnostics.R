# Tests for R/plot_omega_diagnostics.R and R/plot_fixed_effect_diagnostics.R.

test_that("plot_omega_correlation_heatmap() and plot_omega_convergence_heatmap() return NULL for NULL/empty input", {
  expect_null(plot_omega_correlation_heatmap(NULL, c("A", "B")))
  expect_null(plot_omega_correlation_heatmap(data.frame(), c("A", "B")))
  expect_null(plot_omega_convergence_heatmap(NULL, c("A", "B")))
})

.fake_corr_summary <- data.frame(
  class_1 = c("Aminoglycosides", "Aminoglycosides", "Carbapenems"),
  class_2 = c("Carbapenems", "Penicillins", "Penicillins"),
  correlation_mean = c(0.9, 0.85, 0.95), correlation_median = c(0.9, 0.85, 0.95),
  correlation_lower = c(0.85, 0.8, 0.9), correlation_upper = c(0.95, 0.9, 0.98),
  rhat = c(1.005, 1.02, 1.001), ess_bulk = c(500, 200, 800),
  stringsAsFactors = FALSE
)

test_that("Omega correlation heatmap flags high-Rhat cells and never shows empirical-correlation language", {
  p <- plot_omega_correlation_heatmap(.fake_corr_summary, c("Aminoglycosides", "Penicillins", "Carbapenems"), "test")
  expect_false(is.null(p))
  expect_true(grepl("latent", p$labels$subtitle, ignore.case = TRUE))
  expect_false(grepl("empirical", p$labels$subtitle, ignore.case = TRUE) &&
                 !grepl("NOT empirical", p$labels$subtitle))
  # The 0.85 correlation (Aminoglycosides x Penicillins) has rhat = 1.02 > 1.01, must be flagged.
  flagged_row <- p$data[p$data$rhat > 1.01 & !is.na(p$data$rhat), ]
  expect_true(nrow(flagged_row) > 0L)
  expect_true(all(grepl("\\*$", flagged_row$label)))
})

test_that("Omega convergence heatmap surfaces rhat as the fill and includes all classes symmetrically", {
  p <- plot_omega_convergence_heatmap(.fake_corr_summary, c("Aminoglycosides", "Penicillins", "Carbapenems"), "test")
  expect_false(is.null(p))
  expect_identical(rlang::as_label(p$mapping$fill), "rhat")
  # Full symmetric matrix: 3 classes -> 9 cells (3x3), not just the 3 upper-triangle rows.
  expect_identical(nrow(p$data), 9L)
})

test_that("plot_probit_fixed_effect_diagnostics() returns NULL when draws/X_event are missing", {
  expect_null(plot_probit_fixed_effect_diagnostics(list(draws = NULL, X_event = NULL)))
})
