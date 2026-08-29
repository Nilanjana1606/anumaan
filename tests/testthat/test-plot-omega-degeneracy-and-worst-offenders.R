# Tests for R/plot_omega_diagnostics.R's degeneracy/summary-table additions
# and R/plot_worst_offender_diagnostics.R.

.omega_draws_2x2 <- function(rho, n_draws = 50L, n_chains = 2L) {
  arr <- array(0, dim = c(n_draws, n_chains, 4),
              dimnames = list(NULL, NULL, c("Omega[1,1]", "Omega[2,1]", "Omega[1,2]", "Omega[2,2]")))
  arr[, , "Omega[1,1]"] <- 1
  arr[, , "Omega[2,2]"] <- 1
  arr[, , "Omega[2,1]"] <- rho
  arr[, , "Omega[1,2]"] <- rho
  posterior::as_draws_array(arr)
}

test_that(".omega_degeneracy_stats() computes the analytic smallest eigenvalue for a 2x2 correlation matrix (1 - rho)", {
  fit <- list(residual_structure = "correlated", draws = .omega_draws_2x2(0.6), class_cols = c("A", "B"))
  stats_out <- .omega_degeneracy_stats(fit, c("A", "B"))
  expect_false(is.null(stats_out))
  expect_true(all(abs(stats_out$lambda_min - 0.4) < 1e-6))
})

test_that(".omega_degeneracy_stats()/plot_omega_degeneracy_diagnostic() return NULL for identity fits or missing Omega", {
  identity_fit <- list(residual_structure = "identity", draws = .omega_draws_2x2(0.6), class_cols = c("A", "B"))
  expect_null(.omega_degeneracy_stats(identity_fit, c("A", "B")))
  expect_null(plot_omega_degeneracy_diagnostic(identity_fit, c("A", "B")))

  no_draws_fit <- list(residual_structure = "correlated", draws = NULL, class_cols = c("A", "B"))
  expect_null(.omega_degeneracy_stats(no_draws_fit, c("A", "B")))
})

test_that("plot_omega_degeneracy_diagnostic() flags near-degeneracy at a threshold above 1 - rho, not below it", {
  fit <- list(residual_structure = "correlated", draws = .omega_draws_2x2(0.9), class_cols = c("A", "B"))
  p <- plot_omega_degeneracy_diagnostic(fit, c("A", "B"), degenerate_threshold = 0.5)
  expect_false(is.null(p))
  expect_true(grepl("100\\.0% of draws", p$labels$subtitle))
  p_low <- plot_omega_degeneracy_diagnostic(fit, c("A", "B"), degenerate_threshold = 0.05)
  expect_true(grepl("0\\.0% of draws", p_low$labels$subtitle))
})

test_that("plot_omega_summary_table() returns NULL for empty input and includes ESS tail when present", {
  expect_null(plot_omega_summary_table(NULL, c("A", "B")))
  corr_summary <- data.frame(
    class_1 = "Aminoglycosides", class_2 = "Carbapenems",
    correlation_mean = 0.8, correlation_median = 0.8,
    correlation_lower = 0.7, correlation_upper = 0.9,
    rhat = 1.02, ess_bulk = 300, ess_tail = 500,
    stringsAsFactors = FALSE
  )
  p <- plot_omega_summary_table(corr_summary, c("Aminoglycosides", "Carbapenems"), "test")
  expect_false(is.null(p))
  lines <- ggplot2::ggplot_build(p)$data[[1]]$label
  expect_true(any(grepl("500", lines)))
})

.fake_worst_offender_fit <- (function() {
  vars <- c("beta[1,1]", "beta[2,1]", "Omega[1,2]", "lp__")
  arr <- array(rnorm(200 * length(vars)), dim = c(50, 4, length(vars)),
              dimnames = list(NULL, NULL, vars))
  draws <- posterior::as_draws_array(arr)
  monitored <- data.frame(
    variable = vars,
    parameter_group = c("beta", "beta", "Omega", "lp__"),
    rhat = c(1.08, 1.00, 1.02, 1.001),
    ess_bulk = c(60, 900, 150, 950),
    ess_tail = c(150, 1200, 400, 1300),
    stringsAsFactors = FALSE
  )
  X <- matrix(0, nrow = 2, ncol = 2, dimnames = list(NULL, c("(Intercept)", "Age_normalised")))
  list(
    draws = draws,
    X_event = X,
    class_cols = c("Aminoglycosides", "Carbapenems"),
    diagnostics_detail = list(monitored_parameters = monitored)
  )
})()

test_that("plot_probit_worst_offender_diagnostics() relabels the selected variables and returns trace+rank plots", {
  out <- plot_probit_worst_offender_diagnostics(.fake_worst_offender_fit, .fake_worst_offender_fit$draws, "test", n_worst = 2L)
  expect_false(is.null(out))
  expect_true(all(c("trace", "rank") %in% names(out)))
  expect_false(is.null(out$trace))
  expect_false(is.null(out$rank))
})

test_that("plot_probit_worst_offender_diagnostics() returns NULL when monitored_parameters is unavailable", {
  broken_fit <- .fake_worst_offender_fit
  broken_fit$diagnostics_detail <- list()
  expect_null(plot_probit_worst_offender_diagnostics(broken_fit, broken_fit$draws, "test"))
})

test_that(".probit_interpretation_text() picks the energy-failure branch over the rhat-mixing branch when both apply", {
  diag <- list(ebfmi_min = 0.1, n_divergent = 0)
  beta_summary <- data.frame(family = "Hospital", pct_rhat_gt_1_01 = 90, worst_rhat = 1.05, min_ess_bulk = 50)
  txt <- .probit_interpretation_text(diag, "fail_energy", "FAIL", beta_summary, NULL, "identity")
  expect_true(any(grepl("geometry/energy problem", txt)))
  expect_false(any(grepl("slow/inconsistent mixing", txt)))
})

test_that(".probit_interpretation_text() picks the rhat-mixing branch and mentions Omega only for correlated fits", {
  diag <- list(ebfmi_min = 0.8, n_divergent = 0)
  beta_summary <- data.frame(family = "Hospital", pct_rhat_gt_1_01 = 40, worst_rhat = 1.03, min_ess_bulk = 200)
  degeneracy_stats <- list(pct_near_degenerate = 5, med_condition_number = 10, degenerate_threshold = 0.05)

  txt_corr <- .probit_interpretation_text(diag, "warning_rhat", "WARNING", beta_summary, degeneracy_stats, "correlated")
  expect_true(any(grepl("slow/inconsistent mixing", txt_corr)))
  expect_true(any(grepl("^Omega", txt_corr)))

  txt_identity <- .probit_interpretation_text(diag, "warning_rhat", "WARNING", beta_summary, NULL, "identity")
  expect_false(any(grepl("^Omega", txt_identity)))
})

test_that(".probit_interpretation_text() names the specific worst-converging Omega pair when one exceeds the Rhat threshold", {
  diag <- list(ebfmi_min = 0.8, n_divergent = 0)
  beta_summary <- data.frame(family = "Hospital", pct_rhat_gt_1_01 = 5, worst_rhat = 1.005, min_ess_bulk = 500)
  degeneracy_stats <- list(pct_near_degenerate = 2, med_condition_number = 8, degenerate_threshold = 0.05)
  corr_summary_bad <- data.frame(
    class_1 = c("Aminoglycosides", "Penicillins"), class_2 = c("Fluoroquinolones", "Carbapenems"),
    rhat = c(1.05, 1.001), stringsAsFactors = FALSE
  )
  txt <- .probit_interpretation_text(diag, "pass", "PASS", beta_summary, degeneracy_stats, "correlated",
                                     corr_summary = corr_summary_bad)
  omega_line <- txt[grepl("^Omega", txt)]
  expect_true(grepl("least-converged Omega element is AMG x FQ", omega_line))

  corr_summary_ok <- corr_summary_bad
  corr_summary_ok$rhat <- c(1.005, 1.001)
  txt_ok <- .probit_interpretation_text(diag, "pass", "PASS", beta_summary, degeneracy_stats, "correlated",
                                        corr_summary = corr_summary_ok)
  omega_line_ok <- txt_ok[grepl("^Omega", txt_ok)]
  expect_false(grepl("least-converged Omega element", omega_line_ok))
})

test_that(".probit_interpretation_text() reports a clean pass with no concerns", {
  diag <- list(ebfmi_min = 0.8, n_divergent = 0)
  beta_summary <- data.frame(family = c("Hospital", "Age"), pct_rhat_gt_1_01 = c(0, 0),
                             worst_rhat = c(1.0, 1.0), min_ess_bulk = c(900, 900))
  txt <- .probit_interpretation_text(diag, "pass", "PASS", beta_summary, NULL, "identity")
  expect_true(any(grepl("No convergence, energy, or divergence concerns", txt)))
})
