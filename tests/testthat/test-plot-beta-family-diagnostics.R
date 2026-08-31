# Tests for R/plot_beta_family_diagnostics.R.

.fake_X <- matrix(0, nrow = 4, ncol = 5,
                  dimnames = list(NULL, c("(Intercept)", "Age_normalised", "genderMale",
                                           "center_nameHospital_A", "center_nameHospital_B")))

.fake_monitored <- (function() {
  vars <- sprintf("beta[%d,%d]", rep(1:5, 2), rep(1:2, each = 5))
  rhat <- c(1.00, 1.00, 1.00, 1.06, 1.00,   # class 1: only Hospital_A dummy (k=4) fails
            1.00, 1.00, 1.00, 1.00, 1.00)   # class 2: all pass
  data.frame(
    variable = c(vars, "Omega[1,2]", "lp__"),
    parameter_group = c(rep("beta", 10), "Omega", "lp__"),
    rhat = c(rhat, 1.02, 1.001),
    ess_bulk = c(rep(300, 10), 150, 900),
    ess_tail = c(rep(600, 10), 400, 1200),
    stringsAsFactors = FALSE
  )
})()

.fake_fit <- list(
  X_event = .fake_X,
  class_cols = c("Aminoglycosides", "Carbapenems"),
  diagnostics_detail = list(monitored_parameters = .fake_monitored)
)

test_that("plot_probit_beta_family_diagnostics() returns NULL when X_event/monitored_parameters missing", {
  expect_null(plot_probit_beta_family_diagnostics(list(X_event = NULL)))
  expect_null(plot_probit_beta_family_diagnostics(list(X_event = .fake_X, diagnostics_detail = list())))
})

test_that("plot_probit_beta_family_diagnostics() splits Hospital from Intercept/Age/Sex, not lumped into one beta group", {
  out <- plot_probit_beta_family_diagnostics(.fake_fit, "test")
  expect_false(is.null(out))
  expect_true(all(c("table", "bar") %in% names(out)))
  bar_data <- out$bar$data
  expect_setequal(as.character(bar_data$family), c("Intercept", "Age", "Sex", "Hospital"))
  # Only the Hospital family contains the one failing coefficient (rhat 1.06 out of 2 hospital
  # dummies x 2 classes = 4 hospital parameters -> 25% > 1.01); every other family is clean.
  hosp_row <- bar_data[bar_data$family == "Hospital", ]
  expect_equal(hosp_row$pct_rhat_gt_1_01, 25)
  other_rows <- bar_data[bar_data$family != "Hospital", ]
  expect_true(all(other_rows$pct_rhat_gt_1_01 == 0))
})

test_that("plot_probit_worst_parameters() shows only failing parameters with human-readable labels, not raw beta[k,d]/Omega[i,j]", {
  p <- plot_probit_worst_parameters(.fake_fit, "test", rhat_threshold = 1.01)
  expect_false(is.null(p))
  table_lines <- ggplot2::ggplot_build(p)$data[[1]]$label
  # The failing beta[4,1] parameter must be relabeled to covariate x class, not left as "beta[4,1]".
  expect_false(any(grepl("^beta\\[", table_lines)))
  # The Omega[1,2] parameter (rhat 1.02, also > 1.01) must be relabeled as a correlation pair.
  expect_true(any(grepl("^Correlation: ", table_lines)))
  expect_false(any(grepl("^Omega\\[", table_lines)))
  expect_true(any(grepl("Hospital A", table_lines, fixed = TRUE)))
})

test_that("plot_probit_worst_parameters() returns NULL when nothing exceeds the threshold", {
  clean_fit <- .fake_fit
  clean_fit$diagnostics_detail$monitored_parameters$rhat <- 1.00
  expect_null(plot_probit_worst_parameters(clean_fit, "test", rhat_threshold = 1.01))
})
