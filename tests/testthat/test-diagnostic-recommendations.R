# Tests for R/plot_diagnostic_recommendations.R -- the failure-mode-specific
# interpretation logic that replaced a single generic "increase adapt_delta
# or simplify the model" line (see that file's header comment for why the
# old behaviour was a category error, not just imprecise wording).

.fake_grouped <- function(beta_pct = 0, omega_pct = 0, lp_rhat = 1.00) {
  data.frame(
    parameter_group = c("beta", "Omega", "lp__"),
    n_parameters = c(88, 32, 1),
    pct_rhat_above_1_01 = c(beta_pct, omega_pct, if (lp_rhat > 1.01) 100 else 0),
    max_rhat = c(1.0 + beta_pct / 1000, 1.0 + omega_pct / 1000, lp_rhat),
    min_ess_bulk = c(300, 130, 120),
    stringsAsFactors = FALSE
  )
}

.fake_diag <- function(n_divergent = 0, n_treedepth_sat = 0, ebfmi_min = 0.8,
                       max_rhat_structural = 1.005) {
  list(n_divergent = n_divergent, n_treedepth_sat = n_treedepth_sat,
      ebfmi_min = ebfmi_min, max_rhat_structural = max_rhat_structural)
}

test_that("1. divergences > 0 -> adapt_delta recommendation IS allowed, as primary", {
  issues <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 5, ebfmi_min = 0.8),
    grouped_diag = .fake_grouped(), degeneracy_stats = NULL, residual_structure = "identity"
  )
  expect_identical(issues$bucket, "divergence")
  expect_true(any(grepl("adapt_delta", issues$primary, fixed = TRUE)))
  expect_false(any(grepl("adapt_delta", issues$not_indicated, fixed = TRUE)))
})

test_that("2. divergences == 0 and E-BFMI < 0.30 -> adapt_delta must NOT be the primary recommendation", {
  issues <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, ebfmi_min = 0.10),
    grouped_diag = .fake_grouped(), degeneracy_stats = NULL, residual_structure = "identity"
  )
  expect_identical(issues$bucket, "energy_or_omega_geometry")
  expect_false(any(grepl("adapt_delta", issues$primary, fixed = TRUE)))
  expect_true(any(grepl("adapt_delta", issues$not_indicated, fixed = TRUE)))
  expect_true(any(grepl("energy exploration|E-BFMI", issues$primary)))
})

test_that("3. treedepth saturation > 0 -> treedepth-specific recommendation (when it is the worst issue)", {
  issues <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, n_treedepth_sat = 12, ebfmi_min = 0.8),
    grouped_diag = .fake_grouped(), degeneracy_stats = NULL, residual_structure = "identity"
  )
  expect_identical(issues$bucket, "treedepth")
  expect_true(any(grepl("[Tt]reedepth was saturated", issues$primary)))
  expect_false(any(grepl("max_treedepth -- not indicated", issues$not_indicated, fixed = TRUE)))
})

test_that("4. fixed-effect Rhat failures only -> fixed-effect/sparsity recommendation, not generic 'simplify model'", {
  issues <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, ebfmi_min = 0.8),
    grouped_diag = .fake_grouped(beta_pct = 25, omega_pct = 0),
    degeneracy_stats = NULL, residual_structure = "identity",
    beta_summary = data.frame(family = c("Hospital", "Age"), pct_rhat_gt_1_01 = c(25, 0),
                              worst_rhat = c(1.03, 1.0), min_ess_bulk = c(200, 900))
  )
  expect_identical(issues$bucket, "beta_rhat")
  expect_true(any(grepl("Hospital", issues$primary)))
  expect_false(any(grepl("^simplify the model$", issues$primary)))
})

test_that("5. Omega Rhat failures + near-singular Omega -> correlation-geometry recommendation", {
  degeneracy <- list(pct_near_degenerate = 99.4, med_condition_number = 95.1, degenerate_threshold = 0.05)
  issues <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, ebfmi_min = 0.8),
    grouped_diag = .fake_grouped(beta_pct = 0, omega_pct = 95.2),
    degeneracy_stats = degeneracy, residual_structure = "correlated"
  )
  expect_identical(issues$bucket, "energy_or_omega_geometry")
  expect_true(any(grepl("residual correlation geometry", issues$primary)))
  expect_true(any(grepl("near-singular", issues$primary)))
})

test_that("6. convergence passes but PPC fails -> handled by the caller as model-misspecification, not by this classifier (out of scope: PPC unavailable at fit-diagnostic time)", {
  issues <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, ebfmi_min = 0.8, max_rhat_structural = 1.005),
    grouped_diag = .fake_grouped(beta_pct = 0, omega_pct = 0),
    degeneracy_stats = NULL, residual_structure = "identity"
  )
  expect_identical(issues$bucket, "clean")
  expect_true(any(grepl("No divergence, energy, treedepth, or Rhat concerns", issues$primary)))
})

test_that("7. PPC good but geometry fails -> the geometry classification still fires regardless (this classifier only sees sampler diagnostics, which is correct: PPC quality never overrides geometry status)", {
  issues <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, ebfmi_min = 0.10),
    grouped_diag = .fake_grouped(), degeneracy_stats = NULL, residual_structure = "identity"
  )
  expect_identical(issues$bucket, "energy_or_omega_geometry")
  # This is the assertion the interpretation-page text carries: good PPC never
  # substitutes for trustworthy sampling. See .probit_interpretation_text()'s
  # "energy_or_omega_geometry" next_step text.
})

test_that("8. Exp-17-like fixture: primary = energy/Omega geometry, adapt_delta not indicated as primary remedy", {
  degeneracy <- list(pct_near_degenerate = 99.4, med_condition_number = 95.1, degenerate_threshold = 0.05)
  grouped <- .fake_grouped(beta_pct = 23.9, omega_pct = 95.2, lp_rhat = 1.0505)
  issues <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, n_treedepth_sat = 0, ebfmi_min = 0.109, max_rhat_structural = 1.0505),
    grouped_diag = grouped, degeneracy_stats = degeneracy, residual_structure = "correlated",
    beta_summary = data.frame(family = c("Intercept", "Hospital", "Age", "Sex"),
                              pct_rhat_gt_1_01 = c(50, 25, 0, 0),
                              worst_rhat = c(1.0117, 1.0125, 1.0, 1.0),
                              min_ess_bulk = c(307, 320, 1363, 1352))
  )
  expect_identical(issues$bucket, "energy_or_omega_geometry")
  expect_true(any(grepl("energy exploration", issues$primary)))
  expect_true(any(grepl("near-singular residual correlation geometry", issues$primary)))
  # adapt_delta may be MENTIONED in the primary text (to explain why it won't
  # help), but must never be RECOMMENDED there -- and must show up explicitly
  # in not_indicated.
  expect_false(any(grepl("increas(e|ing) adapt_delta may help", issues$primary)))
  expect_true(any(grepl("adapt_delta", issues$not_indicated, fixed = TRUE)))
  # Fixed-effect improvement is real (23.9% beta failure) but secondary, not primary.
  expect_true(any(grepl("Intercept", issues$secondary)))
})

test_that(".probit_diagnostic_table_lines() renders the compact table with an Omega row only for correlated fits", {
  diag <- .fake_diag(n_divergent = 0, ebfmi_min = 0.109, max_rhat_structural = 1.0505)
  degeneracy <- list(pct_near_degenerate = 99.4, med_condition_number = 95.1, degenerate_threshold = 0.05)

  lines_corr <- .probit_diagnostic_table_lines(diag, degeneracy, "correlated")
  expect_true(any(grepl("Omega near-singular draws", lines_corr)))
  expect_true(any(grepl("99.4%", lines_corr)))
  expect_true(any(grepl("Fail", lines_corr)))

  lines_id <- .probit_diagnostic_table_lines(diag, NULL, "identity")
  expect_false(any(grepl("Omega", lines_id)))
})

test_that("thresholds are configurable and actually change the classification", {
  # A 0.15 E-BFMI would normally fail against the default 0.30 threshold;
  # lowering the threshold to 0.10 should flip it to not being an issue.
  issues_default <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, ebfmi_min = 0.15),
    grouped_diag = .fake_grouped(), degeneracy_stats = NULL, residual_structure = "identity"
  )
  expect_identical(issues_default$bucket, "energy_or_omega_geometry")

  issues_relaxed <- .probit_classify_fitting_issues(
    diag = .fake_diag(n_divergent = 0, ebfmi_min = 0.15),
    grouped_diag = .fake_grouped(), degeneracy_stats = NULL, residual_structure = "identity",
    thresholds = .probit_diagnostic_thresholds(ebfmi_min = 0.10)
  )
  expect_identical(issues_relaxed$bucket, "clean")
})
