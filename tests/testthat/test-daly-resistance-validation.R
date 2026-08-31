# Tests for R/daly_resistance_validation.R, focused on the correlated-residual
# fix to validate_complete_profile_calibration(): that function used to
# compute full-profile probabilities as the independent product of Phi(mu_d)
# regardless of residual_structure, which is wrong once Omega has non-zero
# off-diagonal correlation (P(Y_1=1,...,Y_D=1) != prod_d Phi(mu_d) in that
# case). The fix estimates the model-implied profile probability via Monte
# Carlo simulation of Z = mu + L_Omega %*% eps, Y = I(Z > 0) -- the same
# mechanism .ppc_generate_correlated() already uses.
#
# Part 1 tests the core simulation property directly (Omega = I reduces to
# the independent-product formula; rho = 0.8 diverges from it in the
# expected direction) without needing a full fitted_model or cmdstanr.
# Part 2 is an end-to-end check against a real small correlated fit.

.sim_profile_freqs <- function(mu, Omega, M, seed) {
  set.seed(seed)
  D <- ncol(mu); n_ev <- nrow(mu)
  L <- t(chol(Omega))
  eps <- matrix(stats::rnorm(D * n_ev * M), nrow = D, ncol = n_ev * M)
  z <- t(L %*% eps) + mu[rep(seq_len(n_ev), each = M), , drop = FALSE]
  y <- z > 0
  labels <- apply(y, 1L, function(row) paste(ifelse(row, "R", "S"), collapse = ""))
  as.list(prop.table(table(factor(labels, levels = c("SS", "SR", "RS", "RR")))))
}

.indep_profile_freqs <- function(mu) {
  p <- stats::pnorm(mu)
  list(
    SS = mean((1 - p[, 1]) * (1 - p[, 2])),
    SR = mean((1 - p[, 1]) * p[, 2]),
    RS = mean(p[, 1] * (1 - p[, 2])),
    RR = mean(p[, 1] * p[, 2])
  )
}

test_that("correlated-profile simulation with Omega = I reduces to the independent-product formula", {
  set.seed(1)
  mu <- cbind(stats::rnorm(500, 0.2, 0.3), stats::rnorm(500, -0.1, 0.3))
  sim <- .sim_profile_freqs(mu, diag(2), M = 3000L, seed = 42)
  ind <- .indep_profile_freqs(mu)
  for (lbl in c("SS", "SR", "RS", "RR")) {
    expect_equal(sim[[lbl]], ind[[lbl]], tolerance = 0.01)
  }
})

test_that("correlated-profile simulation with rho = 0.8 diverges from independence in the expected direction", {
  set.seed(2)
  mu <- cbind(stats::rnorm(500, 0.2, 0.3), stats::rnorm(500, -0.1, 0.3))
  ind <- .indep_profile_freqs(mu)
  sim <- .sim_profile_freqs(mu, matrix(c(1, 0.8, 0.8, 1), 2, 2), M = 3000L, seed = 43)

  # Positive residual correlation should make CONCORDANT profiles (RR, SS)
  # more likely and DISCORDANT profiles (RS, SR) less likely than under
  # independence.
  expect_gt(sim$RR, ind$RR)
  expect_gt(sim$SS, ind$SS)
  expect_lt(sim$RS + sim$SR, ind$RS + ind$SR)
})

.validation_cmdstan_available <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({ cmdstanr::cmdstan_path(); TRUE }, error = function(e) FALSE))
}

test_that("validate_complete_profile_calibration() uses Omega for a real correlated fit, not the independent product", {
  skip_if_not(.validation_cmdstan_available(), "cmdstanr/CmdStan not available")

  # Simulate D=2 classes sharing a strong common latent factor, so the fitted
  # Omega should have real, non-trivial correlation -- if the fix is wired
  # correctly, validate_complete_profile_calibration()'s model_frequency_mean
  # for RR/SS should differ materially from what the independent-product
  # formula (Phi(mu_1) * Phi(mu_2), etc.) would give on the SAME fitted mu.
  set.seed(7)
  n <- 300L
  shared <- stats::rnorm(n, sd = 1.2)
  age <- stats::rnorm(n)
  wide <- data.frame(
    center_name = "H1", pathogen = "bug", Age_normalised = age,
    classA = stats::rbinom(n, 1, stats::pnorm(0.1 + 0.1 * age + shared)),
    classB = stats::rbinom(n, 1, stats::pnorm(-0.1 + 0.1 * age + shared)),
    event_id = paste0("e", seq_len(n))
  )

  fit <- suppressWarnings(fit_bayesian_multivariate_probit(
    event_class_data = wide, class_cols = c("classA", "classB"),
    fixed_effects = "Age_normalised", random_effects = list(),
    pathogen = "bug", pathogen_col = "pathogen", event_id_col = "event_id",
    profile_group_col = "center_name",
    residual_structure = "correlated",
    prior_config = list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0),
    sampler_config = list(chains = 2L, iter_warmup = 150L, iter_sampling = 150L,
                          seed = 7L, parallel_chains = 2L, adapt_delta = 0.9),
    show_messages = FALSE
  ))

  fixed_result <- validate_complete_profile_calibration(
    fit, n_posterior_draws_for_validation = 100L, seed = 1L,
    min_complete_events = 30L, n_mc_profile_replicates = 200L
  )
  fixed_result <- fixed_result[fixed_result$status == "evaluated", ]
  expect_true(nrow(fixed_result) > 0L)

  # Reconstruct what the OLD (buggy) independent-product code would have
  # produced for the exact same fit/draws, to confirm the fix actually
  # changed the answer rather than accidentally reproducing it.
  setup <- anumaan:::.probit_validation_draws_setup(fit, 100L, 1L)
  old_indep <- vapply(seq_len(setup$S), function(s) {
    p <- stats::pnorm(setup$mu_all_for_draw(s))
    mean(p[, 1] * p[, 2])
  }, numeric(1L))
  old_rr_mean <- mean(old_indep)
  new_rr_mean <- fixed_result$model_frequency_mean[fixed_result$profile_delta == "RR"][1]

  expect_false(isTRUE(all.equal(old_rr_mean, new_rr_mean, tolerance = 0.01)))
})
