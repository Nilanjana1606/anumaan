# Tests for R/probit_prior_predictive.R: simulate_probit_prior_predictive(),
# compute_prior_predictive_status(), compute_prior_predictive_fingerprint().
# All require a real (small/fast) cmdstanr fit as the design/prior-config
# source, but the prior predictive SIMULATION itself never runs NUTS -- only
# a fixed_param LKJ sampler plus base-R RNG draws.

.priorppc_cmdstan_available <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({ cmdstanr::cmdstan_path(); TRUE }, error = function(e) FALSE))
}

.priorppc_make_fit <- function(seed = 1L) {
  set.seed(seed)
  n_hosp <- 3L; n_per_hosp <- 40L
  hosp_eff <- stats::rnorm(n_hosp, sd = 0.3)
  wide <- do.call(rbind, lapply(seq_len(n_hosp), function(h) {
    n <- n_per_hosp
    age <- stats::rnorm(n)
    mu1 <- 0.2 + 0.3 * age + hosp_eff[h]
    mu2 <- -0.1 + 0.2 * age + hosp_eff[h]
    data.frame(
      center_name = paste0("H", h), pathogen = "bug", Age_normalised = age,
      classA = stats::rbinom(n, 1, stats::pnorm(mu1)),
      classB = stats::rbinom(n, 1, stats::pnorm(mu2)),
      event_id = paste0("H", h, "_", seq_len(n))
    )
  }))
  suppressWarnings(fit_bayesian_multivariate_probit(
    event_class_data = wide, class_cols = c("classA", "classB"),
    fixed_effects = "Age_normalised", random_effects = c("center_name"),
    pathogen = "bug", pathogen_col = "pathogen", event_id_col = "event_id",
    residual_structure = "identity",
    prior_config = list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0),
    sampler_config = list(chains = 2L, iter_warmup = 150L, iter_sampling = 150L,
                          seed = seed, parallel_chains = 2L, adapt_delta = 0.9),
    show_messages = FALSE
  ))
}

test_that("LKJ Cholesky-factor prior draws are valid across a range of eta", {
  skip_if_not(.priorppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  for (eta in c(1.0, 2.0, 5.0)) {
    arr <- .draw_lkj_cholesky_prior(D = 3L, eta = eta, n_matrices = 2L, n_states = 4L, seed = 1L)
    expect_equal(dim(arr), c(4L, 2L, 3L, 3L))
    for (s in 1:4) for (m in 1:2) {
      L <- matrix(arr[s, m, , ], nrow = 3, ncol = 3)
      expect_true(all(L[upper.tri(L)] == 0))                    # lower-triangular
      Om <- L %*% t(L)
      expect_equal(diag(Om), rep(1, 3), tolerance = 1e-6)         # unit diagonal
      ev <- eigen(Om, symmetric = TRUE, only.values = TRUE)$values
      expect_true(all(ev > -1e-8))                                # positive semi-definite
    }
  }
})

test_that("D = 1 LKJ shortcut returns a trivial [1] correlation without invoking Stan", {
  arr <- .draw_lkj_cholesky_prior(D = 1L, eta = 2.0, n_matrices = 3L, n_states = 5L, seed = 1L)
  expect_equal(dim(arr), c(5L, 3L, 1L, 1L))
  expect_true(all(arr == 1))
})

test_that("simulate_probit_prior_predictive() is reproducible under a fixed seed", {
  skip_if_not(.priorppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .priorppc_make_fit()
  d1 <- simulate_probit_prior_predictive(fit, n_states = 15L, seed = 42L)
  d2 <- simulate_probit_prior_predictive(fit, n_states = 15L, seed = 42L)
  st1 <- d1$generate_state(3L)
  st2 <- d2$generate_state(3L)
  expect_identical(st1$Y_rep_complete, st2$Y_rep_complete)

  d3 <- simulate_probit_prior_predictive(fit, n_states = 15L, seed = 99L)
  st3 <- d3$generate_state(3L)
  expect_false(identical(st1$Y_rep_complete, st3$Y_rep_complete))
})

test_that("simulate_probit_prior_predictive() respects the observation mask and never conditions on observed AST", {
  skip_if_not(.priorppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .priorppc_make_fit()
  draws <- simulate_probit_prior_predictive(fit, n_states = 5L, seed = 1L)
  st <- draws$generate_state(1L)
  expect_equal(unname(is.na(st$Y_rep)), unname(is.na(draws$setup$obs_ast_mat)))
  expect_false(anyNA(st$Y_rep_complete))
  expect_true(inherits(draws, "amr_ppc_draws"))
  expect_true(inherits(draws, "amr_prior_predictive_draws"))
})

test_that("compute_prior_predictive_status() flags an artificially extreme prior as implausible", {
  skip_if_not(.priorppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .priorppc_make_fit()
  extreme <- simulate_probit_prior_predictive(fit, n_states = 20L, seed = 3L,
                                               prior_config_override = list(beta_sd = 30))
  status <- compute_prior_predictive_status(extreme)
  expect_equal(status$status, "fail_implausible_prior_predictive")
  expect_true(status$summary$fraction_probability_lt_0.001 +
              status$summary$fraction_probability_gt_0.999 > 0.5)
})

test_that("compute_prior_predictive_status() on default (project) priors does not fail outright", {
  skip_if_not(.priorppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .priorppc_make_fit()
  reasonable <- simulate_probit_prior_predictive(fit, n_states = 30L, seed = 3L)
  status <- compute_prior_predictive_status(reasonable)
  expect_false(identical(status$status, "fail_implausible_prior_predictive"))
})

test_that("compute_prior_predictive_status() reports insufficient_prior_check for a non-draws input", {
  res <- compute_prior_predictive_status(NULL)
  expect_equal(res$status, "insufficient_prior_check")
  res2 <- compute_prior_predictive_status(list(not = "a draws object"))
  expect_equal(res2$status, "insufficient_prior_check")
})

test_that("compute_prior_predictive_fingerprint() is stable under no change and sensitive to prior/config changes", {
  skip_if_not(.priorppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .priorppc_make_fit()
  fp1 <- compute_prior_predictive_fingerprint(fit)
  fp2 <- compute_prior_predictive_fingerprint(fit)
  expect_identical(fp1$fingerprint, fp2$fingerprint)

  fit_modified <- fit
  fit_modified$prior_config_used$beta_sd <- fit$prior_config_used$beta_sd + 1
  fp3 <- compute_prior_predictive_fingerprint(fit_modified)
  expect_false(identical(fp1$fingerprint, fp3$fingerprint))

  fit_residual_changed <- fit
  fit_residual_changed$residual_structure <- "correlated"
  fp4 <- compute_prior_predictive_fingerprint(fit_residual_changed)
  expect_false(identical(fp1$fingerprint, fp4$fingerprint))
})

test_that("compute_probit_ppc_statistics() accepts prior predictive draws directly (shared amr_ppc_draws shape)", {
  skip_if_not(.priorppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .priorppc_make_fit()
  prior_draws <- simulate_probit_prior_predictive(fit, n_states = 10L, seed = 1L)
  stats_tbl <- compute_probit_ppc_statistics(prior_draws, statistics = c("marginal", "pairwise"))
  expect_true(is.data.frame(stats_tbl))
  expect_true(all(c("statistic_name", "stratum", "ppc_tail_probability") %in% names(stats_tbl)))
})
