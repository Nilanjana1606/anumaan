# Tests for R/probit_mixed_predictive.R: simulate_probit_mixed_predictive().
# Per the predictive-checking specification (Part 13), mixed predictive
# checking is scoped here to a working, tested CORE generative step and
# full API/argument-validation coverage -- NOT the full Part 5/6 discrepancy-
# statistic battery, performance benchmarking, or synthetic-recovery testing
# (explicitly deferred to a later phase).

.mixedppc_cmdstan_available <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({ cmdstanr::cmdstan_path(); TRUE }, error = function(e) FALSE))
}

.mixedppc_make_fit <- function(seed = 1L) {
  set.seed(seed)
  n_hosp <- 4L; n_per_hosp <- 30L
  hosp_eff <- stats::rnorm(n_hosp, sd = 0.4)
  wide <- do.call(rbind, lapply(seq_len(n_hosp), function(h) {
    n <- n_per_hosp
    age <- stats::rnorm(n)
    mu1 <- 0.1 + 0.2 * age + hosp_eff[h]
    mu2 <- -0.1 + 0.1 * age + hosp_eff[h]
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

test_that("simulate_probit_mixed_predictive() validates blocks_to_replicate against declared blocks", {
  skip_if_not(.mixedppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .mixedppc_make_fit()
  expect_error(
    simulate_probit_mixed_predictive(fit, blocks_to_replicate = "not_a_real_block"),
    "block_names"
  )
  expect_error(simulate_probit_mixed_predictive(fit, blocks_to_replicate = character(0)), "non-empty")
})

test_that("unreplicated_block_behavior = 'error' rejects fits with un-listed declared blocks", {
  # Argument validation happens before any posterior-draws extraction, so a
  # minimal hand-built object with a real two-block random_effects_prep
  # (no cmdstanr fit needed) is sufficient to exercise both branches.
  re_data <- data.frame(hosp = c("H1", "H2"), adm = c("A1", "A2"))
  fake_fit <- list(random_effects_prep = prepare_random_effects(
    re_data, list(list(name = "hospital", group_col = "hosp", terms = "intercept"),
                  list(name = "admission", group_col = "adm", terms = "intercept"))
  ))

  expect_error(
    simulate_probit_mixed_predictive(fake_fit, blocks_to_replicate = "hospital",
                                      unreplicated_block_behavior = "error"),
    "admission"
  )
  # Replicating BOTH declared blocks leaves nothing unhandled -- validation
  # should pass under "error" too (it fails later, deeper in the function,
  # for unrelated reasons since fake_fit has no real draws -- we only assert
  # the error is NOT the unreplicated-block validation message).
  err <- tryCatch(
    simulate_probit_mixed_predictive(fake_fit, blocks_to_replicate = c("hospital", "admission"),
                                      unreplicated_block_behavior = "error"),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("unreplicated_block_behavior", err))
})

test_that("simulate_probit_mixed_predictive() core generative step: new levels differ from fitted levels", {
  skip_if_not(.mixedppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .mixedppc_make_fit()

  mixed_draws <- simulate_probit_mixed_predictive(
    fit, blocks_to_replicate = "center_name", n_states = 10L, seed = 7L, n_new_levels_per_block = 5L
  )
  expect_s3_class(mixed_draws, "amr_mixed_predictive_draws")
  expect_true(inherits(mixed_draws, "amr_ppc_draws"))
  expect_identical(mixed_draws$blocks_to_replicate, "center_name")
  expect_identical(mixed_draws$unreplicated_blocks, character(0))

  st <- mixed_draws$generate_state(1L)
  expect_false(anyNA(st$Y_rep_complete))
  expect_equal(dim(st$Y_rep_complete), c(mixed_draws$setup$N_ev, mixed_draws$setup$D))

  # The generated mu must NOT equal what the SAME posterior state's FITTED
  # (observed-hospital) random effects would produce -- proving genuine new-
  # group generation rather than accidental reuse of fitted levels.
  draws_mat <- posterior::as_draws_matrix(fit$draws)
  re_prep <- fit$random_effects_prep
  ev_idx_obs <- which(fit$event_metadata$.event_idx %in% fit$data_long$ev_idx)
  beta_1 <- matrix(draws_mat[1, grep("^beta\\[", colnames(draws_mat))], nrow = 2, ncol = 2)
  re_eff_1 <- matrix(draws_mat[1, grep("^re_effect\\[", colnames(draws_mat))],
                      nrow = 2, ncol = re_prep$total_re_levels)
  mu_fitted <- (fit$X_event[ev_idx_obs, , drop = FALSE] %*% beta_1) +
    re_contribution(re_eff_1, fit$event_re_idx[ev_idx_obs, , drop = FALSE])

  expect_gt(mean(abs(st$mu - mu_fitted)), 1e-6)
})

test_that("simulate_probit_mixed_predictive() is reproducible under a fixed seed", {
  skip_if_not(.mixedppc_cmdstan_available(), "cmdstanr/CmdStan not available")
  fit <- .mixedppc_make_fit()
  d1 <- simulate_probit_mixed_predictive(fit, blocks_to_replicate = "center_name",
                                          n_states = 8L, seed = 11L)
  d2 <- simulate_probit_mixed_predictive(fit, blocks_to_replicate = "center_name",
                                          n_states = 8L, seed = 11L)
  expect_identical(d1$generate_state(2L)$Y_rep_complete, d2$generate_state(2L)$Y_rep_complete)
})
