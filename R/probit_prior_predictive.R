# probit_prior_predictive.R
#
# Prior predictive checking for the multivariate probit resistance-profile
# model. Answers: "what AMR datasets are implied by the PRIORS, before
# conditioning on any observed AST outcomes?" -- a distinct question from
# posterior predictive checking (probit_posterior_predictive.R, "can the
# FITTED model reproduce data resembling what was observed?"). See that
# file's header for the full four-status-field discipline this package
# maintains (diagnostic_status / posterior_predictive_status /
# profile_validation_status / prior_predictive_status, never merged).
#
# Methodological reference: Stan User's Guide, "Posterior and Prior
# Predictive Checks":
# https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html
#
#   theta_sim ~ p(theta)
#   y_sim     ~ p(y | theta_sim)
#
# No NUTS sampling is used here -- parameter states are drawn directly and
# independently from the declared priors (Stan fixed_param-style prior
# predictive simulation), respecting EXACTLY the same priors as the fitted
# model: beta ~ N(0, beta_sd); per random-effect block, tau_re[r] ~
# half-N(0, tau_sd), L_corr_re[r] ~ LKJCholesky(lkj_eta); residual L_Omega ~
# LKJCholesky(lkj_eta) when residual_structure == "correlated" (see
# .amr_probit_stan_generic_correlated()/.amr_probit_stan_generic_identity()
# in daly_resistance_profiles.R for the fitting model these mirror). LKJ
# correlation-matrix draws use Stan's own lkj_corr_cholesky_rng() via a
# minimal dedicated fixed_param program (compiled once, reusing cmdstanr's
# content-hash cache) -- NOT an ad-hoc R-side correlation-matrix sampler --
# guaranteeing bit-for-bit distributional consistency with the LKJ prior the
# fitting model itself uses.
#
# The purpose of prior predictive checking is NOT to make the simulated
# distribution match the observed data closely -- it is to flag obviously
# implausible/extreme implications of the CHOSEN priors before conditioning
# on any outcomes.

# ---------------------------------------------------------------------------
# Internal: LKJ correlation-matrix Cholesky-factor prior sampling
# ---------------------------------------------------------------------------

#' @keywords internal
.amr_probit_stan_prior_lkj <- function() {
  r"(
data {
  int<lower=1> D;
  int<lower=1> n_matrices;
  real<lower=1> lkj_eta;
}
parameters {
}
model {
}
generated quantities {
  array[n_matrices] cholesky_factor_corr[D] L;
  for (m in 1:n_matrices) {
    L[m] = lkj_corr_cholesky_rng(D, lkj_eta);
  }
}
)"
}

#' Draw \code{n_states * n_matrices} iid LKJ-distributed correlation-matrix
#' Cholesky factors via Stan's own \code{lkj_corr_cholesky_rng()}, in a
#' single \code{fixed_param} sample() call (no NUTS, no warmup).
#'
#' @return \code{array(dim = c(n_states, n_matrices, D, D))}.
#' @keywords internal
.draw_lkj_cholesky_prior <- function(D, eta, n_matrices, n_states, seed) {
  if (D == 1L) {
    # A 1x1 correlation "matrix" is trivially [1]; no need to invoke Stan.
    return(array(1, dim = c(n_states, n_matrices, 1L, 1L)))
  }
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("Package 'cmdstanr' is required for prior-predictive LKJ correlation-matrix sampling.",
         call. = FALSE)
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required (installed with cmdstanr).", call. = FALSE)

  stan_file <- cmdstanr::write_stan_file(.amr_probit_stan_prior_lkj())
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = list(D = as.integer(D), n_matrices = as.integer(n_matrices), lkj_eta = as.numeric(eta)),
    fixed_param = TRUE, chains = 1L, iter_warmup = 0L, iter_sampling = as.integer(n_states),
    seed = as.integer(seed), refresh = 0L
  )
  draws_mat <- posterior::as_draws_matrix(fit$draws(variables = "L"))

  arr <- array(0, dim = c(n_states, n_matrices, D, D))
  for (m in seq_len(n_matrices)) {
    for (a in seq_len(D)) for (b in seq_len(D)) {
      col <- sprintf("L[%d,%d,%d]", m, a, b)
      arr[, m, a, b] <- draws_mat[, col]
    }
  }
  arr
}

# ---------------------------------------------------------------------------
# Internal: prior-predictive setup and per-state parameter/data draw
# ---------------------------------------------------------------------------

#' Static (non-stochastic) per-event setup for prior predictive simulation --
#' mirrors the STATIC portion of \code{.probit_predictive_draws_setup()}
#' (posterior draw arrays are irrelevant here; a prior predictive draw needs
#' only the real design structure -- X_event, grouping indices, class
#' panel -- as conditioning predictors, per Part 10).
#' @keywords internal
.prior_predictive_draws_setup <- function(fitted_model, prior_config_override) {
  class_cols   <- fitted_model$class_cols
  event_meta   <- fitted_model$event_metadata
  upper_re_col <- fitted_model$upper_re_col
  pathogen_col <- fitted_model$pathogen_col
  re_prep      <- fitted_model$random_effects_prep
  residual_structure <- .null_default(fitted_model$residual_structure, "identity")

  D <- length(fitted_model$index_maps$class_levels)
  X_event_sim  <- if (!is.null(fitted_model$X_event)) fitted_model$X_event else fitted_model$X_design
  event_re_idx <- fitted_model$event_re_idx
  K <- ncol(X_event_sim)

  stopifnot(".event_idx" %in% names(event_meta))
  has_obs <- event_meta$.event_idx %in% fitted_model$data_long$ev_idx
  event_meta_obs <- event_meta[has_obs, , drop = FALSE]
  N_ev <- nrow(event_meta_obs)
  ev_row_idx <- event_meta_obs$.event_idx
  X_event <- X_event_sim[ev_row_idx, , drop = FALSE]
  flat_re_idx_obs <- event_re_idx[ev_row_idx, , drop = FALSE]

  obs_ast_mat <- as.matrix(event_meta_obs[, class_cols, drop = FALSE])
  storage.mode(obs_ast_mat) <- "double"
  obs_mask <- !is.na(obs_ast_mat)

  default_pc <- list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0)
  pc <- utils::modifyList(default_pc, .null_default(fitted_model$prior_config_used, list()))
  pc <- utils::modifyList(pc, .null_default(prior_config_override, list()))

  list(
    N_ev = N_ev, D = D, K = K, R = re_prep$R,
    n_levels = re_prep$n_levels, level_start = re_prep$level_start,
    total_re_levels = re_prep$total_re_levels,
    class_cols = class_cols, event_meta_obs = event_meta_obs,
    obs_ast_mat = obs_ast_mat, obs_mask = obs_mask,
    flat_re_idx_obs = flat_re_idx_obs,
    random_effects_prep = re_prep,
    upper_re_col = upper_re_col, pathogen_col = pathogen_col,
    residual_structure = residual_structure,
    eligibility_report = fitted_model$eligibility_report,
    X_event = X_event, prior_config = pc
  )
}

#' Draw one prior predictive parameter state and reconstruct mu, mirroring
#' the Stan fitting model's transformed-parameters block EXACTLY:
#' \code{re_effect[, lo:hi] = diag_pre_multiply(tau_re[r], L_corr_re[r]) * z_re[, lo:hi]},
#' \code{mu = X_event beta + re_contribution(re_effect, flat_re_idx)}.
#' @keywords internal
.prior_draw_state <- function(setup, lkj_arr, s) {
  K <- setup$K; D <- setup$D; R <- setup$R
  pc <- setup$prior_config

  beta <- matrix(stats::rnorm(K * D, mean = 0, sd = pc$beta_sd), nrow = K, ncol = D)

  re_effect <- matrix(0, nrow = D, ncol = setup$total_re_levels)
  for (r in seq_len(R)) {
    lo <- setup$level_start[r]; hi <- lo + setup$n_levels[r] - 1L
    tau_r <- abs(stats::rnorm(D, mean = 0, sd = pc$tau_sd))     # half-normal
    L_r <- matrix(lkj_arr[s, r, , ], nrow = D, ncol = D)
    z_r <- matrix(stats::rnorm(D * setup$n_levels[r]), nrow = D, ncol = setup$n_levels[r])
    re_effect[, lo:hi] <- (diag(tau_r, nrow = D) %*% L_r) %*% z_r
  }

  L_Omega <- if (identical(setup$residual_structure, "correlated")) {
    matrix(lkj_arr[s, R + 1L, , ], nrow = D, ncol = D)
  } else NULL

  mu <- (setup$X_event %*% beta) + re_contribution(re_effect, setup$flat_re_idx_obs)
  list(beta = beta, re_effect = re_effect, L_Omega = L_Omega, mu = mu)
}

# ---------------------------------------------------------------------------
# simulate_probit_prior_predictive()
# ---------------------------------------------------------------------------

#' Simulate prior predictive replicate datasets for a probit model
#'
#' Draws \code{n_states} independent parameter states directly from the
#' declared priors (no NUTS/no conditioning on outcomes) and generates a
#' complete replicate outcome matrix from each, using the REAL fixed-effect
#' design matrix and random-effect grouping structure of \code{fitted_model}
#' as conditioning predictors (per Stan User's Guide-style prior predictive
#' simulation). Reuses the identical generative mechanism as
#' \code{\link{simulate_probit_posterior_predictive}} (Bernoulli(Phi(mu)) for
#' identity residual, \code{Z ~ MVN(mu, Omega); Y = I(Z>0)} for correlated
#' residual) -- only the source of \code{theta} differs (prior vs. posterior
#' draws).
#'
#' Returns an object of class \code{c("amr_prior_predictive_draws",
#' "amr_ppc_draws")} -- structurally compatible with
#' \code{\link{compute_probit_ppc_statistics}} (which accepts any
#' \code{"amr_ppc_draws"}-classed object), so the same AMR-specific
#' discrepancy-statistic machinery can summarise the prior predictive
#' distribution's own properties (see
#' \code{\link{compute_prior_predictive_status}}, which additionally derives
#' the prior-specific plausibility summaries in Part 11 of the predictive-
#' checking specification: fraction of event-class probabilities near 0/1,
#' fraction of degenerate all-resistant/all-susceptible profiles, and
#' hospital-level spread).
#'
#' @param fitted_model A fitted model object as returned by
#'   \code{\link{fit_bayesian_multivariate_probit}} -- supplies the real
#'   \code{X_event} design matrix, random-effect grouping structure, class
#'   panel, and \code{prior_config_used} (the priors to replicate exactly).
#'   Only design/configuration fields are read; no posterior draws are used.
#' @param n_states Integer; number of independent prior states to draw.
#'   Default 1000.
#' @param seed Integer seed (deterministic).
#' @param prior_config_override Optional named list overriding any of
#'   \code{beta_sd}, \code{tau_sd}, \code{lkj_eta} from
#'   \code{fitted_model$prior_config_used} -- e.g. to prior-predictive-check
#'   a CANDIDATE prior before fitting.
#' @param preserve_observation_mask,return_replicates See
#'   \code{\link{simulate_probit_posterior_predictive}}.
#' @return An \code{"amr_prior_predictive_draws"}/\code{"amr_ppc_draws"}
#'   object with the same shape as
#'   \code{\link{simulate_probit_posterior_predictive}}'s return value.
#' @references
#' Stan Development Team. "Posterior and Prior Predictive Checks." Stan
#' User's Guide. \url{https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html}
#' @export
simulate_probit_prior_predictive <- function(
    fitted_model,
    n_states = 1000L,
    seed = 123L,
    prior_config_override = NULL,
    preserve_observation_mask = TRUE,
    return_replicates = FALSE
) {
  setup <- .prior_predictive_draws_setup(fitted_model, prior_config_override)
  S <- as.integer(n_states)
  if (is.na(S) || S < 1L) stop("`n_states` must resolve to a positive integer.", call. = FALSE)

  n_matrices <- setup$R + if (identical(setup$residual_structure, "correlated")) 1L else 0L
  lkj_arr <- .draw_lkj_cholesky_prior(setup$D, setup$prior_config$lkj_eta, n_matrices, S, seed)

  generate_state <- local({
    cache <- new.env(parent = emptyenv())
    function(s) {
      key <- as.character(as.integer(s))
      cached <- cache[[key]]
      if (!is.null(cached)) return(cached)

      val <- .ppc_with_local_seed(.ppc_state_seed(seed, s), {
        theta <- .prior_draw_state(setup, lkj_arr, s)
        Y_complete <- if (identical(setup$residual_structure, "correlated")) {
          .ppc_generate_correlated(theta$mu, theta$L_Omega)
        } else {
          .ppc_generate_identity(theta$mu)
        }
        Y_masked <- Y_complete
        Y_masked[!setup$obs_mask] <- NA_integer_
        list(Y_rep_complete = Y_complete, Y_rep = Y_masked, mu = theta$mu)
      })
      cache[[key]] <- val
      val
    }
  })

  out <- list(
    setup = setup,
    generate_state = generate_state,
    n_states_used = S,
    seed_used = as.integer(seed),
    preserve_observation_mask = isTRUE(preserve_observation_mask),
    residual_structure = setup$residual_structure,
    return_replicates = isTRUE(return_replicates),
    Y_rep_array = NULL,
    Y_rep_complete_array = NULL
  )

  if (isTRUE(return_replicates)) {
    n_elements <- as.double(S) * setup$N_ev * setup$D
    max_elements <- 5e7
    if (n_elements > max_elements) {
      stop(sprintf(
        paste(
          "return_replicates = TRUE would materialize %.0f array elements",
          "(n_states=%d x N_events=%d x D=%d), above the safety ceiling of %.0f.",
          "Reduce n_states, or leave return_replicates = FALSE."
        ), n_elements, S, setup$N_ev, setup$D, max_elements
      ), call. = FALSE)
    }
    Y_rep_array          <- array(NA_integer_, dim = c(S, setup$N_ev, setup$D))
    Y_rep_complete_array <- array(NA_integer_, dim = c(S, setup$N_ev, setup$D))
    for (s in seq_len(S)) {
      st <- generate_state(s)
      Y_rep_complete_array[s, , ] <- st$Y_rep_complete
      Y_rep_array[s, , ] <- if (isTRUE(preserve_observation_mask)) st$Y_rep else st$Y_rep_complete
    }
    out$Y_rep_array <- Y_rep_array
    out$Y_rep_complete_array <- Y_rep_complete_array
  }

  class(out) <- c("amr_prior_predictive_draws", "amr_ppc_draws")
  out
}

# ---------------------------------------------------------------------------
# compute_prior_predictive_status()
# ---------------------------------------------------------------------------

#' Classify plausibility of the prior predictive distribution
#'
#' Unlike posterior predictive checking, the goal here is NOT to make the
#' simulated distribution match the observed data closely -- it is to flag
#' obviously implausible/extreme implications of the CHOSEN priors, before
#' conditioning on any outcomes. Computes, per prior state and averaged
#' across states: the fraction of event-class probabilities
#' \eqn{\Phi(\mu_{ed})} below 0.001 or above 0.999 (near-deterministic
#' implied resistance), the fraction of generated events that are
#' degenerate (all classes resistant or all susceptible), and the spread of
#' per-"hospital" (the model's primary declared random-effect grouping)
#' mean implied probability -- a very large spread implies an implausibly
#' extreme facility-to-facility prior.
#'
#' @param prior_draws An \code{"amr_prior_predictive_draws"} object from
#'   \code{\link{simulate_probit_prior_predictive}}.
#' @param thresholds Named list overriding any of the defaults:
#'   \code{max_fraction_extreme_probability = 0.10},
#'   \code{severe_fraction_extreme_probability = 0.50},
#'   \code{max_fraction_degenerate_profiles = 0.50},
#'   \code{severe_fraction_degenerate_profiles = 0.90},
#'   \code{max_hospital_spread_sd = 0.35}. All thresholds are project
#'   decisions, not universal truths -- override freely and document why.
#' @return List with \code{status} (one of \code{"pass"},
#'   \code{"warning_extreme_prior_predictions"},
#'   \code{"warning_excessive_hospital_heterogeneity"},
#'   \code{"warning_degenerate_profiles"},
#'   \code{"fail_implausible_prior_predictive"},
#'   \code{"insufficient_prior_check"}), \code{reasons}, \code{thresholds_used},
#'   and \code{summary} (the averaged plausibility fractions).
#' @export
compute_prior_predictive_status <- function(prior_draws, thresholds = list()) {
  th_defaults <- list(
    max_fraction_extreme_probability = 0.10,
    severe_fraction_extreme_probability = 0.50,
    max_fraction_degenerate_profiles = 0.50,
    severe_fraction_degenerate_profiles = 0.90,
    max_hospital_spread_sd = 0.35
  )
  th <- utils::modifyList(th_defaults, thresholds)

  if (is.null(prior_draws) || !inherits(prior_draws, "amr_prior_predictive_draws"))
    return(list(status = "insufficient_prior_check", reasons = character(0),
                thresholds_used = th, summary = list()))

  setup <- prior_draws$setup
  S <- prior_draws$n_states_used
  if (is.null(S) || S < 1L)
    return(list(status = "insufficient_prior_check", reasons = character(0),
                thresholds_used = th, summary = list()))

  hospitals <- sort(unique(setup$event_meta_obs[[setup$upper_re_col]]))
  frac_lt <- numeric(S); frac_gt <- numeric(S)
  frac_all_r <- numeric(S); frac_all_s <- numeric(S)
  hosp_spread <- numeric(S)

  for (s in seq_len(S)) {
    st <- prior_draws$generate_state(s)
    p <- stats::pnorm(st$mu)
    p_tested <- p[setup$obs_mask]
    frac_lt[s] <- if (length(p_tested) > 0L) mean(p_tested < 0.001) else NA_real_
    frac_gt[s] <- if (length(p_tested) > 0L) mean(p_tested > 0.999) else NA_real_

    Yc <- st$Y_rep_complete
    frac_all_r[s] <- mean(rowSums(Yc == 1) == setup$D)
    frac_all_s[s] <- mean(rowSums(Yc == 0) == setup$D)

    if (length(hospitals) >= 2L) {
      hosp_means <- vapply(hospitals, function(h) {
        idx <- which(setup$event_meta_obs[[setup$upper_re_col]] == h)
        mean(p[idx, , drop = FALSE])
      }, numeric(1L))
      hosp_spread[s] <- stats::sd(hosp_means)
    } else {
      hosp_spread[s] <- NA_real_
    }
  }

  summary <- list(
    fraction_probability_lt_0.001 = mean(frac_lt, na.rm = TRUE),
    fraction_probability_gt_0.999 = mean(frac_gt, na.rm = TRUE),
    fraction_all_resistant = mean(frac_all_r, na.rm = TRUE),
    fraction_all_susceptible = mean(frac_all_s, na.rm = TRUE),
    hospital_spread_sd_mean = mean(hosp_spread, na.rm = TRUE)
  )

  extreme_frac <- summary$fraction_probability_lt_0.001 + summary$fraction_probability_gt_0.999
  degenerate_frac <- summary$fraction_all_resistant + summary$fraction_all_susceptible

  reasons <- character(0)
  if (!is.na(extreme_frac) && extreme_frac > th$severe_fraction_extreme_probability)
    reasons <- c(reasons, "fail_implausible_prior_predictive", "warning_extreme_prior_predictions")
  else if (!is.na(extreme_frac) && extreme_frac > th$max_fraction_extreme_probability)
    reasons <- c(reasons, "warning_extreme_prior_predictions")

  if (!is.na(degenerate_frac) && degenerate_frac > th$severe_fraction_degenerate_profiles)
    reasons <- c(reasons, "fail_implausible_prior_predictive", "warning_degenerate_profiles")
  else if (!is.na(degenerate_frac) && degenerate_frac > th$max_fraction_degenerate_profiles)
    reasons <- c(reasons, "warning_degenerate_profiles")

  if (!is.na(summary$hospital_spread_sd_mean) && summary$hospital_spread_sd_mean > th$max_hospital_spread_sd)
    reasons <- c(reasons, "warning_excessive_hospital_heterogeneity")

  reasons <- unique(reasons)
  status <- if (any(grepl("^fail_", reasons))) {
    "fail_implausible_prior_predictive"
  } else if (any(grepl("^warning_", reasons))) {
    reasons[grepl("^warning_", reasons)][1]
  } else {
    "pass"
  }

  list(status = status, reasons = reasons, thresholds_used = th, summary = summary)
}

# ---------------------------------------------------------------------------
# compute_prior_predictive_fingerprint() -- Part 12: reuse/certification
# ---------------------------------------------------------------------------

#' Fingerprint the prior-generative configuration of a fitted model
#'
#' Two experiments that share the same prior-generative structure (same
#' residual structure, same beta/tau/LKJ priors, same declared random-effect
#' block structure, same fixed-effect design-column definitions, same
#' \code{D}) can safely REUSE a prior predictive check's results rather than
#' rerunning \code{\link{simulate_probit_prior_predictive}} for every
#' experiment -- prior predictive checking depends only on the generative
#' structure, never on the fitted posterior. If ANY component changes
#' (a prior, a scaling, the declared block structure), the fingerprint
#' changes and any cached result must be invalidated.
#'
#' @param fitted_model A fitted model object as returned by
#'   \code{\link{fit_bayesian_multivariate_probit}}.
#' @return List with \code{fingerprint} (a single character hash via
#'   \code{rlang::hash()}, already an \pkg{anumaan} dependency -- no new
#'   package dependency introduced) and \code{components} (the exact named
#'   list that was hashed, for audit/debugging).
#' @export
compute_prior_predictive_fingerprint <- function(fitted_model) {
  re_prep <- fitted_model$random_effects_prep
  X_event <- if (!is.null(fitted_model$X_event)) fitted_model$X_event else fitted_model$X_design

  components <- list(
    residual_structure = .null_default(fitted_model$residual_structure, "identity"),
    prior_config = .null_default(fitted_model$prior_config_used, list()),
    random_effect_blocks = list(
      block_names = re_prep$block_names,
      group_cols  = re_prep$group_cols,
      n_levels    = as.integer(re_prep$n_levels)
    ),
    fixed_effect_columns = colnames(X_event),
    fixed_effect_scaling = list(
      mean = round(colMeans(X_event), 6L),
      sd   = round(apply(X_event, 2L, stats::sd), 6L)
    ),
    D = length(fitted_model$class_cols)
  )

  list(fingerprint = rlang::hash(components), components = components)
}
