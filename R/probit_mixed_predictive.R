# probit_mixed_predictive.R
#
# Mixed predictive checking (Stan User's Guide terminology) for the
# multivariate probit resistance-profile model: for a requested random-
# effect block, retain the FITTED posterior hyperparameters (tau_r, the
# block's own cross-class correlation R_block[r]) but draw ENTIRELY NEW
# random-effect levels rather than reusing the fitted ones -- answering
# "what would NEW groups drawn from the inferred group-level distribution
# look like?" (e.g. "does the inferred hospital population generate
# realistic new hospitals?"), a distinct question from standard posterior
# predictive checking (probit_posterior_predictive.R), which conditions on
# the fitted random effects exactly as estimated for the OBSERVED groups.
#
# https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html
#
# SCOPE (per the predictive-checking specification, Part 13): the API is
# fully designed and the core generative step is fully implemented and
# tested (one smoke test) for the explicitly scoped hospital-level-block use
# case. The full Part 5/6-style discrepancy-statistic battery, performance
# benchmarking, and synthetic-recovery testing for mixed predictive checks
# are explicitly DEFERRED to a later phase -- not silently dropped.

#' Simulate mixed predictive replicate datasets: new random-effect levels
#' for a requested block, retaining fitted hyperparameters
#'
#' For each requested block \eqn{r} in \code{blocks_to_replicate}, retains
#' the posterior draws of \eqn{\tau_r^{(s)}} and the block's own fitted
#' cross-class correlation \eqn{R_{block[r]}^{(s)}} (a generated quantity
#' already stored in \code{fitted_model$draws}), but draws
#' \code{n_new_levels_per_block} ENTIRELY NEW synthetic level effects
#' \deqn{u_{r,\text{new}}^{(s)} \sim MVN(0, \Sigma_r^{(s)})} per posterior
#' state, where \eqn{\Sigma_r^{(s)} = \mathrm{diag}(\tau_r^{(s)})\,
#' R_{block[r]}^{(s)}\,\mathrm{diag}(\tau_r^{(s)})}. Every real fitted event
#' is deterministically (seeded) reassigned to one of the
#' \code{n_new_levels_per_block} new synthetic levels for each replicated
#' block, so a complete replicate outcome matrix can be generated exactly as
#' in \code{\link{simulate_probit_posterior_predictive}}. Blocks NOT in
#' \code{blocks_to_replicate} use their FITTED posterior random effects
#' (\code{unreplicated_block_behavior = "retain_fitted"}, default) or cause
#' an error (\code{"error"}) if any declared block is left unhandled.
#'
#' This is NOT standard posterior predictive checking (which conditions on
#' the observed groups' own fitted random effects) -- do not conflate the
#' two. \code{blocks_to_replicate} is currently exercised for a single
#' block at a time in this package's own tests (the hospital-level use case
#' explicitly scoped by Part 13 of the predictive-checking specification);
#' multi-block replication is implemented generically but not yet covered
#' by synthetic recovery tests.
#'
#' @param fitted_model A fitted model object as returned by
#'   \code{\link{fit_bayesian_multivariate_probit}}.
#' @param blocks_to_replicate Character vector; a subset of
#'   \code{fitted_model$random_effects_prep$block_names} to draw NEW levels
#'   for.
#' @param n_states Integer; number of posterior states to use. Default 500.
#' @param seed Integer seed (deterministic).
#' @param n_new_levels_per_block Integer; how many new synthetic levels to
#'   draw per replicated block. Default 1.
#' @param unreplicated_block_behavior \code{"retain_fitted"} (default) uses
#'   fitted posterior random effects for every declared block not in
#'   \code{blocks_to_replicate}; \code{"error"} requires every declared
#'   block to be listed in \code{blocks_to_replicate}.
#' @return An object of class \code{c("amr_mixed_predictive_draws",
#'   "amr_ppc_draws")} with \code{setup}, \code{generate_state}, and the
#'   other fields shared with
#'   \code{\link{simulate_probit_posterior_predictive}}'s return value, plus
#'   \code{blocks_to_replicate}, \code{unreplicated_blocks},
#'   \code{n_new_levels_per_block}, \code{unreplicated_block_behavior}.
#' @references
#' Stan Development Team. "Posterior and Prior Predictive Checks." Stan
#' User's Guide. \url{https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html}
#' @export
simulate_probit_mixed_predictive <- function(
    fitted_model,
    blocks_to_replicate,
    n_states = 500L,
    seed = 123L,
    n_new_levels_per_block = 1L,
    unreplicated_block_behavior = c("retain_fitted", "error")
) {
  unreplicated_block_behavior <- match.arg(unreplicated_block_behavior)
  re_prep <- fitted_model$random_effects_prep

  if (is.null(re_prep) || re_prep$R == 0L) {
    stop("No fitted random-effect blocks are available to replicate. Use simulate_probit_posterior_predictive() for fixed-only models.",
         call. = FALSE)
  }

  if (length(blocks_to_replicate) == 0L)
    stop("`blocks_to_replicate` must be non-empty.", call. = FALSE)
  bad_blocks <- setdiff(blocks_to_replicate, re_prep$block_names)
  if (length(bad_blocks) > 0L)
    stop("blocks_to_replicate must be a subset of fitted_model$random_effects_prep$block_names. ",
         "Not found: ", paste(bad_blocks, collapse = ", "), call. = FALSE)

  unreplicated <- setdiff(re_prep$block_names, blocks_to_replicate)
  if (identical(unreplicated_block_behavior, "error") && length(unreplicated) > 0L)
    stop("unreplicated_block_behavior = 'error' but the following declared blocks are not in ",
         "blocks_to_replicate: ", paste(unreplicated, collapse = ", "), call. = FALSE)

  n_new <- as.integer(n_new_levels_per_block)
  if (is.na(n_new) || n_new < 1L) stop("`n_new_levels_per_block` must be a positive integer.", call. = FALSE)

  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required (installed with cmdstanr).", call. = FALSE)

  set.seed(as.integer(seed))

  draws        <- fitted_model$draws
  class_cols   <- fitted_model$class_cols
  D <- length(fitted_model$index_maps$class_levels)
  X_event_sim  <- if (!is.null(fitted_model$X_event)) fitted_model$X_event else fitted_model$X_design
  event_re_idx <- fitted_model$event_re_idx
  K <- ncol(X_event_sim)
  residual_structure <- .null_default(fitted_model$residual_structure, "identity")
  R <- re_prep$R

  draws_mat <- posterior::as_draws_matrix(draws)
  n_total   <- nrow(draws_mat)
  S <- min(as.integer(n_states), n_total)
  if (is.na(S) || S < 1L) stop("`n_states` must resolve to a positive integer.", call. = FALSE)
  draw_idx  <- if (S < n_total) sort(sample.int(n_total, S)) else seq_len(n_total)
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  .arr2 <- function(prefix, d1, d2) {
    cols <- as.vector(outer(seq_len(d1), seq_len(d2), function(a, b) sprintf("%s[%d,%d]", prefix, a, b)))
    array(draws_mat[, cols, drop = FALSE], dim = c(S, d1, d2))
  }
  .arr3 <- function(prefix, d1, d2, d3) {
    arr <- array(0, dim = c(S, d1, d2, d3))
    for (a in seq_len(d1)) for (b in seq_len(d2)) for (cc in seq_len(d3)) {
      arr[, a, b, cc] <- draws_mat[, sprintf("%s[%d,%d,%d]", prefix, a, b, cc)]
    }
    arr
  }

  beta_arr    <- .arr2("beta", K, D)
  re_eff_arr  <- .arr2("re_effect", D, re_prep$total_re_levels)
  tau_re_arr  <- .arr2("tau_re", R, D)
  R_block_arr <- .arr3("R_block", R, D, D)
  L_omega_arr <- if (identical(residual_structure, "correlated")) .arr2("L_Omega", D, D) else NULL

  event_meta <- fitted_model$event_metadata
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

  replicate_idx <- match(blocks_to_replicate, re_prep$block_names)
  unrep_idx     <- match(unreplicated, re_prep$block_names)

  # Deterministic (seeded) round-robin assignment of every real event to one
  # of the n_new synthetic new levels, per replicated block.
  new_level_of_event <- lapply(seq_along(replicate_idx), function(i) {
    .ppc_with_local_seed(.ppc_state_seed(seed, 900000L + replicate_idx[i]), {
      ((sample.int(N_ev) - 1L) %% n_new) + 1L
    })
  })

  generate_state <- local({
    cache <- new.env(parent = emptyenv())
    function(s) {
      key <- as.character(as.integer(s))
      cached <- cache[[key]]
      if (!is.null(cached)) return(cached)

      val <- .ppc_with_local_seed(.ppc_state_seed(seed, s), {
        beta_s   <- matrix(beta_arr[s, , ],   nrow = K, ncol = D)
        re_eff_s <- matrix(re_eff_arr[s, , ], nrow = D, ncol = re_prep$total_re_levels)
        mu <- X_event %*% beta_s

        if (length(unrep_idx) > 0L) {
          mu <- mu + re_contribution(re_eff_s, flat_re_idx_obs[, unrep_idx, drop = FALSE])
        }

        for (i in seq_along(replicate_idx)) {
          r <- replicate_idx[i]
          tau_r <- tau_re_arr[s, r, ]
          Om_r  <- matrix(R_block_arr[s, r, , ], nrow = D, ncol = D)
          Om_r  <- (Om_r + t(Om_r)) / 2
          L_r <- tryCatch(t(chol(Om_r)), error = function(e) t(chol(Om_r + diag(1e-6, D))))
          z_new <- matrix(stats::rnorm(D * n_new), nrow = D, ncol = n_new)
          u_new <- (diag(tau_r, nrow = D) %*% L_r) %*% z_new     # D x n_new
          event_level <- new_level_of_event[[i]]                 # N_ev-length in 1..n_new
          mu <- mu + t(u_new[, event_level, drop = FALSE])
        }

        Y_complete <- if (identical(residual_structure, "correlated")) {
          .ppc_generate_correlated(mu, matrix(L_omega_arr[s, , ], nrow = D, ncol = D))
        } else {
          .ppc_generate_identity(mu)
        }
        Y_masked <- Y_complete
        Y_masked[!obs_mask] <- NA_integer_
        list(Y_rep_complete = Y_complete, Y_rep = Y_masked, mu = mu)
      })
      cache[[key]] <- val
      val
    }
  })

  setup <- list(
    N_ev = N_ev, D = D, class_cols = class_cols, event_meta_obs = event_meta_obs,
    obs_ast_mat = obs_ast_mat, obs_mask = obs_mask,
    upper_re_col = fitted_model$upper_re_col, pathogen_col = fitted_model$pathogen_col,
    residual_structure = residual_structure, eligibility_report = fitted_model$eligibility_report,
    random_effects_prep = re_prep, flat_re_idx_obs = flat_re_idx_obs
  )

  out <- list(
    setup = setup, generate_state = generate_state,
    n_states_used = S, seed_used = as.integer(seed),
    blocks_to_replicate = blocks_to_replicate, unreplicated_blocks = unreplicated,
    n_new_levels_per_block = n_new, unreplicated_block_behavior = unreplicated_block_behavior,
    residual_structure = residual_structure, preserve_observation_mask = TRUE,
    return_replicates = FALSE, Y_rep_array = NULL, Y_rep_complete_array = NULL
  )
  class(out) <- c("amr_mixed_predictive_draws", "amr_ppc_draws")
  out
}
