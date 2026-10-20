# Pathway 1 (convex QP) resistance-profile estimation

agg_marginals <- function() {
  m <- data.frame(
    pathogen = rep(c("K", "E"), each = 2),
    antibiotic_class = rep(c("A", "B"), 2),
    n_tested = 100L,
    n_resistant = c(30L, 50L, 20L, 40L)
  )
  m$marginal_resistance <- m$n_resistant / m$n_tested
  m
}

line_list <- function(n = 80, seed = 1) {
  set.seed(seed)
  data.frame(
    isolate_id = rep(sprintf("I%03d", seq_len(n)), each = 3),
    organism_name = "K",
    org_group = "Enterobacterales",
    antibiotic_class = rep(c("A", "B", "C"), n),
    antibiotic_value = sample(c("S", "R"), 3 * n, replace = TRUE)
  )
}

test_that("enumerate_binary_profiles covers all 2^n S/R combinations", {
  p <- enumerate_binary_profiles(c("A", "B", "C"))
  expect_equal(nrow(p), 8L)
  expect_setequal(p$profile_delta, c("SSS", "RSS", "SRS", "RRS", "SSR", "RSR", "SRR", "RRR"))
  expect_equal(colSums(p[, c("A", "B", "C")]), c(A = 4, B = 4, C = 4))
})

test_that("build_constraint_matrix encodes marginals and independence fallback pairs", {
  p <- enumerate_binary_profiles(c("A", "B"))
  cm <- build_constraint_matrix(p, c(A = 0.3, B = 0.5))
  expect_equal(cm$constraint_names, c("marg_A", "marg_B", "pair_A_B"))
  expect_equal(unname(cm$v), c(0.3, 0.5, 0.15))
  expect_equal(cm$fallback_pairs, "pair_A_B")
  # pair row is 1 only for the RR profile
  expect_equal(cm$M[3, ], as.numeric(p$profile_delta == "RR"))
})

test_that("validate_aggregate_inputs accepts well-formed marginals and rejects bad rates", {
  expect_no_error(suppressMessages(validate_aggregate_inputs(agg_marginals())))
  bad <- agg_marginals()
  bad$marginal_resistance[1] <- 1.5
  expect_error(suppressMessages(validate_aggregate_inputs(bad)))
})

test_that("estimate_profiles_convex recovers the exact 2-class solution with both solvers", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("osqp")
  skip_if_not_installed("Matrix")
  pm <- list(K = c("A", "B"), E = c("A", "B"))
  for (solver in c("osqp", "quadprog")) {
    out <- suppressMessages(estimate_profiles_convex(agg_marginals(), panel_map = pm, solver = solver))
    k <- out[out$pathogen == "K", ]
    expect_true(all(k$convergence_flag), info = solver)
    probs <- setNames(k$profile_probability, k$profile_delta)
    # independence fallback: P(RR) = 0.3 * 0.5
    expect_equal(unname(probs[c("SS", "RS", "SR", "RR")]), c(0.35, 0.15, 0.35, 0.15),
      tolerance = 1e-3, info = solver
    )
  }
})

test_that("compute_marginal_resistance -> pairwise -> profiles chain is internally consistent", {
  skip_if_not_installed("quadprog")
  ast <- line_list()
  m <- suppressMessages(compute_marginal_resistance(ast))
  expect_equal(nrow(m$marginal), 3L)
  expect_true(all(m$marginal$n_tested == 80L))
  expect_equal(m$marginal$marginal_resistance, m$marginal$n_resistant / m$marginal$n_tested)

  co <- suppressMessages(compute_pairwise_coresistance(m))
  prev <- co[["K"]]$prevalence
  expect_equal(dim(prev), c(3L, 3L))
  expect_true(isSymmetric(unname(prev)))

  rp <- suppressMessages(compute_resistance_profiles(m, co, exclude_near_zero = FALSE))
  probs <- rp[["K"]]$profiles$probability
  expect_equal(sum(probs), 1, tolerance = 1e-8)
  expect_true(all(probs >= 0))

  checks <- check_profile_constraints(rp,
    marginals = m$marginal,
    pathogen_col = "organism_name", class_col = "antibiotic_class",
    rate_col = "marginal_resistance", tolerance = 1e-3
  )
  expect_true(all(checks$pass[checks$constraint_type %in% c("nonneg", "sum_to_one")]))
})

test_that("compute_marginal_resistance drops classes below min_n_tested", {
  ast <- line_list(n = 20)
  m <- suppressMessages(compute_marginal_resistance(ast, min_n_tested = 30))
  expect_equal(nrow(m$marginal), 0L)
})

test_that("bootstrap_profiles_convex is reproducible and brackets the point estimate", {
  skip_if_not_installed("quadprog")
  a <- suppressMessages(bootstrap_profiles_convex(agg_marginals(), B = 15L, seed = 7L, pathogen_col = "pathogen"))
  b <- suppressMessages(bootstrap_profiles_convex(agg_marginals(), B = 15L, seed = 7L, pathogen_col = "pathogen"))
  expect_identical(a, b)
  k <- a[["K"]]
  expect_true(all(k$lower <= k$probability_median + 1e-12))
  expect_true(all(k$upper >= k$probability_median - 1e-12))
  expect_equal(k$convergence_rate, rep(1, nrow(k)))
})

test_that("bootstrap_profiles_convex leaves the caller's RNG state untouched", {
  skip_if_not_installed("quadprog")
  set.seed(123)
  before <- .Random.seed
  suppressMessages(bootstrap_profiles_convex(agg_marginals(), B = 3L, seed = 1L, pathogen_col = "pathogen"))
  expect_identical(.Random.seed, before)
})

test_that("validate_profile_inputs and preprocess_for_profiles run on long AST data", {
  ast_long <- data.frame(
    isolate_id = rep(sprintf("ISO%02d", 1:6), each = 2),
    patient_id = rep(sprintf("P%02d", 1:6), each = 2),
    pathogen = rep(c("Klebsiella pneumoniae", "Escherichia coli"), each = 6),
    antibiotic_name = rep(c("Meropenem", "Ceftriaxone"), 6),
    ast_value = rep(c("R", "S", "S", "R"), 3),
    date_of_culture = as.Date("2024-01-01") + rep(0:5, each = 2),
    specimen_type = "Blood", state = "Maharashtra", age = 40
  )
  rep <- suppressWarnings(suppressMessages(validate_profile_inputs(ast_long)))
  expect_type(rep, "list")

  # duplicated isolate x drug row with conflicting results: worst phenotype wins
  dup <- rbind(ast_long, transform(ast_long[2, ], ast_value = "R"))
  res <- suppressWarnings(suppressMessages(preprocess_for_profiles(dup,
    panel_map = list(
      "Klebsiella pneumoniae" = c("Carbapenems", "Third-generation-cephalosporins"),
      "Escherichia coli" = c("Carbapenems", "Third-generation-cephalosporins")
    )
  )))
  log <- res$preprocessing_log
  expect_equal(log$n_rows_out[log$step == "4b_dedup"], nrow(ast_long))
  expect_true(all(c("data_wide", "preprocessing_log") %in% names(res)))
})

wide_from_long <- function(n = 60, seed = 2) {
  set.seed(seed)
  ast_long <- data.frame(
    isolate_id = rep(sprintf("ISO%03d", seq_len(n)), each = 2),
    patient_id = rep(sprintf("P%03d", seq_len(n)), each = 2),
    pathogen = "Klebsiella pneumoniae",
    antibiotic_name = rep(c("Meropenem", "Ceftriaxone"), n),
    ast_value = sample(c("R", "S"), 2 * n, replace = TRUE),
    date_of_culture = as.Date("2024-01-01"),
    specimen_type = "Blood", state = "MH", age = 40
  )
  suppressWarnings(suppressMessages(preprocess_for_profiles(ast_long, panel_map = kp_panel)))
}
kp_panel <- list("Klebsiella pneumoniae" = c("Carbapenems", "Third-generation-cephalosporins"))

test_that("compute_marginals_from_data keeps classes whose names contain hyphens", {
  res <- wide_from_long()
  mg <- suppressMessages(compute_marginals_from_data(res$data_wide, res$col_map_resolved, kp_panel))
  expect_setequal(mg$antibiotic_class, kp_panel[[1]])
  expect_true(all(mg$n_tested == 60L))
})

test_that("compute_pairwise_from_data and estimate_profiles_convex consume the data-driven marginals", {
  skip_if_not_installed("quadprog")
  res <- wide_from_long()
  mg <- suppressMessages(compute_marginals_from_data(res$data_wide, res$col_map_resolved, kp_panel))
  pw <- suppressMessages(compute_pairwise_from_data(res$data_wide, mg, res$col_map_resolved, kp_panel))
  expect_equal(nrow(pw), 1L)
  expect_equal(pw$n_co_tested, 60L)

  out <- suppressMessages(estimate_profiles_convex(mg, pw, kp_panel))
  expect_true(all(out$convergence_flag))
  expect_equal(sum(out$profile_probability), 1, tolerance = 1e-8)
  # marginals are reproduced: P(R in class 1) = RS + RR
  probs <- setNames(out$profile_probability, out$profile_delta)
  expect_equal(unname(probs["RS"] + probs["RR"]),
    mg$marginal_resistance[mg$antibiotic_class == "Carbapenems"],
    tolerance = 1e-4
  )

  via_dispatch <- suppressMessages(estimate_resistance_profiles(mg,
    method = "convex",
    pairwise = pw, panel_map = kp_panel
  ))
  expect_equal(via_dispatch$profiles$profile_probability, out$profile_probability, tolerance = 1e-8)
  expect_equal(via_dispatch$diagnostics$pct_converged, 1)
})

test_that("select_resistance_class keeps the highest-RR resistant class per event", {
  d <- data.frame(
    event_id = rep(c("E1", "E2"), each = 3),
    antibiotic_class = rep(c("Carbapenems", "3GC", "FQ"), 2),
    antibiotic_value = c("R", "R", "S", "S", "R", "R"),
    rr_value = rep(c(1.5, 1.2, 1.1), 2)
  )
  out <- suppressMessages(select_resistance_class(d))
  expect_equal(nrow(out), 2L)
  expect_equal(out$antibiotic_class[out$event_id == "E1"], "Carbapenems")
  expect_equal(out$antibiotic_class[out$event_id == "E2"], "3GC")
})
