# Tests for the Stage-1 generic random-effect architecture
# (R/random_effects_generic.R): prepare_random_effects(), re_contribution(),
# and the normalize/validate helpers. No Stan fitting here -- see
# test-daly-resistance-profiles-bayesian.R for end-to-end fit tests using
# these structures.

test_that("legacy character vector and named-block spec produce the same structure", {
  d <- data.frame(hosp = rep(c("H1", "H2", "H3"), each = 4), adm = paste0("A", 1:12))
  re_legacy <- prepare_random_effects(d, c("hosp", "adm"))
  re_named  <- prepare_random_effects(d, list(
    list(name = "hospital", group_col = "hosp", terms = "intercept"),
    list(name = "admission", group_col = "adm", terms = "intercept")
  ))
  expect_equal(re_legacy$R, re_named$R)
  expect_identical(re_legacy$group_cols, re_named$group_cols)
  expect_identical(re_legacy$block_names, c("hosp", "adm"))
  expect_identical(re_named$block_names, c("hospital", "admission"))
})

test_that("an empty block list is a canonical zero-random-effect representation", {
  d <- data.frame(profile_group = c("H1", "H2", "H1"))
  re <- prepare_random_effects(d, list())

  expect_identical(re$R, 0L)
  expect_identical(re$total_re_levels, 0L)
  expect_identical(dim(re$group_index), c(3L, 0L))
  expect_identical(dim(re$flat_group_index), c(3L, 0L))
  expect_length(re$blocks, 0L)
})

test_that("anonymous numbered blocks are rejected", {
  expect_error(
    .normalize_random_effects_spec(list(list(group_col = "hosp"))),
    "missing a .name."
  )
})

test_that("random slopes are rejected in Stage 1", {
  expect_error(
    .normalize_random_effects_spec(list(list(name = "hospital", group_col = "hosp", terms = c("intercept", "slope_x")))),
    "random slopes"
  )
})

test_that("duplicate block names are rejected", {
  d <- data.frame(hosp = c("H1", "H2"), pat = c("P1", "P2"))
  expect_error(
    prepare_random_effects(d, list(list(name = "x", group_col = "hosp"), list(name = "x", group_col = "pat"))),
    "unique"
  )
})

test_that("missing/empty group ids are rejected", {
  d <- data.frame(hosp = c("H1", "", NA))
  expect_error(prepare_random_effects(d, "hosp"), "missing/empty")
})

test_that("nesting detection is data-driven, not order-dependent", {
  # Declared order: hospital -> admission -> patient, but patient actually
  # has multiple admissions -- admission must be detected as nested WITHIN
  # patient regardless of declaration order.
  d <- data.frame(
    hosp = "H1",
    patient = rep(c("P1", "P2"), each = 4),
    admission = c("P1_A1", "P1_A1", "P1_A2", "P1_A2", "P2_A1", "P2_A1", "P2_A2", "P2_A2")
  )
  re <- prepare_random_effects(d, list(
    list(name = "hospital", group_col = "hosp", terms = "intercept"),
    list(name = "admission", group_col = "admission", terms = "intercept"),
    list(name = "patient", group_col = "patient", terms = "intercept")
  ))
  rel <- re$pairwise_relationships
  adm_pat <- rel[rel$child_block == "admission" & rel$parent_block == "patient", ]
  expect_equal(adm_pat$relationship, "nested")
})

test_that("pairwise_relationships correctly identifies identical_partition", {
  # Two blocks that partition events identically (1:1 correspondence)
  d <- data.frame(a = rep(c("x1", "x2", "x3"), each = 3), b = rep(c("y1", "y2", "y3"), each = 3))
  re <- prepare_random_effects(d, c("a", "b"))
  rel <- re$pairwise_relationships
  ab <- rel[rel$child_block == "a" & rel$parent_block == "b", ]
  expect_equal(ab$relationship, "identical_partition")
})

test_that("min_repeated_levels is checked independently of singleton_threshold (regression test for the nesting bug)", {
  # 10 levels: 8 singleton, 2 repeated -> singleton_fraction = 80%, BELOW the
  # default 90% threshold. Before the fix, min_repeated_levels was only ever
  # checked INSIDE the singleton_fraction > 0.9 branch, so this case would
  # never trigger the min_repeated_levels warning at all.
  d <- data.frame(
    hosp = rep(c("H1", "H2", "H3"), each = 4),
    grp = c(paste0("g", 1:8), "g9", "g9", "g10", "g10")
  )
  warnings_seen <- character(0)
  withCallingHandlers(
    prepare_random_effects(d, c("hosp", "grp"), min_repeated_levels = 5L),
    warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  grp_warning <- warnings_seen[grepl("grp", warnings_seen)]
  expect_length(grp_warning, 1L)
  expect_match(grp_warning, "only 2 of its 10 levels have >=2 observations \\(need >= 5\\)")
})

test_that("singleton_threshold is configurable", {
  d <- data.frame(hosp = rep(c("H1", "H2", "H3"), each = 4), grp = paste0("g", 1:12))
  w <- character(0)
  withCallingHandlers(
    prepare_random_effects(d, c("hosp", "grp"), singleton_threshold = 0.5),
    warning = function(w2) { w <<- c(w, conditionMessage(w2)); invokeRestart("muffleWarning") }
  )
  grp_w <- w[grepl("grp", w)]
  expect_match(grp_w, "above the 50% threshold")
})

test_that("re_contribution() sums the correct blocks per event", {
  D <- 2; total_levels <- 5
  re_effect <- matrix(1:(D * total_levels), nrow = D)
  fgi <- matrix(c(1, 3, 2, 4, 5, 1), ncol = 2, byrow = TRUE)  # 3 events x 2 blocks
  out <- re_contribution(re_effect, fgi)
  expect_equal(out[1, ], re_effect[, 1] + re_effect[, 3])
  expect_equal(out[2, ], re_effect[, 2] + re_effect[, 4])
  expect_equal(out[3, ], re_effect[, 5] + re_effect[, 1])
})

test_that("re_contribution() handles a single event vector input", {
  D <- 2; total_levels <- 4
  re_effect <- matrix(1:(D * total_levels), nrow = D)
  out <- re_contribution(re_effect, c(1L, 3L))
  expect_equal(out, re_effect[, 1] + re_effect[, 3])
})

test_that("re_contribution() is an N by D zero matrix for zero blocks", {
  re_effect <- matrix(numeric(), nrow = 2L, ncol = 0L)
  flat_group_index <- matrix(integer(), nrow = 4L, ncol = 0L)
  out <- re_contribution(re_effect, flat_group_index)
  expect_identical(dim(out), c(4L, 2L))
  expect_true(all(out == 0))
})

test_that("print.amr_random_effects does not error", {
  d <- data.frame(hosp = c("H1", "H2"), adm = c("A1", "A2"))
  re <- prepare_random_effects(d, c("hosp", "adm"))
  expect_output(print(re), "amr_random_effects")
})
