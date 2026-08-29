# Tests for R/plot_display_labels.R (generic display-label helpers used by
# both the PPC report and, eventually, the validation report) and the
# structural properties of the Phase A redesigned PPC plotting functions in
# R/probit_posterior_predictive.R.

test_that("hospital_display_label() replaces underscores with spaces, nothing else", {
  expect_identical(hospital_display_label("AIIMS_trauma_center"), "AIIMS trauma center")
  expect_identical(hospital_display_label("Assam_Medical_College"), "Assam Medical College")
  expect_identical(hospital_display_label("SGRH"), "SGRH")
})

test_that("class_display_label() looks up known classes and falls back gracefully for unknown ones", {
  expect_identical(class_display_label("Third_generation_cephalosporins"), "Third-generation cephalosporins")
  expect_identical(class_display_label("Aminoglycosides"), "Aminoglycosides")
  # Unmapped class: falls back to underscore-to-space, never errors.
  expect_identical(class_display_label("Some_Future_Class"), "Some Future Class")
})

test_that("class_short_label() looks up known abbreviations and falls back gracefully", {
  expect_identical(class_short_label("Third_generation_cephalosporins"), "3GC")
  expect_identical(class_short_label("Aminoglycosides"), "AMG")
  expect_identical(class_short_label("Carbapenems"), "CARB")
  expect_identical(class_short_label("Fluoroquinolones"), "FQ")
  expect_identical(class_short_label("Penicillins"), "PEN")
  # Unmapped: 4-char upper-cased fallback, never errors.
  expect_identical(class_short_label("Zzznotarealclass"), "ZZZN")
})

test_that("class_pair_label() combines two classes with 'x' using short labels by default", {
  expect_identical(class_pair_label("Aminoglycosides", "Carbapenems"), "AMG x CARB")
  expect_identical(
    class_pair_label("Aminoglycosides", "Carbapenems", short = FALSE),
    "Aminoglycosides x Carbapenems"
  )
})

# --- Structural properties of the Phase A redesigned PPC plots -------------

.fake_ppc_marginal <- data.frame(
  statistic_name = "marginal_resistance",
  stratum = c("hospital:AIIMS_trauma_center|class:Carbapenems", "hospital:SGRH|class:Carbapenems"),
  observed_value = c(0.6, 0.7), replicated_mean = c(0.55, 0.65),
  replicated_q025 = c(0.4, 0.5), replicated_q975 = c(0.7, 0.8),
  support_status = "supported", stringsAsFactors = FALSE
)

test_that("marginal PPC plot facets by class and shows a proper resistance-proportion x-axis, never 'Value'", {
  p <- .ppc_plot_marginal(.fake_ppc_marginal, "test")
  expect_false(is.null(p))
  expect_identical(p$labels$x, "Resistance proportion")
  expect_false(identical(p$labels$x, "Value"))
  expect_identical(p$labels$y, "Hospital")
  # Faceting variable is the parsed class column, not the raw stratum string.
  expect_true("class_display" %in% names(p$data))
  expect_true(all(c("AIIMS trauma center", "SGRH") %in% p$data$hospital_display))
})

.fake_ppc_pairwise <- data.frame(
  statistic_name = c("pairwise_RR", "pairwise_RS"),
  stratum = c("hospital:AIIMS_trauma_center|Aminoglycosides_x_Carbapenems",
              "hospital:AIIMS_trauma_center|Aminoglycosides_x_Carbapenems"),
  observed_value = c(0.3, 0.1), replicated_mean = c(0.28, 0.12),
  replicated_q025 = c(0.2, 0.05), replicated_q975 = c(0.35, 0.2),
  support_status = "supported", stringsAsFactors = FALSE
)

test_that("pairwise PPC plot returns one page per RR/RS/SR/SS cell, faceted by class pair with short labels", {
  plots <- .ppc_plot_pairwise(.fake_ppc_pairwise, "test")
  expect_true(all(c("RR", "RS") %in% names(plots)))
  expect_true(is.null(plots$SR)) # no SR/SS rows in this fixture
  expect_true("AMG x CARB" %in% plots$RR$data$class_pair)
  expect_identical(plots$RR$labels$x, "Observed/replicated RR proportion")
  expect_false(grepl("hospital:|_x_", plots$RR$labels$x))
})

.fake_ppc_resistant_count <- data.frame(
  statistic_name = "resistant_count_proportion",
  stratum = c("hospital:AFMC|pathogen:bug|panel_size:2|C=0", "hospital:AFMC|pathogen:bug|panel_size:2|C=1"),
  observed_value = c(0.4, 0.6), replicated_mean = c(0.45, 0.55),
  replicated_q025 = c(0.3, 0.4), replicated_q975 = c(0.6, 0.7),
  support_status = "supported", stringsAsFactors = FALSE
)

test_that("resistant-count distribution plot never conflates different panel widths into one facet", {
  p <- .ppc_plot_resistant_count_distribution(.fake_ppc_resistant_count, "test")
  expect_false(is.null(p))
  expect_true(grepl("panel size = 2", p$data$facet_label[1]))
  expect_identical(p$labels$x, "Number of resistant classes (C)")
  expect_identical(p$labels$y, "Proportion of fully-observed events")
})

.fake_ppc_hh <- data.frame(
  statistic_name = c("hospital_heterogeneity_sd", "hospital_heterogeneity_range"),
  stratum = c("class:Carbapenems", "class:Carbapenems"),
  observed_value = c(0.1, 0.3), replicated_mean = c(0.12, 0.28),
  replicated_q025 = c(0.05, 0.2), replicated_q975 = c(0.15, 0.35),
  support_status = "supported", stringsAsFactors = FALSE
)

test_that("hospital heterogeneity plot facets by metric and never shows 'Value' as the axis label", {
  p <- .ppc_plot_hospital_heterogeneity(.fake_ppc_hh, "test")
  expect_false(is.null(p))
  expect_false(identical(p$labels$x, "Value"))
  expect_true(all(c("Standard deviation", "Range") %in% levels(p$data$metric_display)))
})

.fake_ppc_profile <- data.frame(
  statistic_name = c("profile_all_resistant_frequency", "profile_n_distinct", "profile_shannon_entropy"),
  stratum = rep("hospital:AFMC|pathogen:bug", 3),
  observed_value = c(0.5, 8, 1.2), replicated_mean = c(0.45, 7.5, 1.1),
  replicated_q025 = c(0.3, 5, 0.9), replicated_q975 = c(0.6, 10, 1.4),
  support_status = "supported", stringsAsFactors = FALSE
)

test_that("complete-profile statistics are split by unit, never mixed on one 'Value' axis", {
  plots <- .ppc_plot_complete_profile_statistics(.fake_ppc_profile, "test")
  expect_true(all(c("probability", "n_distinct", "entropy") %in% names(plots)))
  expect_identical(plots$probability$labels$x, "Proportion / concentration")
  expect_identical(plots$n_distinct$labels$x, "Number of distinct resistance profiles")
  expect_identical(plots$entropy$labels$x, "Profile entropy (nats)")
  # None of the three should claim the "Value" placeholder axis.
  expect_false(any(vapply(plots, function(p) identical(p$labels$x, "Value"), logical(1))))
  # The panel caveat about "all modelled classes resistant" != "pan-resistant"
  # must be present as a caption, not silently omitted.
  expect_true(grepl("pan-resistant", plots$probability$labels$caption, ignore.case = TRUE))
})
