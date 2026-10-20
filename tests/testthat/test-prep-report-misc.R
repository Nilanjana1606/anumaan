# Pipeline orchestration, reporting, and remaining prep_* / plot helpers

raw_amr <- function() {
  data.frame(
    PID = c("pt_001", "pt_001", "pt_002", "pt_003"),
    Date.of.admission = c("2021-01-01", "2021-01-01", "2021-01-05", "2021-01-10"),
    CultureDate = c("2021-01-03", "2021-01-03", "2021-01-08", "2021-01-11"),
    Final.outcome = c("alive", "alive", "expired", "Discharged Alive"),
    Organism = c("E. coli", "E. coli", "Klebsiella pneumoniae", "Staphylococcus aureus"),
    Sample_type1_name = c("Blood", "Blood", "Urine", "Blood"),
    Antibiotic = c("Amikacin", "Ciprofloxacin", "Meropenem", "Vancomycin"),
    Result = c("Resistant", "Susceptible", "R", "S"),
    Gender = c("M", "M", "female", "F"),
    Age = c(34, 34, 61, 47)
  )
}

run_quiet <- function(...) suppressWarnings(suppressMessages(run_preprocess(...)))

test_that("run_preprocess standardises, enriches and derives in one call", {
  r <- run_quiet(raw_amr(), verbose = FALSE, generate_report = FALSE)
  expect_s3_class(r, "amr_result")
  expect_named(r, c("data", "config", "log", "report", "metadata"))
  expect_equal(nrow(r$data), 4L)
  expect_true(all(c("organism_normalized", "antibiotic_class", "Age_bin", "infection_type") %in% names(r$data)))
  expect_equal(r$data$organism_normalized[1], "escherichia coli") # normalised key is lower-case
  expect_equal(r$data$antibiotic_class[r$data$antibiotic_name == "Meropenem"], "Carbapenems")
})

test_that("run_preprocess honours partial phases", {
  r <- run_quiet(raw_amr(), phases = "standardize", verbose = FALSE, generate_report = FALSE)
  expect_false("Age_bin" %in% names(r$data))
  expect_true("organism_normalized" %in% names(r$data))
})

test_that("summary.amr_result prints without error", {
  r <- run_quiet(raw_amr(), verbose = FALSE, generate_report = FALSE)
  expect_output(summary(r))
})

test_that("preprocessing report summarises raw data and exports to every supported format", {
  raw <- raw_amr()
  r <- run_quiet(raw, verbose = FALSE, generate_report = FALSE)
  rep <- suppressMessages(generate_preprocessing_report(raw, r$data, r$log, r$config, include_plots = FALSE))
  expect_s3_class(rep, "amr_preprocessing_report")
  expect_equal(rep$raw_data_summary$n_rows, 4L)
  expect_equal(rep$raw_data_summary$n_cols, ncol(raw))
  expect_output(print(rep))

  for (ext in c("html", "txt", "rds")) {
    f <- tempfile(fileext = paste0(".", ext))
    suppressMessages(export_report(rep, f))
    expect_true(file.exists(f), info = ext)
    expect_gt(file.size(f), 100)
  }
  expect_identical(readRDS(local({
    f <- tempfile(fileext = ".rds")
    suppressMessages(export_report(rep, f))
    f
  })), rep)
  expect_error(export_report(rep, tempfile(fileext = ".pdf")), "LaTeX")
  expect_error(export_report(rep, tempfile(fileext = ".xyz")), "Unsupported")
})

test_that("report building blocks return the expected summaries", {
  raw <- raw_amr()
  r <- run_quiet(raw, verbose = FALSE, generate_report = FALSE)
  s <- summarize_raw_data(raw)
  expect_equal(s$column_names, names(raw))
  expect_equal(nrow(s$missing_data), ncol(raw))

  expect_type(summarize_column_mapping(raw, r$data, r$log$column_mapping), "list")
  expect_type(summarize_transformations(raw, r$data, r$log), "list")
  q <- summarize_data_quality(raw, r$data, include_plots = FALSE)
  expect_type(q, "list")
  st <- compute_processing_statistics(raw, r$data, r$log)
  expect_type(st, "list")
})

test_that("prep_assign_facility_type maps case-insensitively with a default", {
  d <- data.frame(center_name = c("AIIMS Delhi", "CMC vellore", "Unknown"))
  out <- suppressMessages(prep_assign_facility_type(d,
    mapping = c("aiims delhi" = "Tertiary", "CMC Vellore" = "Tertiary"), default = "Other"
  ))
  expect_equal(out$facility_type, c("Tertiary", "Tertiary", "Other"))
})

test_that("prep_classify_mortality applies the mortality window", {
  m <- data.frame(
    final_outcome = c("Died", "Died", "Discharged"),
    date_of_culture = as.Date("2024-01-01"),
    date_of_final_outcome = as.Date(c("2024-01-05", "2024-02-20", "2024-01-10"))
  )
  out <- suppressMessages(prep_classify_mortality(m, window = 14))
  expect_equal(out$mortality_infection, c("Yes", "No", "No"))
  expect_equal(out$gap_days, c(4, 50, 9))
})

test_that("prep_filter_fungal drops fungal isolates", {
  f <- data.frame(
    organism_group = c("Enterobacterales", "Fungi", "Yeast"),
    organism_normalized = c("Escherichia coli", "Candida albicans", "Candida auris")
  )
  out <- suppressMessages(prep_filter_fungal(f))
  expect_equal(out$organism_normalized, "Escherichia coli")
})

test_that("prep_flag_contaminants flags blood-culture commensals only", {
  d <- data.frame(
    organism_normalized = c("Staphylococcus epidermidis", "Escherichia coli"),
    specimen_type = c("Blood", "Blood")
  )
  out <- suppressMessages(prep_flag_contaminants(d))
  expect_equal(out$is_contaminant, c(TRUE, FALSE))
  expect_true(all(c("contaminant_confidence", "contaminant_method") %in% names(out)))
})

test_that("prep_classify_mdr_xdr follows the Magiorakos thresholds", {
  cls <- c(
    "Aminoglycosides", "Carbapenems", "Third-generation-cephalosporins",
    "Fluoroquinolones", "Penicillins", "Tetracyclines", "Polymyxins"
  )
  d <- data.frame(
    event_id = rep(c("E1", "E2", "E3"), each = length(cls)),
    organism_normalized = "Escherichia coli", org_group = "Enterobacterales",
    antibiotic_class = rep(cls, 3),
    class_result_event = c(
      rep("R", 2), rep("S", 5), # 2 categories: not MDR
      rep("R", 3), rep("S", 4), # 3 categories: MDR
      rep("R", 6), "S"
    ) # all but 1: XDR
  )
  out <- unique(suppressMessages(prep_classify_mdr_xdr(d))[, c("event_id", "mdr", "xdr")])
  expect_equal(out$mdr, c(FALSE, TRUE, TRUE))
  expect_equal(out$xdr, c(FALSE, FALSE, TRUE))
})

test_that("prep_collapse_class_level and prep_create_resistance_profile summarise drug rows", {
  ast <- data.frame(
    event_id = rep(c("E1", "E2"), each = 3),
    organism_normalized = rep(c("Escherichia coli", "Klebsiella pneumoniae"), each = 3),
    antibiotic_normalized = rep(c("Meropenem", "Imipenem", "Ciprofloxacin"), 2),
    antibiotic_class = rep(c("Carbapenems", "Carbapenems", "Fluoroquinolones"), 2),
    antibiotic_value = c("S", "R", "R", "S", "S", "R")
  )
  collapsed <- suppressMessages(prep_collapse_class_level(ast))
  carb <- collapsed[collapsed$antibiotic_class == "Carbapenems", ]
  # any R within a class makes the class R
  expect_equal(carb$class_resistance[carb$event_id == "E1"], "R")
  expect_equal(carb$class_resistance[carb$event_id == "E2"], "S")
  expect_equal(carb$n_drugs_in_class, c(2L, 2L))

  prof <- unique(suppressMessages(prep_create_resistance_profile(ast))[, c("event_id", "resistance_profile", "n_resistant")])
  expect_equal(prof$resistance_profile, c("Ciprofloxacin; Imipenem", "Ciprofloxacin"))
  expect_equal(prof$n_resistant, c(2L, 1L))
  expect_error(prep_create_resistance_profile(ast, format = "nope"), "format must be one of")
})

test_that("plot_infection_type_by_location builds both styles", {
  set.seed(3)
  d <- data.frame(
    PatientInformation_id = sprintf("P%02d", 1:40),
    center_name = rep(c("H1", "H2"), 20),
    location = rep(c("ICU", "Ward"), each = 20),
    infection_type = sample(c("HAI", "CAI"), 40, replace = TRUE)
  )
  for (st in c("bar", "heatmap")) {
    p <- suppressMessages(plot_infection_type_by_location(d, style = st))
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("plot_abx_complex_heatmap returns a ComplexHeatmap object", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(3)
  d <- expand.grid(
    PatientInformation_id = sprintf("P%02d", 1:20),
    antibiotic_name = c("Amikacin", "Meropenem", "Ciprofloxacin"),
    stringsAsFactors = FALSE
  )
  d$center_name <- rep(c("H1", "H2", "H3", "H4"), length.out = nrow(d))
  d$organism_name <- "E. coli"
  d$antibiotic_value <- sample(c("R", "S"), nrow(d), replace = TRUE)
  expect_s4_class(suppressMessages(plot_abx_complex_heatmap(d)), "Heatmap")
})

test_that("print methods for config and result objects", {
  expect_output(print(amr_config()), "hai_cutoff|HAI|config", ignore.case = TRUE)
  r <- run_quiet(raw_amr(), verbose = FALSE, generate_report = FALSE)
  expect_output(print(r))
  cfg <- amr_config()
  capture.output(res <- print(cfg))
  expect_identical(res, cfg)
})

test_that("get_antimicrobial_categories returns group lists and falls back with a warning", {
  expect_true("Carbapenems" %in% get_antimicrobial_categories("Enterobacterales"))
  expect_true("Glycopeptides" %in% get_antimicrobial_categories("Staphylococcus aureus"))
  expect_warning(fb <- get_antimicrobial_categories("Martian bacteria"), "No category list")
  expect_identical(fb, get_antimicrobial_categories("Enterobacterales"))
})

test_that("load_icd_reference reads a user file and errors clearly without one", {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(description3 = "sepsis", icd_code_who_eq = "A41.9"), f, row.names = FALSE)
  expect_equal(load_icd_reference(f)$icd_code_who_eq, "A41.9")
  expect_error(load_icd_reference(tempfile()), "File not found")
})

test_that("legacy aliases point at the prep_* polymicrobial functions", {
  expect_identical(flag_polymicrobial, prep_flag_polymicrobial)
  expect_identical(compute_polymicrobial_weight, prep_compute_poly_weights)
})
