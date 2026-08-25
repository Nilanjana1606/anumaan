# Regression test for the stale parameter-name matching bug in
# R/plot_probit_diagnostics.R: the trace/rank-plot parameter selection used
# to search for pre-generic-architecture names (hospital_effect[, R_hospital[,
# etc.) that never appear in the current package's Stan output, so re_vars/
# corr_vars were silently always empty. This does not exercise the full PDF
# rendering path (which needs bayesplot + a real fit); it isolates just the
# parameter-selection regexes, mirroring the code in plot_probit_diagnostics().

test_that("diagnostic parameter selection finds current generic re_effect[]/R_block[] names", {
  all_vars <- c(
    "lp__",
    "tau_re[1]", "tau_re[2]",
    "beta[1,1]", "beta[1,2]", "beta[2,1]",
    "re_effect[1,1]", "re_effect[1,2]", "re_effect[2,1]",
    "R_block[1,1,1]", "R_block[1,1,2]", "R_block[1,2,2]"
  )

  re_vars <- grep("^re_effect\\[", all_vars, value = TRUE)
  corr_vars <- grep("^R_block\\[", all_vars, value = TRUE)

  expect_identical(re_vars, c("re_effect[1,1]", "re_effect[1,2]", "re_effect[2,1]"))
  expect_identical(corr_vars, c("R_block[1,1,1]", "R_block[1,1,2]", "R_block[1,2,2]"))
})

test_that("diagnostic parameter selection ignores stale pre-generic parameter names", {
  all_vars <- c(
    "lp__",
    "hospital_effect[1]", "patient_effect[1]", "admission_effect[1]",
    "R_hospital[1,1,1]", "R_patient[1,1,1]", "R_admission[1,1,1]",
    "re_effect[1,1]", "R_block[1,1,1]"
  )

  re_vars <- grep("^re_effect\\[", all_vars, value = TRUE)
  corr_vars <- grep("^R_block\\[", all_vars, value = TRUE)

  expect_identical(re_vars, "re_effect[1,1]")
  expect_identical(corr_vars, "R_block[1,1,1]")
  expect_false(any(grepl("hospital_effect|patient_effect|admission_effect", re_vars)))
  expect_false(any(grepl("R_hospital|R_patient|R_admission", corr_vars)))
})
