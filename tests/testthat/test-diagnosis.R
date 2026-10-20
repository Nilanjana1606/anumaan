# Diagnosis text -> ICD -> infectious syndrome mapping

icd_ref <- function() {
  data.frame(
    description3 = c("sepsis", "pneumonia", "urinary tract infection"),
    icd_code_who_eq = c("A41.9", "J18.9", "N39.0")
  )
}

test_that("prep_diagnosis_text trims text and falls back when the primary field is empty", {
  d <- data.frame(dx = c("  Sepsis ", NA), fb = c(NA, "Pneumonia"))
  out <- suppressMessages(prep_diagnosis_text(d, "dx", fallback_col = "fb"))
  expect_equal(out$diagnosis_text, c("Sepsis", "Pneumonia"))
})

test_that("prep_map_diagnosis_to_icd exact matching returns one top hit per text", {
  t <- data.frame(diagnosis_text = c("sepsis", "pneumonia", "no such thing"))
  out <- suppressMessages(prep_map_diagnosis_to_icd(t, reference = icd_ref(), method = "exact"))
  expect_equal(out$icd_code, c("A41.9", "J18.9"))
  expect_true(all(out$icd_score == 1))
})

test_that("prep_map_diagnosis_to_icd fuzzy matching ranks the closest description first", {
  skip_if_not_installed("stringdist")
  t <- data.frame(diagnosis_text = "urinary tract infectn")
  out <- suppressMessages(prep_map_diagnosis_to_icd(t, reference = icd_ref(), method = "fuzzy"))
  top <- out[out$icd_rank == 1, ]
  expect_equal(top$icd_code, "N39.0")
  expect_gt(top$icd_score, 0.9)
  expect_false(is.unsorted(rev(out$icd_score)))
})

test_that("prep_map_icd_to_syndrome maps codes and labels unmatched ones", {
  m <- data.frame(icd_code = c("A41.9", "Z99"))
  ref <- data.frame(icd_code = "A41.9", infectious_syndrome = "Bloodstream infections")
  out <- suppressMessages(prep_map_icd_to_syndrome(m, icd_to_syndrome_ref = ref, unmatched_label = "Unmapped"))
  expect_equal(out$syndrome, c("Bloodstream infections", "Unmapped"))
})

test_that("prep_assign_patient_syndrome selects the highest-priority syndrome per patient", {
  a <- data.frame(
    patient_id = c("P1", "P1", "P2"),
    syndrome = c("Bloodstream infections", "Meningitis", "Diarrhoea"),
    icd_score = c(0.9, 0.8, 0.7)
  )
  out <- suppressMessages(prep_assign_patient_syndrome(a))
  sel <- out[out$syndrome_selected, ]
  # Meningitis ranks above Bloodstream infections in the GBD hierarchy
  expect_equal(sel$syndrome[sel$patient_id == "P1"], "Meningitis")
  expect_equal(nrow(sel), 2L)
})

test_that("infer_patient_syndrome_wide picks one syndrome per patient from indicator columns", {
  w <- data.frame(
    patient_id = c("P1", "P2", "P3"),
    Meningitis = c(1, 0, 0),
    `Bloodstream infections` = c(1, 0, 0),
    `Upper respiratory infections` = c(0, 1, 0),
    check.names = FALSE
  )
  out <- suppressMessages(infer_patient_syndrome_wide(w,
    syndrome_cols = c("Meningitis", "Bloodstream infections", "Upper respiratory infections")
  ))
  expect_equal(out$patient_id, c("P1", "P2")) # P3 has no positive syndrome
  expect_equal(out$syndrome, c("Meningitis", "Upper respiratory infections"))
  expect_equal(unname(out$contributes_to_amr_burden), c(TRUE, FALSE))

  burden_only <- suppressMessages(infer_patient_syndrome_wide(w,
    syndrome_cols = c("Meningitis", "Bloodstream infections", "Upper respiratory infections"),
    keep_only_burden_syndromes = TRUE
  ))
  expect_equal(burden_only$patient_id, "P1")
})
