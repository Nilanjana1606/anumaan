# Test sp. -> spp. Normalization Fix

test_that("prep_standardize_organisms converts sp/sp. variants to spp.", {
  test_data <- data.frame(
    patient_id = 1:10,
    organism_name = c(
      "Acinetobacter sp.",
      "Acinetobacter sp",
      "acinetobacter sp.",
      "ACINETOBACTER SP",
      "Pseudomonas sp.",
      "Pseudomonas sp",
      "Enterococcus sp.",
      "Proteus sp",
      "Acinetobacter baumannii",
      "Enterobacter sp."
    ),
    stringsAsFactors = FALSE
  )

  result <- prep_standardize_organisms(test_data, organism_col = "organism_name")

  # All sp/sp. variants (excluding full species names) should become spp.
  sp_rows <- grepl("\\bsp\\.?$", test_data$organism_name, ignore.case = TRUE)
  expect_true(all(grepl("spp\\.", result$organism_normalized[sp_rows])))

  # Full species name should be preserved
  baumannii_row <- test_data$organism_name == "Acinetobacter baumannii"
  expect_true(grepl("baumannii", result$organism_normalized[baumannii_row]))
})
