test_that("Gibbs stability is not applicable to analytic identity completion", {
  fit <- list(residual_structure = "identity")
  out <- assess_gibbs_profile_stability(fit)
  expect_identical(out$status, "not_applicable")
  expect_identical(out$summary$gibbs_stability_status, "not_applicable")
  expect_length(out$posterior_draw_indices, 0L)
})

test_that("Gibbs stability validates schedules before attempting any model work", {
  expect_error(
    assess_gibbs_profile_stability(list(residual_structure = "identity"), schedules = list(baseline = list(burnin = 1))),
    "baseline, medium, and long"
  )
})
