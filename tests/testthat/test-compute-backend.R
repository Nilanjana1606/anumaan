skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  testthat::skip_if(!nzchar(path), "CmdStan is not configured")
}

testthat::test_that("compute backend defaults to CPU", {
  skip_if_no_cmdstan()
  cfg <- anumaan::validate_compute_backend(list())
  testthat::expect_identical(cfg$backend, "cpu")
  testthat::expect_false(cfg$allow_cpu_fallback)
  testthat::expect_null(cfg$opencl_platform_id)
  testthat::expect_null(cfg$opencl_device_id)
})

testthat::test_that("opencl backend requires ids", {
  skip_if_no_cmdstan()
  testthat::expect_error(
    anumaan::validate_compute_backend(list(backend = "opencl")),
    "requires both"
  )
})

testthat::test_that("cpu backend ignores supplied opencl ids", {
  cfg <- anumaan:::.resolve_compute_config(list(
    backend = "cpu",
    opencl_platform_id = 0L,
    opencl_device_id = 0L
  ))
  testthat::expect_identical(cfg$backend, "cpu")
  testthat::expect_null(cfg$opencl_platform_id)
  testthat::expect_null(cfg$opencl_device_id)
})

testthat::test_that("cache keys differ by backend", {
  skip_if_no_cmdstan()
  stan_code <- "parameters {real y;} model {y ~ normal(0,1);}"
  cpu <- anumaan::validate_compute_backend(list(backend = "cpu"))
  gpu <- anumaan:::.resolve_compute_config(list(
    backend = "opencl",
    opencl_platform_id = 0L,
    opencl_device_id = 0L
  ))
  key_cpu <- anumaan:::.amr_compute_cache_key(stan_code, "identity", cpu, list())
  key_gpu <- anumaan:::.amr_compute_cache_key(stan_code, "identity", gpu, list(stan_opencl = TRUE))
  testthat::expect_false(identical(key_cpu, key_gpu))
})
