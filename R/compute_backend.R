#' Validate a selectable Stan compute backend configuration
#'
#' Validates a compute configuration for CPU or OpenCL execution. This checks
#' only package/toolchain availability and the structure of the requested
#' configuration; it does not inspect machine-specific GPU inventories.
#'
#' @param compute Named list. Supported fields are \code{backend},
#'   \code{opencl_platform_id}, \code{opencl_device_id}, and
#'   \code{allow_cpu_fallback}.
#'
#' @return A canonical validated list with fields
#'   \code{backend}, \code{opencl_platform_id}, \code{opencl_device_id}, and
#'   \code{allow_cpu_fallback}.
#' @export
validate_compute_backend <- function(compute = list()) {
  cfg <- .resolve_compute_config(compute)

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(paste0(
      "Package 'cmdstanr' is required for Bayesian multivariate-probit fitting.\n",
      "  install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))"
    ), call. = FALSE)
  }

  cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  if (!nzchar(cmdstan_path)) {
    stop("CmdStan is not configured. Run cmdstanr::install_cmdstan() or cmdstanr::set_cmdstan_path().",
         call. = FALSE)
  }

  invisible(cfg)
}

.resolve_compute_config <- function(compute = list()) {
  defaults <- list(
    backend = "cpu",
    opencl_platform_id = NULL,
    opencl_device_id = NULL,
    allow_cpu_fallback = FALSE
  )
  if (is.null(compute)) {
    compute <- list()
  }
  if (!is.list(compute)) {
    stop("`compute` must be a named list.", call. = FALSE)
  }

  cfg <- defaults
  for (nm in names(compute)) cfg[[nm]] <- compute[[nm]]

  if (!is.character(cfg$backend) || length(cfg$backend) != 1L || is.na(cfg$backend)) {
    stop("`compute$backend` must be one of 'cpu' or 'opencl'.", call. = FALSE)
  }
  cfg$backend <- tolower(cfg$backend)
  if (!cfg$backend %in% c("cpu", "opencl")) {
    stop(sprintf("Unknown compute backend '%s'. Use 'cpu' or 'opencl'.", cfg$backend),
         call. = FALSE)
  }

  if (!is.logical(cfg$allow_cpu_fallback) || length(cfg$allow_cpu_fallback) != 1L || is.na(cfg$allow_cpu_fallback)) {
    stop("`compute$allow_cpu_fallback` must be TRUE or FALSE.", call. = FALSE)
  }

  .validate_opencl_id <- function(value, field) {
    if (is.null(value)) return(NULL)
    if (length(value) != 1L || is.na(value) || !is.numeric(value) || value < 0 || value != as.integer(value)) {
      stop(sprintf("`compute$%s` must be a non-negative integer scalar.", field),
           call. = FALSE)
    }
    as.integer(value)
  }

  cfg$opencl_platform_id <- .validate_opencl_id(cfg$opencl_platform_id, "opencl_platform_id")
  cfg$opencl_device_id <- .validate_opencl_id(cfg$opencl_device_id, "opencl_device_id")

  if (identical(cfg$backend, "opencl")) {
    if (is.null(cfg$opencl_platform_id) || is.null(cfg$opencl_device_id)) {
      stop("OpenCL backend requires both `compute$opencl_platform_id` and `compute$opencl_device_id`.",
           call. = FALSE)
    }
  } else {
    cfg$opencl_platform_id <- NULL
    cfg$opencl_device_id <- NULL
  }

  cfg
}

.amr_compute_cache_key <- function(stan_code, residual_structure, compute_config, cpp_options = list()) {
  payload <- paste(
    stan_code,
    residual_structure,
    compute_config$backend,
    paste(names(cpp_options), unlist(cpp_options), sep = "=", collapse = "|"),
    sep = "\n--\n"
  )
  tmp <- tempfile("amr_compute_key_", fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(payload, tmp, useBytes = TRUE)
  unname(tools::md5sum(tmp)[[1L]])
}

.amr_compile_stan_backend <- function(stan_code,
                                      residual_structure,
                                      compute_config,
                                      cpp_options = list(),
                                      compile = TRUE,
                                      quiet = FALSE) {
  compute_config <- .resolve_compute_config(compute_config)
  merged_cpp_options <- cpp_options
  if (identical(compute_config$backend, "opencl")) {
    merged_cpp_options$stan_opencl <- TRUE
  }
  cache_key <- .amr_compute_cache_key(
    stan_code = stan_code,
    residual_structure = residual_structure,
    compute_config = compute_config,
    cpp_options = merged_cpp_options
  )
  basename <- sprintf("mvprobit_%s_%s_%s", residual_structure, compute_config$backend, cache_key)
  stan_file <- cmdstanr::write_stan_file(stan_code, basename = basename)
  model <- cmdstanr::cmdstan_model(
    stan_file,
    compile = compile,
    cpp_options = merged_cpp_options,
    quiet = quiet
  )
  list(
    model = model,
    stan_file = stan_file,
    basename = basename,
    cache_key = cache_key,
    cpp_options = merged_cpp_options,
    compute_config = compute_config
  )
}

.amr_sample_with_backend <- function(model, sample_args, compute_config) {
  compute_config <- .resolve_compute_config(compute_config)
  if (identical(compute_config$backend, "opencl")) {
    sample_args$opencl_ids <- c(
      compute_config$opencl_platform_id,
      compute_config$opencl_device_id
    )
  }
  do.call(model$sample, sample_args)
}
