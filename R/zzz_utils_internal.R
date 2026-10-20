# zzz_utils_internal.R
# Internal utility helpers (not exported)

#' Parse Age Bin Labels to Breaks
#'
#' Parses age-bin labels into numeric breakpoints for interval binning.
#'
#' @param labels Character vector of age-bin labels, such as \code{"<1"},
#'   \code{"1-5"}, or \code{"85+"}.
#'
#' @return A list with \code{breaks} (numeric breakpoints) and \code{labels}
#'   (original labels).
#'
#' @keywords internal
parse_age_bin_labels <- function(labels) {
  clean_labels <- as.character(labels)
  breaks <- numeric(0)

  for (i in seq_along(clean_labels)) {
    label <- clean_labels[i]

    if (is.na(label) || !nzchar(trimws(label))) {
      stop(sprintf("Cannot parse age bin label: %s", label))
    }

    if (grepl("^<\\s*[0-9.]+\\s*$", label)) {
      # Handle "<1" format.
      upper <- as.numeric(gsub("^<\\s*", "", label))
      breaks <- c(breaks, -Inf, upper)
    } else if (grepl("^\\s*[0-9.]+\\s*-\\s*[0-9.]+\\s*$", label)) {
      # Handle "1-5" format.
      parts <- strsplit(gsub("\\s+", "", label), "-", fixed = TRUE)[[1]]
      lower <- as.numeric(parts[1])
      upper <- as.numeric(parts[2])
      breaks <- c(breaks, lower, upper)
    } else if (grepl("^\\s*[0-9.]+\\s*\\+\\s*$", label)) {
      # Handle "85+" format.
      lower <- as.numeric(gsub("\\+\\s*$", "", gsub("\\s+", "", label)))
      breaks <- c(breaks, lower, Inf)
    } else {
      stop(sprintf("Cannot parse age bin label: %s", label))
    }
  }

  list(breaks = unique(breaks), labels = clean_labels)
}

#' Largest-Remainder Rounding
#'
#' Rounds a numeric vector so that the individual rounded values sum exactly to
#' a target total. Uses the largest-remainder method (Hamilton method).
#' Input is coerced with \code{as.numeric()}; non-numeric values become
#' \code{NA}.
#'
#' @param x Numeric vector to round.
#' @param target Integer target sum. Default is
#'   \code{round(sum(x, na.rm = TRUE))}.
#'
#' @return Integer vector of the same length as \code{x} whose sum equals
#'   \code{target} across non-missing entries; missing values are preserved as
#'   \code{NA}.
#' @export
#'
#' @examples
#' round_to_sum(c(3.3, 3.3, 3.4), target = 10)
round_to_sum <- function(x, target = round(sum(x, na.rm = TRUE))) {
  x_num <- suppressWarnings(as.numeric(x))
  out <- rep(NA_integer_, length(x_num))
  idx <- which(!is.na(x_num))

  if (!length(idx)) {
    return(out)
  }

  floors <- floor(x_num[idx])
  fracs <- x_num[idx] - floors
  n <- length(floors)
  remainder <- as.integer(round(target - sum(floors)))

  if (remainder > 0L) {
    full <- remainder %/% n
    extra <- remainder %% n
    floors <- floors + full
    if (extra > 0L) {
      ranks <- order(fracs, decreasing = TRUE)
      floors[ranks[seq_len(extra)]] <- floors[ranks[seq_len(extra)]] + 1L
    }
  } else if (remainder < 0L) {
    take <- -remainder
    full <- take %/% n
    extra <- take %% n
    floors <- floors - full
    if (extra > 0L) {
      ranks <- order(fracs, decreasing = FALSE)
      floors[ranks[seq_len(extra)]] <- floors[ranks[seq_len(extra)]] - 1L
    }
  }

  out[idx] <- as.integer(floors)
  out
}


#' Shorten Antibiotic Class Names
#'
#' Maps long antibiotic class names to common abbreviations used in GBD-style
#' figures.
#'
#' @param x Character vector of antibiotic class names.
#'
#' @return Character vector of shortened names.
#' @export
#'
#' @examples
#' shorten_drug_class(c("Carbapenems", "Third-generation-cephalosporins"))
shorten_drug_class <- function(x) {
  x <- sub("_R$", "", x)
  dplyr::case_when(
    x == "Third-generation-cephalosporins" ~ "3GC",
    x == "Fourth-generation-cephalosporins" ~ "4GC",
    x == "Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal" ~ "BL/BLI-AP",
    x == "Beta-lactam/beta-lactamase-inhibitor" ~ "BL/BLI",
    x == "Trimethoprim-sulfamethoxazole" ~ "TMP-SMX",
    x == "Anti-pseudomonal-penicillins" ~ "Anti-pseudo-PCN",
    TRUE ~ x
  )
}
normalize_join <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9 ]", "", x)
  x <- gsub("\\s+", " ", x)
  x
}

# ponytail: print() to console can't be suppressed; route tables through message()
.msg_print <- function(x, ...) {
  message(paste(utils::capture.output(print(x, ...)), collapse = "\n"))
  invisible(x)
}

# Seed the RNG for the caller and restore the user's .Random.seed when the
# caller exits (CRAN: functions must not leave global RNG state changed).
.local_seed <- function(seed, envir = parent.frame()) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv())
  }
  do.call(on.exit, list(bquote(.restore_seed(.(old))), add = TRUE, after = FALSE),
    envir = envir
  )
  set.seed(seed)
}

.restore_seed <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", old, envir = globalenv())
  }
}

# osqp >= 1.0 renamed `polish` to `polishing`; pass whichever this version takes.
.osqp_settings <- function(...) {
  polish_arg <- intersect(c("polishing", "polish"), names(formals(osqp::osqpSettings)))[1]
  args <- list(...)
  if (!is.na(polish_arg)) args[[polish_arg]] <- TRUE
  do.call(osqp::osqpSettings, args)
}

# Canonical internal key for organism and drug names: lower case, trimmed,
# single-spaced. Every comparison against reference tables or user maps goes
# through this so "Klebsiella pneumoniae", "klebsiella  pneumoniae " match.
.name_key <- function(x) {
  gsub("\\s+", " ", tolower(trimws(as.character(x))))
}

# Re-spell `x` (map names or a column) to match `reference` case-insensitively,
# so matching ignores case but outputs keep the data's own labels.
.align_names <- function(x, reference) {
  hit <- match(.name_key(x), .name_key(reference))
  ifelse(is.na(hit), x, as.character(reference)[hit])
}

.align_panel_map <- function(panel_map, values) {
  if (!is.null(panel_map)) names(panel_map) <- .align_names(names(panel_map), unique(values))
  panel_map
}

.align_col <- function(df, col, values) {
  if (!is.null(df) && !is.null(col) && col %in% names(df)) {
    df[[col]] <- .align_names(df[[col]], unique(values))
  }
  df
}
