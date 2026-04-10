# zzz_utils_internal.R
# Internal utility helpers (not exported)

#' Parse Age Bin Labels to Breaks
#' @keywords internal
parse_age_bin_labels <- function(labels) {
  breaks <- numeric(length(labels) + 1)
  clean_labels <- character(length(labels))

  for (i in seq_along(labels)) {
    label <- labels[i]

    if (grepl("\\+$", label)) {
      # Handle "85+" format
      lower <- as.numeric(gsub("\\+", "", label))
      breaks[i] <- lower
      breaks[i + 1] <- Inf
      clean_labels[i] <- label
    } else if (grepl("^<", label)) {
      # Handle "<1" format -- lower = -Inf to capture ages like -1
      upper <- as.numeric(gsub("^<", "", label))
      breaks[i] <- -Inf
      breaks[i + 1] <- upper
      clean_labels[i] <- label
    } else if (grepl("-", label)) {
      # Handle "1-5" format
      parts <- strsplit(label, "-")[[1]]
      lower <- as.numeric(parts[1])
      upper <- as.numeric(parts[2])
      breaks[i] <- lower
      breaks[i + 1] <- upper
      clean_labels[i] <- label
    } else {
      stop(sprintf("Cannot parse age bin label: %s", label))
    }
  }

  # Remove duplicate breaks
  breaks <- unique(breaks)

  return(list(breaks = breaks, labels = clean_labels))
}
# utils.R
# Small utility functions used across modules

#' Largest-Remainder Rounding
#'
#' Rounds a numeric vector so that the individual rounded values sum exactly to
#' a target total.  Uses the largest-remainder method (Hamilton method).
#'
#' @param x Numeric vector to round.
#' @param target Integer target sum. Default is \code{round(sum(x))}.
#'
#' @return Integer vector of the same length as \code{x} whose sum equals
#'   \code{target}.
#' @export
#'
#' @examples
#' round_to_sum(c(3.3, 3.3, 3.4), target = 10)
round_to_sum <- function(x, target = round(sum(x))) {
  floors <- floor(x)
  remainder <- target - sum(floors)
  ranks <- order(x - floors, decreasing = TRUE)
  floors[ranks[seq_len(remainder)]] <- floors[ranks[seq_len(remainder)]] + 1L
  floors
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
