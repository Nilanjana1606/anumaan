# plot_display_labels.R
# Generic human-readable display-label helpers shared by every plotting
# function in this package (PPC report, validation report, diagnostic
# plots). Never alter the underlying data values these labels are derived
# from -- display-only transforms.

.class_display_labels_cache <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      path <- system.file("extdata", "class_display_labels.csv", package = "anumaan")
      cache <<- if (nzchar(path) && file.exists(path)) {
        utils::read.csv(path, stringsAsFactors = FALSE)
      } else {
        data.frame(class_name = character(0), display_label = character(0),
                   short_label = character(0), stringsAsFactors = FALSE)
      }
    }
    cache
  }
})

#' Human-readable antimicrobial class label
#'
#' Looks up \code{inst/extdata/class_display_labels.csv}; falls back to
#' replacing underscores with spaces for any class not in the table (e.g. a
#' pathogen-specific class not yet added there), so this never errors on an
#' unmapped class -- it just degrades to a readable-if-generic label.
#'
#' @param x Character vector of raw class names (e.g. \code{"Third_generation_cephalosporins"}).
#' @return Character vector of display labels, same length as \code{x}.
#' @keywords internal
class_display_label <- function(x) {
  tbl <- .class_display_labels_cache()
  idx <- match(x, tbl$class_name)
  out <- ifelse(is.na(idx), gsub("_", " ", x), tbl$display_label[idx])
  out
}

#' Short antimicrobial class label (for dense facet strips)
#'
#' Same lookup/fallback behaviour as \code{class_display_label()}, but
#' returns the abbreviated form (e.g. \code{"3GC"}). Falls back to the first
#' 4 characters of the underscore-cleaned name, upper-cased, for any class
#' not in the reference table.
#'
#' @param x Character vector of raw class names.
#' @return Character vector of short labels, same length as \code{x}.
#' @keywords internal
class_short_label <- function(x) {
  tbl <- .class_display_labels_cache()
  idx <- match(x, tbl$class_name)
  fallback <- toupper(substr(gsub("_", " ", x), 1L, 4L))
  ifelse(is.na(idx), fallback, tbl$short_label[idx])
}

#' Human-readable class-pair label
#'
#' @param class_1,class_2 Character vectors of raw class names.
#' @param short Logical; use \code{class_short_label()} instead of
#'   \code{class_display_label()}. Default \code{TRUE} (pairwise facet strips
#'   get crowded fast with full names).
#' @return Character vector, e.g. \code{"AMG x CARB"}.
#' @keywords internal
class_pair_label <- function(class_1, class_2, short = TRUE) {
  f <- if (isTRUE(short)) class_short_label else class_display_label
  paste(f(class_1), f(class_2), sep = " x ")
}

#' Human-readable hospital/site display label
#'
#' Replaces underscores with spaces only -- does not alter the underlying
#' `center_name`/hospital identifier used anywhere else (joins, grouping,
#' file paths all keep the raw value; this is a plotting-display transform
#' only).
#'
#' @param x Character vector of raw hospital/site names (e.g. \code{"AIIMS_trauma_center"}).
#' @return Character vector, e.g. \code{"AIIMS trauma center"}.
#' @keywords internal
hospital_display_label <- function(x) {
  gsub("_", " ", x)
}
