# random_effects_generic.R
# Stage 1 of the generic random-effect architecture: an arbitrary number of
# random-INTERCEPT blocks (random slopes are explicitly out of scope for this
# commit -- see .normalize_random_effects_spec()). Random effects are
# user-declared blocks (name + group_col [+ terms]), never inferred from
# column names that merely look like IDs, and never exposed to the user as
# anonymous numbered compartments (re_1/re_2/re_3).

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Normalize a random_effects specification into the canonical block-list
#' form used throughout the package.
#'
#' Accepts EITHER:
#'  - the legacy character vector, e.g. c("center_name", "readmission_id")
#'    (temporary backward compatibility -- each column becomes its own block,
#'    using the COLUMN NAME itself as the block name, never a numbered label); or
#'  - the new list-of-blocks form, e.g.
#'      list(list(name = "hospital", group_col = "center_name", terms = "intercept"),
#'           list(name = "admission", group_col = "admission_key", terms = "intercept"))
#'
#' @keywords internal
.normalize_random_effects_spec <- function(random_effects) {
  if (is.character(random_effects)) {
    if (anyDuplicated(random_effects)) {
      stop("`random_effects` column names must be unique.", call. = FALSE)
    }
    return(lapply(random_effects, function(col) {
      list(name = col, group_col = col, terms = "intercept")
    }))
  }

  if (!is.list(random_effects)) {
    stop("`random_effects` must be a character vector or a list of blocks (name/group_col/terms).",
      call. = FALSE
    )
  }

  for (i in seq_along(random_effects)) {
    b <- random_effects[[i]]
    if (!is.list(b)) {
      stop(sprintf("random_effects[[%d]] must be a list with name/group_col/terms.", i), call. = FALSE)
    }
    if (is.null(b$name) || !nzchar(b$name)) {
      stop(sprintf("random_effects[[%d]] is missing a `name`. Anonymous numbered blocks (re_1, re_2, ...) are not supported -- every block must have an explicit, meaningful name.", i),
        call. = FALSE
      )
    }
    if (is.null(b$group_col) || !nzchar(b$group_col)) {
      stop(sprintf("random_effects block '%s' is missing a `group_col`.", b$name), call. = FALSE)
    }
    terms <- b$terms %||% "intercept"
    terms_chr <- if (is.list(terms)) unlist(terms) else terms
    if (!identical(terms_chr, "intercept")) {
      stop(sprintf(
        "random_effects block '%s': terms = %s is not supported. Stage 1 of the generic random-effect architecture supports random-intercept-only blocks (terms: [intercept]); random slopes are deferred to a later commit.",
        b$name, paste(terms_chr, collapse = ", ")
      ), call. = FALSE)
    }
  }
  random_effects
}

#' Build the generic random-effect representation for an arbitrary number of blocks
#'
#' Flattens an arbitrary number of declared random-intercept blocks into the
#' representation used by the generic Stan models and by every downstream
#' mu-reconstruction site (profile generation, DALY aggregation, validation).
#'
#' @param data Data frame containing every declared block's group_col.
#' @param random_effects Character vector (legacy) or list-of-blocks (see
#'   .normalize_random_effects_spec()).
#' @param min_repeated_levels Optional integer. If supplied, a block whose
#'   number of levels with >=2 observations falls below this is checked
#'   INDEPENDENTLY of singleton_threshold below (previously this was bugged:
#'   nested inside the singleton_fraction check, so it was silently skipped
#'   whenever singleton_fraction fell at or below the threshold even if
#'   n_repeated_levels was itself too low). Triggers a stop (via
#'   on_mostly_singleton = "stop") or warning. NULL (default) never checks
#'   this -- eligibility gating belongs in the analysis layer (see
#'   anumaan-analysis's random_effect_eligibility.R), not hardcoded into the
#'   package.
#' @param singleton_threshold Numeric in (0, 1]. A block whose fraction of
#'   singleton levels (exactly 1 observation) exceeds this triggers a
#'   separate warning/stop, evaluated independently of min_repeated_levels.
#'   Default 0.90.
#' @param on_mostly_singleton "warn" (default) or "stop" -- applies to BOTH
#'   the min_repeated_levels check and the singleton_threshold check.
#' @return An object of class "amr_random_effects": block metadata, level
#'   maps, per-event flattened group indices, nesting diagnostics.
#' @export
prepare_random_effects <- function(data, random_effects,
                                    min_repeated_levels = NULL,
                                    singleton_threshold = 0.90,
                                    on_mostly_singleton = c("warn", "stop")) {
  on_mostly_singleton <- match.arg(on_mostly_singleton)
  blocks <- .normalize_random_effects_spec(random_effects)

  block_names <- vapply(blocks, function(b) b$name, character(1L))
  if (anyDuplicated(block_names)) {
    stop("random_effects block names must be unique: ",
      paste(block_names[duplicated(block_names)], collapse = ", "),
      call. = FALSE
    )
  }

  group_cols <- vapply(blocks, function(b) b$group_col, character(1L))
  missing_cols <- setdiff(group_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("random_effects group_col(s) not found in data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  R <- length(blocks)
  n_events <- nrow(data)

  level_maps <- vector("list", R)
  n_levels <- integer(R)
  group_index <- matrix(NA_integer_, nrow = n_events, ncol = R) # 1-indexed, local to block
  n_obs_per_level <- vector("list", R)

  for (r in seq_len(R)) {
    v <- data[[group_cols[r]]]
    n_missing <- sum(is.na(v) | (is.character(v) & !nzchar(v)))
    if (n_missing > 0L) {
      stop(sprintf(
        "random_effects block '%s' (group_col '%s') has %d missing/empty group id(s).",
        block_names[r], group_cols[r], n_missing
      ), call. = FALSE)
    }
    lvls <- sort(unique(v))
    level_maps[[r]] <- lvls
    n_levels[r] <- length(lvls)
    idx <- match(v, lvls)
    group_index[, r] <- idx
    n_obs_per_level[[r]] <- as.integer(table(factor(idx, levels = seq_along(lvls))))
  }

  level_start <- cumsum(c(1L, n_levels))[seq_len(R)]
  level_end <- level_start + n_levels - 1L
  total_re_levels <- sum(n_levels)

  flat_group_index <- matrix(NA_integer_, nrow = n_events, ncol = R)
  for (r in seq_len(R)) flat_group_index[, r] <- group_index[, r] + level_start[r] - 1L

  # Full pairwise relationship table between EVERY ordered pair of declared
  # blocks (not just consecutive ones), determined from the actual data --
  # never assumed from declaration order. Declaring blocks as
  # hospital/admission/patient (in that order) does NOT mean patient is
  # meaningfully related only to admission; a patient with several
  # admissions is nested-within-patient from admission's perspective
  # regardless of where "patient" appears in the declaration list.
  #   nested:              every level of child maps to exactly one parent level
  #   identical_partition: nested AND every parent level maps to exactly one
  #                        child level (the two blocks partition events identically)
  #   crossed_or_non_nested: otherwise
  pairwise_relationships <- if (R >= 2L) {
    rows <- list()
    for (child in seq_len(R)) {
      for (parent in seq_len(R)) {
        if (child == parent) next
        parents_per_child_level <- tapply(group_index[, parent], group_index[, child], function(x) length(unique(x)))
        n_child_one_parent  <- sum(parents_per_child_level == 1L)
        n_child_multi_parent <- sum(parents_per_child_level > 1L)
        is_nested <- n_child_multi_parent == 0L
        children_per_parent_level <- tapply(group_index[, child], group_index[, parent], function(x) length(unique(x)))
        is_reverse_nested <- all(children_per_parent_level == 1L)
        relationship <- if (is_nested && is_reverse_nested) "identical_partition"
                        else if (is_nested) "nested"
                        else "crossed_or_non_nested"
        rows[[length(rows) + 1L]] <- data.frame(
          child_block = block_names[child],
          parent_block = block_names[parent],
          n_child_levels = n_levels[child],
          n_child_levels_mapping_to_one_parent = n_child_one_parent,
          n_child_levels_mapping_to_multiple_parents = n_child_multi_parent,
          relationship = relationship,
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  } else {
    data.frame(child_block = character(0), parent_block = character(0),
               n_child_levels = integer(0),
               n_child_levels_mapping_to_one_parent = integer(0),
               n_child_levels_mapping_to_multiple_parents = integer(0),
               relationship = character(0), stringsAsFactors = FALSE)
  }

  # Compact per-block summary retained for print()/backward-compat callers --
  # block r's relationship to block r-1 specifically. This is NOT a
  # substitute for pairwise_relationships above, which covers every pair.
  nesting <- character(R)
  if (R > 0L) nesting[1L] <- "root"
  if (R >= 2L) {
    for (r in 2:R) {
      rel <- pairwise_relationships$relationship[
        pairwise_relationships$child_block == block_names[r] &
        pairwise_relationships$parent_block == block_names[r - 1L]]
      nesting[r] <- if (length(rel) == 1L) rel else NA_character_
    }
  }

  singleton_fraction <- vapply(n_obs_per_level, function(x) mean(x == 1L), numeric(1L))
  .signal <- function(msg) {
    if (identical(on_mostly_singleton, "stop")) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  for (r in seq_len(R)) {
    n_repeated <- sum(n_obs_per_level[[r]] >= 2L)

    # Two INDEPENDENT checks -- previously min_repeated_levels was nested
    # inside the singleton_fraction condition, so a block could fall below
    # min_repeated_levels while its singleton_fraction was still <= the
    # threshold and the check would silently never fire.
    if (!is.null(min_repeated_levels) && n_repeated < min_repeated_levels) {
      .signal(sprintf(
        "random_effects block '%s': only %d of its %d levels have >=2 observations (need >= %d) -- this block's variance may be weakly identified.",
        block_names[r], n_repeated, n_levels[r], min_repeated_levels))
    }

    if (singleton_fraction[r] > singleton_threshold) {
      .signal(sprintf(
        "random_effects block '%s': %.1f%% of its %d levels are singletons (exactly 1 observation), above the %.0f%% threshold -- this block's variance may be weakly identified.",
        block_names[r], 100 * singleton_fraction[r], n_levels[r], 100 * singleton_threshold))
    }
  }

  structure(list(
    blocks             = blocks,
    block_names        = block_names,
    group_cols         = group_cols,
    R                  = R,
    n_levels           = n_levels,
    level_start        = level_start,
    level_end          = level_end,
    total_re_levels    = total_re_levels,
    level_maps         = level_maps,
    group_index        = group_index,
    flat_group_index   = flat_group_index,
    n_obs_per_level    = n_obs_per_level,
    nesting            = nesting,
    pairwise_relationships = pairwise_relationships,
    singleton_fraction = singleton_fraction
  ), class = "amr_random_effects")
}

#' @export
print.amr_random_effects <- function(x, ...) {
  cat(sprintf("<amr_random_effects> %d block(s), %d total levels\n", x$R, x$total_re_levels))
  for (r in seq_len(x$R)) {
    cat(sprintf(
      "  [%d] %-15s group_col=%-15s n_levels=%-6d nesting=%s\n",
      r, x$block_names[r], x$group_cols[r], x$n_levels[r], x$nesting[r]
    ))
  }
  invisible(x)
}

#' Compute each event's total random-effect contribution from a flattened
#' re_effect[D, total_re_levels] matrix (as emitted by the generic Stan
#' models' generated re_effect, or reconstructed R-side from posterior
#' draws of z_re/tau_re/L_corr_re). This is the SINGLE generic helper every
#' downstream mu-reconstruction site should call instead of hand-summing
#' hospital_effect/patient_effect/admission_effect.
#'
#' @param re_effect D x total_re_levels numeric matrix (one posterior draw).
#' @param flat_group_index n_events x R integer matrix (see
#'   prepare_random_effects()$flat_group_index), or a single event's R-length
#'   integer vector.
#' @return If flat_group_index is a matrix: n_events x D matrix of summed RE
#'   contributions. If a vector: length-D vector for one event.
#' @export
re_contribution <- function(re_effect, flat_group_index) {
  if (is.matrix(flat_group_index)) {
    n_events <- nrow(flat_group_index)
    R <- ncol(flat_group_index)
    D <- nrow(re_effect)
    out <- matrix(0, nrow = n_events, ncol = D)
    for (r in seq_len(R)) out <- out + t(re_effect[, flat_group_index[, r], drop = FALSE])
    out
  } else {
    R <- length(flat_group_index)
    if (R == 0L) return(rep(0, nrow(re_effect)))
    Reduce(`+`, lapply(seq_len(R), function(r) re_effect[, flat_group_index[r]]))
  }
}
