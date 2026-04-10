# prep_stewardship_join.R
# Automated joining pipeline for ICMR stewardship multi-centre data.
#
# Covers the full workflow extracted from Data_merging_all_centres.Rmd:
#   1. Pre-join sanity checks  (column presence, types, key quality)
#   2. Date detection & repair (Excel serial, ISO string, encrypted/NULL)
#   3. Per-centre join         (patient -> abxinfo -> audit -> HAI -> organism bundle)
#   4. Lookup enrichment       (district+state, sample_type, organisms, antibiotics, center)
#   5. Duplicate-column reconciliation
#   6. Post-join cleaning      (placeholders -> NA, fungal exclusion, analysis-ready filter)
#   7. Multi-centre binding    (harmonise columns, cross-centre QC)
#
# Centres handled out-of-the-box: SRGH, KGMU, JIPMER, KMC, HINDUJA, PGIMER.
# Any additional centre follows the same pattern via prep_join_stewardship_centre().


# ---------------------------------------------------------------------------
# Section 1: Pre-join sanity checks
# ---------------------------------------------------------------------------

#' Check required columns exist and report types
#'
#' Runs before any table is merged. Reports missing columns, type mismatches,
#' and blank/duplicate column names.  Stops on hard failures; warns on soft ones.
#'
#' @param data Data frame to inspect.
#' @param required Character vector of column names that must be present.
#' @param expected_types Named character vector mapping column names to expected
#'   R classes (e.g. \code{c(PatientInformation_id = "character")}). Only
#'   checked when the column exists.
#' @param table_label Character. Label used in messages (e.g. "patient_KGMU").
#' @param stop_on_missing Logical. Stop if required columns are absent (default TRUE).
#'
#' @return Invisibly returns a tibble summarising every checked column.
#' @export
prep_check_columns <- function(data,
                               required = character(),
                               expected_types = character(),
                               table_label = "table",
                               stop_on_missing = TRUE) {
  if (!is.data.frame(data)) {
    stop(sprintf("[%s] `data` must be a data frame.", table_label))
  }

  # Blank / duplicate column names
  blank_names  <- sum(is.na(names(data)) | trimws(names(data)) == "")
  dup_names    <- sum(duplicated(names(data)))
  if (blank_names > 0L) {
    warning(sprintf("[%s] %d blank column name(s) found.", table_label, blank_names))
  }
  if (dup_names > 0L) {
    warning(sprintf("[%s] %d duplicated column name(s): %s",
                    table_label,
                    dup_names,
                    paste(names(data)[duplicated(names(data))], collapse = ", ")))
  }

  # Missing required columns
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    msg <- sprintf("[%s] Missing required column(s): %s",
                   table_label, paste(missing, collapse = ", "))
    if (stop_on_missing) stop(msg) else warning(msg)
  }

  # Type checks on columns that exist
  report_rows <- list()
  for (col in names(data)) {
    actual_class  <- paste(class(data[[col]]), collapse = "/")
    expected_class <- if (col %in% names(expected_types)) expected_types[[col]] else NA_character_
    type_ok <- if (!is.na(expected_class)) {
      inherits(data[[col]], expected_class)
    } else {
      TRUE
    }
    n_na      <- sum(is.na(data[[col]]))
    n_total   <- nrow(data)
    pct_na    <- if (n_total > 0) round(100 * n_na / n_total, 1) else NA_real_

    report_rows[[length(report_rows) + 1L]] <- data.frame(
      table        = table_label,
      column       = col,
      required     = col %in% required,
      present      = TRUE,
      actual_type  = actual_class,
      expected_type = if (!is.na(expected_class)) expected_class else "",
      type_ok      = type_ok,
      n_total      = n_total,
      n_na         = n_na,
      pct_na       = pct_na,
      stringsAsFactors = FALSE
    )

    if (!type_ok) {
      warning(sprintf(
        "[%s] Column '%s': expected class '%s', found '%s'.",
        table_label, col, expected_class, actual_class
      ))
    }
  }

  # Add rows for required columns that are absent
  for (col in missing) {
    report_rows[[length(report_rows) + 1L]] <- data.frame(
      table        = table_label,
      column       = col,
      required     = TRUE,
      present      = FALSE,
      actual_type  = NA_character_,
      expected_type = if (col %in% names(expected_types)) expected_types[[col]] else "",
      type_ok      = FALSE,
      n_total      = NA_integer_,
      n_na         = NA_integer_,
      pct_na       = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  report <- dplyr::bind_rows(report_rows)

  message(sprintf(
    "[%s] Column check: %d cols | %d required | %d missing | %d type warnings",
    table_label,
    ncol(data),
    length(required),
    length(missing),
    sum(!report$type_ok, na.rm = TRUE)
  ))

  invisible(report)
}


#' Check join key quality
#'
#' Reports missing, blank, duplicate, and placeholder keys in a column before
#' a merge.  Run this on both sides of every planned join.
#'
#' @param data Data frame.
#' @param key_col Character. Name of the key column.
#' @param table_label Character. Label used in messages.
#' @param warn_missing_pct Numeric. Warn when proportion of missing keys exceeds
#'   this threshold (0-100). Default 5.
#'
#' @return Invisibly returns a one-row summary tibble.
#' @export
prep_check_keys <- function(data,
                            key_col,
                            table_label = "table",
                            warn_missing_pct = 5) {
  if (!key_col %in% names(data)) {
    warning(sprintf("[%s] Key column '%s' not found.", table_label, key_col))
    return(invisible(NULL))
  }

  raw <- as.character(data[[key_col]])
  placeholders <- c("", "NULL", "null", "NA", "N/A", "None", "none", "nan", "NaN")
  is_missing <- is.na(raw) | trimws(raw) %in% placeholders
  n_total    <- length(raw)
  n_missing  <- sum(is_missing)
  pct_missing <- round(100 * n_missing / max(n_total, 1), 1)
  n_distinct <- length(unique(raw[!is_missing]))
  n_dup      <- sum(duplicated(raw[!is_missing]))

  if (pct_missing > warn_missing_pct) {
    warning(sprintf(
      "[%s] Key '%s': %.1f%% missing/placeholder (%d of %d rows).",
      table_label, key_col, pct_missing, n_missing, n_total
    ))
  }

  message(sprintf(
    "[%s] Key '%s': %d rows | %d missing (%.1f%%) | %d distinct | %d duplicated",
    table_label, key_col, n_total, n_missing, pct_missing, n_distinct, n_dup
  ))

  invisible(data.frame(
    table          = table_label,
    key_col        = key_col,
    n_total        = n_total,
    n_missing      = n_missing,
    pct_missing    = pct_missing,
    n_distinct     = n_distinct,
    n_duplicated   = n_dup,
    stringsAsFactors = FALSE
  ))
}


#' Detect and decode encrypted or non-standard date columns
#'
#' Inspects a column for date-like content and automatically selects the
#' correct parsing strategy:
#'
#' \describe{
#'   \item{Excel serial}{Numeric values in the range 10 000-99 999 (days since
#'     1899-12-30). Also detects when the column is character but every non-NA
#'     value is numeric-looking in that range.}
#'   \item{Unix timestamp ms}{Numeric > 1e10 -- milliseconds since Unix epoch.}
#'   \item{ISO / DMY / MDY strings}{Parsed via \code{lubridate::parse_date_time}
#'     with 12 format orders.}
#'   \item{Reversed YYYYDDMM / DDMMYYYY}{Detected when standard parsers fail
#'     but swapping day/month produces a plausible date.}
#'   \item{Encrypted / undecodable}{Values that are long (>12 chars), contain
#'     no digit-separator pattern, and cannot be parsed are flagged. The column
#'     is still returned (as NA) with a warning so the pipeline does not crash.}
#' }
#'
#' @param x Vector (character, numeric, or Date-like) to parse.
#' @param col_name Character. Column name used in messages.
#' @param table_label Character. Table label used in messages.
#'
#' @return A Date vector of the same length as \code{x}.
#' @export
prep_parse_date_column <- function(x, col_name = "date", table_label = "table") {
  n <- length(x)

  # Already a Date
  if (inherits(x, "Date")) {
    return(x)
  }
  # POSIXct / POSIXt
  if (inherits(x, c("POSIXct", "POSIXt"))) {
    return(as.Date(x))
  }

  out <- rep(as.Date(NA_character_), n)

  # Coerce to character for inspection
  x_chr <- trimws(as.character(x))
  placeholders <- c("", "NULL", "null", "NA", "N/A", "None", "none", "nan", "NaN")
  x_chr[x_chr %in% placeholders] <- NA_character_

  idx_present <- which(!is.na(x_chr))
  if (length(idx_present) == 0L) {
    message(sprintf("[%s] '%s': all values missing - returning NA Date vector.",
                    table_label, col_name))
    return(out)
  }

  # ---- Numeric probe -------------------------------------------------------
  x_num <- suppressWarnings(as.numeric(x_chr))
  is_numeric <- !is.na(x_num)

  # Excel serial dates: 10 000-99 999  (roughly 1927-2173)
  excel_idx <- idx_present[is_numeric[idx_present] & x_num[idx_present] > 10000 &
                             x_num[idx_present] < 100000]
  if (length(excel_idx) > 0L) {
    out[excel_idx] <- as.Date(x_num[excel_idx], origin = "1899-12-30")
    message(sprintf(
      "[%s] '%s': %d value(s) decoded as Excel serial date (origin 1899-12-30).",
      table_label, col_name, length(excel_idx)
    ))
  }

  # Unix timestamps in milliseconds: > 1e10
  unix_idx <- idx_present[is_numeric[idx_present] & x_num[idx_present] > 1e10]
  if (length(unix_idx) > 0L) {
    out[unix_idx] <- as.Date(as.POSIXct(x_num[unix_idx] / 1000, origin = "1970-01-01"))
    message(sprintf(
      "[%s] '%s': %d value(s) decoded as Unix timestamp (ms).",
      table_label, col_name, length(unix_idx)
    ))
  }

  # ---- String parsing for remaining indices --------------------------------
  remaining <- idx_present[is.na(out[idx_present]) & !is_numeric[idx_present]]

  if (length(remaining) > 0L) {
    parsed <- suppressWarnings(
      lubridate::parse_date_time(
        x_chr[remaining],
        orders = c(
          "Y-m-d", "d-m-Y", "m-d-Y",
          "Y/m/d", "d/m/Y", "m/d/Y",
          "Ymd",   "dmY",   "mdY",
          "Y-m-d H:M:S", "d-m-Y H:M:S", "m-d-Y H:M:S"
        ),
        quiet = TRUE
      )
    )
    parsed_date <- as.Date(parsed)
    succeeded   <- !is.na(parsed_date)
    out[remaining[succeeded]] <- parsed_date[succeeded]

    # Reversed YYYYDDMM / DDMMYYYY attempt on still-failing values
    still_bad <- remaining[!succeeded]
    if (length(still_bad) > 0L) {
      candidates <- x_chr[still_bad]
      digits_only <- gsub("[^0-9]", "", candidates)
      swapped <- ifelse(
        nchar(digits_only) == 8L,
        paste0(substr(digits_only, 1, 4), "-",
               substr(digits_only, 7, 8), "-",
               substr(digits_only, 5, 6)),
        NA_character_
      )
      swap_parsed <- suppressWarnings(as.Date(swapped, format = "%Y-%m-%d"))
      swap_ok <- !is.na(swap_parsed)
      if (any(swap_ok)) {
        out[still_bad[swap_ok]] <- swap_parsed[swap_ok]
        message(sprintf(
          "[%s] '%s': %d value(s) decoded via day/month swap (YYYYDDMM -> YYYYMMDD).",
          table_label, col_name, sum(swap_ok)
        ))
        still_bad <- still_bad[!swap_ok]
      }
    }

    # Flag values that look encrypted (long, no recognisable separator pattern)
    if (length(still_bad) > 0L) {
      encrypted_like <- nchar(x_chr[still_bad]) > 12 &
        !grepl("[/\\-]", x_chr[still_bad])
      n_encrypted <- sum(encrypted_like, na.rm = TRUE)
      n_unparsed  <- length(still_bad)

      if (n_encrypted > 0L) {
        warning(sprintf(
          "[%s] '%s': %d value(s) appear encrypted/undecodable (long string, no date separator). They will be set to NA. Sample: %s",
          table_label, col_name, n_encrypted,
          paste(head(x_chr[still_bad][encrypted_like], 3), collapse = " | ")
        ))
      } else if (n_unparsed > 0L) {
        warning(sprintf(
          "[%s] '%s': %d value(s) could not be parsed as dates and will be set to NA. Sample: %s",
          table_label, col_name, n_unparsed,
          paste(head(x_chr[still_bad], 3), collapse = " | ")
        ))
      }
    }
  }

  # Summary
  n_decoded <- sum(!is.na(out[idx_present]))
  n_failed  <- length(idx_present) - n_decoded
  message(sprintf(
    "[%s] '%s': %d / %d non-missing value(s) successfully parsed as Date (%d failed).",
    table_label, col_name, n_decoded, length(idx_present), n_failed
  ))

  out
}


#' Detect and convert all date-like columns in a table
#'
#' Auto-detects columns with date/time-like names (or uses a supplied list),
#' then calls \code{prep_parse_date_column()} on each.  Run this before any join.
#'
#' @param data Data frame.
#' @param cols Character vector of column names to convert. When \code{NULL},
#'   columns whose names match a date/time pattern are detected automatically.
#' @param table_label Character. Label used in messages.
#'
#' @return Data frame with date columns converted to \code{Date}.
#' @export
prep_coerce_dates <- function(data,
                              cols = NULL,
                              table_label = "table") {
  if (is.null(cols)) {
    cols <- grep(
      "(^date_|_date$|date$|_date_|^fever_date|^date_HAI|_treat$|date_treat)",
      names(data),
      value = TRUE, ignore.case = TRUE
    )
  }
  cols <- intersect(cols, names(data))
  if (length(cols) == 0L) {
    return(data)
  }
  for (col in cols) {
    data[[col]] <- prep_parse_date_column(data[[col]],
                                          col_name    = col,
                                          table_label = table_label)
  }
  data
}


#' Run all pre-join sanity checks for one table
#'
#' Convenience wrapper that calls \code{prep_check_columns()},
#' \code{prep_check_keys()}, and \code{prep_coerce_dates()} in sequence.
#' Returns the date-coerced data along with a check report.
#'
#' @param data Data frame.
#' @param required_cols Character vector of required column names.
#' @param key_col Character. Primary/foreign key column to quality-check.
#' @param expected_types Named character vector of expected column classes.
#' @param date_cols Character vector of date column names to coerce. NULL
#'   triggers auto-detection.
#' @param table_label Character. Label used in all messages.
#' @param stop_on_missing Logical. Passed to \code{prep_check_columns()}.
#'
#' @return A list with \code{data} (date-coerced) and \code{report} (check summary).
#' @export
prep_validate_table <- function(data,
                                required_cols  = character(),
                                key_col        = NULL,
                                expected_types = character(),
                                date_cols      = NULL,
                                table_label    = "table",
                                stop_on_missing = TRUE) {
  col_report <- prep_check_columns(
    data            = data,
    required        = required_cols,
    expected_types  = expected_types,
    table_label     = table_label,
    stop_on_missing = stop_on_missing
  )

  key_report <- NULL
  if (!is.null(key_col)) {
    key_report <- prep_check_keys(
      data         = data,
      key_col      = key_col,
      table_label  = table_label
    )
  }

  data <- prep_coerce_dates(data, cols = date_cols, table_label = table_label)

  list(data = data, col_report = col_report, key_report = key_report)
}


# ---------------------------------------------------------------------------
# Section 2: Internal helpers (rename, standardize, clean)
# ---------------------------------------------------------------------------

# Standardise a key column: trim whitespace, map placeholders -> NA
.stwd_standardize_key <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NULL", "null", "NA", "N/A", "None", "none", "nan", "NaN")] <- NA_character_
  x
}

# Apply all known ICMR rename rules for a given table type
.stwd_apply_renames <- function(data, table_type) {
  rename_map <- switch(
    table_type,
    patient = c(
      id            = "PatientInformation_id",
      user          = "user_id",
      center        = "center_id",
      state         = "state_id",
      district      = "district_id",
      quantification = "quantification_id"
    ),
    antibioticinfo = c(
      id               = "id_antibioticinfo",
      PatientInformation = "PatientInformation_id",
      user             = "user_id",
      antibiotic_name_1 = "antibiotic_name_1_id",
      revised_id       = "revised_id_id",
      treat_id         = "treat_id_id"
    ),
    audit = c(
      id      = "id_audit",
      patient = "patient_id",
      user    = "user_id"
    ),
    hai = c(
      id               = "id_hai",
      PatientInformation = "PatientInformation_id",
      user             = "user_id"
    ),
    organisminfo = c(
      id                      = "id_organisminfo",
      PatientInformation      = "PatientInformation_id",
      user                    = "user_id",
      organism                = "organism_id",
      Sample_type1            = "Sample_type1_id",
      site_of_isolation_treat = "site_of_isolation_treat_id",
      organism_treat          = "organism_treat_id"
    ),
    organisminfo_antibiotic = c(
      id                         = "id_organisminfo_antibiotic",
      PatientInformation         = "PatientInformation_id",
      user                       = "user_id",
      organism_info              = "organism_info_id",
      organism_treat_antibiotic  = "organism_treat_antibiotic_id"
    ),
    district   = c(name = "district_name"),
    state      = c(name = "state_name"),
    organisms  = c(name = "organism_name"),
    center     = c(name = "center_name"),
    antibiotics = c(
      id           = "id_antibiotics_lookup",
      date_created = "date_created_antibiotics"
    ),
    c()  # sample_type and unknown types: no renames
  )

  present <- intersect(names(rename_map), names(data))
  if (length(present) > 0L) {
    names(data)[match(present, names(data))] <- unname(rename_map[present])
  }

  # Drop user_id from secondary tables (not the patient table)
  if (table_type %in% c("antibioticinfo", "audit", "hai",
                         "organisminfo", "organisminfo_antibiotic")) {
    data$user_id <- NULL
  }

  data
}

# Replace placeholder strings with NA in specified columns
.stwd_placeholders_to_na <- function(data, cols) {
  cols <- intersect(cols, names(data))
  placeholders <- c("", "NULL", "null", "NA", "N/A", "None", "none", "nan", "NaN")
  for (col in cols) {
    vals <- trimws(as.character(data[[col]]))
    data[[col]][vals %in% placeholders] <- NA
  }
  data
}

# Coalesce a .x/.y column pair and drop the suffixed originals
.stwd_coalesce_pair <- function(data, var) {
  lx <- paste0(var, ".x")
  ly <- paste0(var, ".y")
  if (all(c(lx, ly) %in% names(data))) {
    data[[var]] <- dplyr::coalesce(data[[lx]], data[[ly]])
    data[[lx]] <- NULL
    data[[ly]] <- NULL
  }
  data
}


# ---------------------------------------------------------------------------
# Section 3: Per-centre joining
# ---------------------------------------------------------------------------

#' Join all tables for one stewardship centre
#'
#' Replicates the per-centre joining pattern from Data_merging_all_centres.Rmd
#' as a single reusable function.  Performs pre-join validation and date
#' conversion before every merge, then reconciles duplicate columns afterwards.
#'
#' @section Join order:
#' \enumerate{
#'   \item Rename & validate all input tables.
#'   \item Pre-join \code{organisminfo_antibiotic} into \code{organisminfo}
#'     on \code{id_organisminfo = organism_info_id} (avoids cartesian explosion).
#'   \item Left-join patient <- antibioticinfo, audit, hai, organism bundle.
#'   \item Left-join lookup tables (district+state, sample_type, organisms,
#'     antibiotics, center).
#'   \item Reconcile duplicate columns.
#'   \item Replace placeholders with NA.
#'   \item Optionally normalise antibiotics and create event IDs.
#' }
#'
#' @param patient Patient main table (data frame).
#' @param antibioticinfo Antibiotic info table or NULL.
#' @param audit Audit table or NULL.
#' @param hai HAI table or NULL.
#' @param organisminfo Organism info table or NULL.
#' @param organisminfo_antibiotic Organism-antibiotic table or NULL.
#' @param lookups Named list of lookup tables. Recognised names:
#'   \code{district}, \code{state}, \code{sample_type}, \code{organisms},
#'   \code{antibiotics}, \code{center}.
#' @param centre_name Character label for this centre.
#' @param patient_id_col Patient identifier column. Default \code{"PatientInformation_id"}.
#' @param kmc_organisminfo_col_names Optional character vector of column names to
#'   supply when the organisminfo sheet is missing its header row (KMC quirk).
#' @param normalize_abx Logical. Run \code{prep_standardize_antibiotics()} after joining.
#' @param create_events Logical. Run \code{prep_create_event_ids()} after joining.
#' @param culture_date_col Date column for event creation. Default \code{"date_treat"}.
#' @param specimen_col Specimen column for event creation. Default \code{"sample_type"}.
#' @param organism_col Organism column for event creation. Default \code{"organism_name"}.
#' @param antibiotic_col Antibiotic column. Default \code{"antibiotic_name"}.
#' @param required_analysis_cols Columns used for completeness checks.
#'
#' @return A list with class \code{"stwd_join_result"} containing:
#'   \describe{
#'     \item{\code{data}}{Joined data frame (all rows, pre-filter).}
#'     \item{\code{qc}}{Named list of check reports and merge summaries.}
#'     \item{\code{centre_name}}{Centre label.}
#'   }
#' @export
prep_join_stewardship_centre <- function(
    patient,
    antibioticinfo              = NULL,
    audit                      = NULL,
    hai                        = NULL,
    organisminfo               = NULL,
    organisminfo_antibiotic    = NULL,
    lookups                    = list(),
    centre_name                = NULL,
    patient_id_col             = "PatientInformation_id",
    kmc_organisminfo_col_names = NULL,
    normalize_abx              = TRUE,
    create_events              = TRUE,
    culture_date_col           = "date_treat",
    specimen_col               = "sample_type",
    organism_col               = "organism_name",
    antibiotic_col             = "antibiotic_name",
    required_analysis_cols     = c("final_outcome", "organism_name",
                                   "antibiotic_name", "antibiotic_value")
) {
  lbl <- centre_name %||% "centre"
  message(sprintf("\n========== prep_join_stewardship_centre: %s ==========", lbl))

  # --- Validate & rename patient -------------------------------------------
  patient <- as.data.frame(patient, stringsAsFactors = FALSE)
  patient <- .stwd_apply_renames(patient, "patient")

  if (!patient_id_col %in% names(patient)) {
    stop(sprintf("[%s] Patient table must contain '%s'.", lbl, patient_id_col))
  }

  val_patient <- prep_validate_table(
    data           = patient,
    required_cols  = c(patient_id_col),
    key_col        = patient_id_col,
    table_label    = paste0("patient_", lbl)
  )
  patient     <- val_patient$data
  check_log   <- list(patient = val_patient$col_report)

  patient[[patient_id_col]] <- .stwd_standardize_key(patient[[patient_id_col]])

  # --- Helper: validate + rename a secondary table -------------------------
  .prep_secondary <- function(df, type) {
    if (is.null(df)) return(NULL)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df <- .stwd_apply_renames(df, type)
    v  <- prep_validate_table(df,
                              key_col     = patient_id_col,
                              table_label = paste0(type, "_", lbl))
    check_log[[type]] <<- v$col_report
    v$data
  }

  antibioticinfo          <- .prep_secondary(antibioticinfo,           "antibioticinfo")
  audit                   <- .prep_secondary(audit,                    "audit")
  hai                     <- .prep_secondary(hai,                      "hai")
  organisminfo            <- .prep_secondary(organisminfo,             "organisminfo")
  organisminfo_antibiotic <- .prep_secondary(organisminfo_antibiotic,  "organisminfo_antibiotic")

  # KMC quirk: header-less organisminfo sheet -- column names already applied at read time
  # (col_names passed to read_excel), so no extra action needed here.

  # --- Pre-join: organisminfo + organisminfo_antibiotic --------------------
  merge_log   <- list()
  orginfo_combined <- NULL

  if (!is.null(organisminfo) && !is.null(organisminfo_antibiotic)) {
    message(sprintf("[%s] Pre-joining organisminfo + organisminfo_antibiotic ...", lbl))

    .stwd_check_prejoin_keys <- function() {
      if (!"id_organisminfo" %in% names(organisminfo)) {
        stop(sprintf("[%s] organisminfo must contain 'id_organisminfo'.", lbl))
      }
      if (!"organism_info_id" %in% names(organisminfo_antibiotic)) {
        stop(sprintf("[%s] organisminfo_antibiotic must contain 'organism_info_id'.", lbl))
      }
      prep_check_keys(organisminfo,           "id_organisminfo",  paste0("organisminfo_", lbl))
      prep_check_keys(organisminfo_antibiotic,"organism_info_id", paste0("orgabx_", lbl))
    }
    .stwd_check_prejoin_keys()

    dt_org    <- data.table::as.data.table(organisminfo)
    dt_orgabx <- data.table::as.data.table(organisminfo_antibiotic)

    orginfo_combined <- merge(
      dt_org, dt_orgabx,
      by.x = "id_organisminfo", by.y = "organism_info_id",
      all.x = TRUE, allow.cartesian = TRUE
    )
    orginfo_combined <- as.data.frame(orginfo_combined, stringsAsFactors = FALSE)

    # Resolve duplicate PatientInformation_id
    pid_x <- paste0(patient_id_col, ".x")
    pid_y <- paste0(patient_id_col, ".y")
    if (all(c(pid_x, pid_y) %in% names(orginfo_combined))) {
      orginfo_combined[[patient_id_col]] <- dplyr::coalesce(
        orginfo_combined[[pid_x]],
        orginfo_combined[[pid_y]]
      )
      orginfo_combined[[pid_x]] <- NULL
      orginfo_combined[[pid_y]] <- NULL
    }

    # Coalesce date_treat: organisminfo (.x) preferred over orgabx (.y)
    if (all(c("date_treat.x", "date_treat.y") %in% names(orginfo_combined))) {
      orginfo_combined$date_treat.x <- prep_parse_date_column(
        orginfo_combined$date_treat.x, "date_treat.x", lbl)
      orginfo_combined$date_treat.y <- prep_parse_date_column(
        orginfo_combined$date_treat.y, "date_treat.y", lbl)
      orginfo_combined$date_treat <- dplyr::coalesce(
        orginfo_combined$date_treat.x,
        orginfo_combined$date_treat.y
      )
      orginfo_combined$date_treat.x <- NULL
      orginfo_combined$date_treat.y <- NULL
    }

    merge_log$orginfo_prejoin <- data.frame(
      step = "orginfo_prejoin",
      rows_orginfo    = nrow(organisminfo),
      rows_orgabx     = nrow(organisminfo_antibiotic),
      rows_combined   = nrow(orginfo_combined),
      stringsAsFactors = FALSE
    )
    message(sprintf("[%s] Organism pre-join: %d + %d -> %d rows.",
                    lbl, nrow(organisminfo), nrow(organisminfo_antibiotic),
                    nrow(orginfo_combined)))

  } else if (!is.null(organisminfo)) {
    orginfo_combined <- organisminfo
  }

  # --- Sequential patient-level joins --------------------------------------
  .left_merge <- function(left, right, by_l, by_r, step_name) {
    if (is.null(right)) return(left)
    rows_before <- nrow(left)
    right[[by_r]] <- .stwd_standardize_key(right[[by_r]])
    left[[by_l]]  <- .stwd_standardize_key(left[[by_l]])

    prep_check_keys(right, by_r, paste0(step_name, "_right_", lbl))

    out <- merge(
      data.table::as.data.table(left),
      data.table::as.data.table(right),
      by.x = by_l, by.y = by_r,
      all.x = TRUE, allow.cartesian = TRUE
    )
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    merge_log[[step_name]] <<- data.frame(
      step         = step_name,
      rows_before  = rows_before,
      rows_right   = nrow(right),
      rows_after   = nrow(out),
      stringsAsFactors = FALSE
    )
    message(sprintf("[%s] %-30s: %d -> %d rows (right: %d).",
                    lbl, step_name, rows_before, nrow(out), nrow(right)))
    out
  }

  patient <- .left_merge(patient, antibioticinfo,   patient_id_col, patient_id_col, "join_antibioticinfo")
  patient <- .left_merge(patient, audit,            patient_id_col, "patient_id",   "join_audit")
  patient <- .left_merge(patient, hai,              patient_id_col, patient_id_col, "join_hai")
  patient <- .left_merge(patient, orginfo_combined, patient_id_col, patient_id_col, "join_organism_bundle")

  # --- Lookup joins --------------------------------------------------------
  district  <- lookups$district
  state     <- lookups$state
  sample_tp <- lookups$sample_type
  organisms <- lookups$organisms
  antibiotics <- lookups$antibiotics
  center    <- lookups$center

  # Rename lookup tables
  if (!is.null(state))      state      <- .stwd_apply_renames(as.data.frame(state,      stringsAsFactors = FALSE), "state")
  if (!is.null(district))   district   <- .stwd_apply_renames(as.data.frame(district,   stringsAsFactors = FALSE), "district")
  if (!is.null(sample_tp))  sample_tp  <- .stwd_apply_renames(as.data.frame(sample_tp,  stringsAsFactors = FALSE), "sample_type")
  if (!is.null(organisms))  organisms  <- .stwd_apply_renames(as.data.frame(organisms,  stringsAsFactors = FALSE), "organisms")
  if (!is.null(antibiotics)) antibiotics <- .stwd_apply_renames(as.data.frame(antibiotics, stringsAsFactors = FALSE), "antibiotics")
  if (!is.null(center))     center     <- .stwd_apply_renames(as.data.frame(center,     stringsAsFactors = FALSE), "center")

  # Standardize lookup IDs to character
  .as_chr_id <- function(df, col) { if (!is.null(df) && col %in% names(df)) df[[col]] <- as.character(df[[col]]); df }
  state      <- .as_chr_id(state,      "id")
  district   <- .as_chr_id(district,   "id")
  district   <- .as_chr_id(district,   "state_id")
  sample_tp  <- .as_chr_id(sample_tp,  "id")
  organisms  <- .as_chr_id(organisms,  "id")
  center     <- .as_chr_id(center,     "id")

  # Build district_with_state
  district_with_state <- NULL
  if (!is.null(district) && !is.null(state) &&
      "state_id" %in% names(district) && "id" %in% names(state)) {
    district_with_state <- merge(
      data.table::as.data.table(district),
      data.table::as.data.table(state),
      by.x = "state_id", by.y = "id", all.x = TRUE
    )
    district_with_state <- as.data.frame(district_with_state, stringsAsFactors = FALSE)
    message(sprintf("[%s] district+state lookup built: %d rows.", lbl, nrow(district_with_state)))
  }

  # Standardize patient FK columns to character before lookup joins
  for (fk in c("district_id", "Sample_type1_id", "organism_id", "center_id")) {
    if (fk %in% names(patient)) {
      patient[[fk]] <- as.character(patient[[fk]])
    }
  }

  if (!is.null(district_with_state) && "district_id" %in% names(patient) && "id" %in% names(district_with_state)) {
    patient <- .left_merge(patient, district_with_state, "district_id", "id", "lookup_district_state")
    patient <- .stwd_coalesce_pair(patient, "state_id")
  }

  if (!is.null(sample_tp) && "Sample_type1_id" %in% names(patient) && "id" %in% names(sample_tp)) {
    patient <- .left_merge(patient, sample_tp, "Sample_type1_id", "id", "lookup_sample_type")
  }

  if (!is.null(organisms) && "organism_id" %in% names(patient) && "id" %in% names(organisms)) {
    patient <- .left_merge(patient, organisms, "organism_id", "id", "lookup_organisms")
  }

  if (!is.null(antibiotics) && antibiotic_col %in% names(patient) && "name" %in% names(antibiotics)) {
    antibiotics$name <- trimws(as.character(antibiotics$name))
    patient[[antibiotic_col]] <- trimws(as.character(patient[[antibiotic_col]]))
    patient <- .left_merge(patient, antibiotics, antibiotic_col, "name", "lookup_antibiotics")
  }

  if (!is.null(center) && "center_id" %in% names(patient) && "id" %in% names(center)) {
    patient <- .left_merge(patient, center, "center_id", "id", "lookup_center")
  }

  # --- Reconcile duplicate columns -----------------------------------------
  patient <- .stwd_coalesce_pair(patient, "type_of_infection")
  patient <- .stwd_coalesce_pair(patient, "date_of_antibiotics")

  # Drop duplicated date_created columns
  for (dc in c("date_created.x", "date_created.y")) {
    patient[[dc]] <- NULL
  }

  # --- Replace placeholders with NA ----------------------------------------
  patient <- .stwd_placeholders_to_na(
    patient,
    cols = intersect(required_analysis_cols, names(patient))
  )

  # Add centre label if not already present
  if (!"center_name" %in% names(patient) && !is.null(centre_name)) {
    patient$center_name <- centre_name
  }

  # --- Optional: normalize antibiotics ------------------------------------
  if (normalize_abx && antibiotic_col %in% names(patient)) {
    patient <- tryCatch(
      prep_standardize_antibiotics(patient, antibiotic_col = antibiotic_col,
                                   add_class = TRUE, add_aware = TRUE),
      error = function(e) {
        warning(sprintf("[%s] prep_standardize_antibiotics() failed: %s", lbl, conditionMessage(e)))
        patient
      }
    )
  }

  # --- Optional: create event IDs -----------------------------------------
  if (create_events &&
      all(c(patient_id_col, culture_date_col, organism_col) %in% names(patient))) {
    patient <- tryCatch(
      prep_create_event_ids(
        patient,
        patient_col  = patient_id_col,
        date_col     = culture_date_col,
        organism_col = organism_col,
        specimen_col = specimen_col,
        antibiotic_col = antibiotic_col
      ),
      error = function(e) {
        warning(sprintf("[%s] prep_create_event_ids() failed: %s", lbl, conditionMessage(e)))
        patient
      }
    )
  }

  # --- Centre summary ------------------------------------------------------
  n_pts  <- dplyr::n_distinct(patient[[patient_id_col]], na.rm = TRUE)
  message(sprintf(
    "[%s] Join complete: %d rows | %d unique patients | %d columns.",
    lbl, nrow(patient), n_pts, ncol(patient)
  ))

  out <- list(
    data        = patient,
    qc          = list(col_checks = check_log, merge_log = dplyr::bind_rows(merge_log)),
    centre_name = centre_name %||% NA_character_
  )
  class(out) <- c("stwd_join_result", "list")
  out
}


# ---------------------------------------------------------------------------
# Section 4: Post-join cleaning & analysis-ready filter
# ---------------------------------------------------------------------------

#' Filter joined centre data to analysis-ready patients
#'
#' Applies the attrition logic from Data_merging_all_centres.Rmd:
#'   1. Replace placeholder strings with NA in analytical columns.
#'   2. Exclude fungal isolates (optional).
#'   3. Keep only patients for whom all required analytical columns have at
#'      least one non-NA value across their rows.
#'
#' @param data Joined data frame (output of \code{prep_join_stewardship_centre()$data}).
#' @param patient_id_col Patient identifier column.
#' @param required_cols Required analytical columns.
#' @param organism_group_col Column identifying organism groups.
#' @param fungal_label Label used for fungal isolates in \code{organism_group_col}.
#' @param verbose Logical. Print attrition counts. Default TRUE.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{data}}{Filtered data frame.}
#'     \item{\code{patient_flags}}{Per-patient Present/NA flag table.}
#'     \item{\code{combo_counts}}{Unique flag combinations and patient counts.}
#'     \item{\code{attrition}}{Attrition summary tibble.}
#'   }
#' @export
prep_filter_analysis_ready <- function(
    data,
    patient_id_col  = "PatientInformation_id",
    required_cols   = c("final_outcome", "organism_name", "antibiotic_name", "antibiotic_value"),
    organism_group_col = "org_group",
    fungal_label    = "Fungal isolates",
    verbose         = TRUE
) {
  if (!patient_id_col %in% names(data)) {
    stop(sprintf("Missing patient column '%s'.", patient_id_col))
  }

  data <- .stwd_placeholders_to_na(data, cols = intersect(required_cols, names(data)))

  attrition <- tibble::tibble(
    step            = "input",
    rows            = nrow(data),
    unique_patients = dplyr::n_distinct(data[[patient_id_col]], na.rm = TRUE)
  )

  # Fungal exclusion
  n_fungal_pts <- 0L
  if (organism_group_col %in% names(data)) {
    fungal_ids <- unique(data[[patient_id_col]][
      !is.na(data[[organism_group_col]]) & data[[organism_group_col]] == fungal_label
    ])
    n_fungal_pts <- length(fungal_ids)
    data <- data[is.na(data[[organism_group_col]]) |
                   data[[organism_group_col]] != fungal_label, , drop = FALSE]
    attrition <- dplyr::bind_rows(
      attrition,
      tibble::tibble(
        step            = "after_fungal_exclusion",
        rows            = nrow(data),
        unique_patients = dplyr::n_distinct(data[[patient_id_col]], na.rm = TRUE)
      )
    )
    if (verbose) {
      message(sprintf("Fungal exclusion: removed %d unique patient(s).", n_fungal_pts))
    }
  }

  # Patient-level completeness flags
  used_cols <- intersect(required_cols, names(data))
  if (length(used_cols) == 0L) {
    stop("None of the required analysis columns were found in `data`.")
  }

  flags <- data %>%
    dplyr::group_by(.data[[patient_id_col]]) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(used_cols),
        ~ ifelse(all(is.na(.x)), "NA", "Present"),
        .names = "{.col}_flag"
      ),
      .groups = "drop"
    )

  combo_counts <- flags %>%
    dplyr::count(dplyr::across(dplyr::everything()),
                 name = "unique_patients", sort = TRUE)

  flag_cols <- paste0(used_cols, "_flag")
  keep_ids  <- flags
  for (fc in flag_cols) {
    keep_ids <- keep_ids[keep_ids[[fc]] == "Present", , drop = FALSE]
  }

  filtered <- data[data[[patient_id_col]] %in% keep_ids[[patient_id_col]], , drop = FALSE]

  attrition <- dplyr::bind_rows(
    attrition,
    tibble::tibble(
      step            = "analysis_ready",
      rows            = nrow(filtered),
      unique_patients = dplyr::n_distinct(filtered[[patient_id_col]], na.rm = TRUE)
    )
  )

  if (verbose) {
    message(sprintf(
      "Analysis-ready filter: %d -> %d rows | %d -> %d unique patients.",
      nrow(data), nrow(filtered),
      dplyr::n_distinct(data[[patient_id_col]], na.rm = TRUE),
      dplyr::n_distinct(filtered[[patient_id_col]], na.rm = TRUE)
    ))
    print(attrition)
  }

  list(
    data          = filtered,
    patient_flags = flags,
    combo_counts  = combo_counts,
    attrition     = attrition
  )
}


# ---------------------------------------------------------------------------
# Section 5: Multi-centre binding
# ---------------------------------------------------------------------------

#' Bind multiple stewardship centre results into one dataset
#'
#' Accepts a named list of \code{stwd_join_result} objects, plain data frames,
#' or \code{filter_icmr_analysis_ready} results.  Harmonises columns,
#' re-coerces dates, and produces a cross-centre QC summary.
#'
#' @param centres Named list. Each element is either a \code{stwd_join_result},
#'   a list with a \code{$data} element, or a plain data frame.
#' @param keep_common_only Logical. Keep only columns present in all centres.
#' @param patient_id_col Patient identifier column.
#' @param verbose Logical. Print binding summary. Default TRUE.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{data}}{Combined data frame.}
#'     \item{\code{qc}}{Cross-centre QC: centre summary, common columns,
#'       duplicate patient IDs across centres.}
#'   }
#' @export
prep_bind_stewardship_centres <- function(
    centres,
    keep_common_only = TRUE,
    patient_id_col   = "PatientInformation_id",
    verbose          = TRUE
) {
  if (!is.list(centres) || length(centres) == 0L) {
    stop("`centres` must be a non-empty named list.")
  }

  centre_names <- names(centres)
  if (is.null(centre_names) || any(centre_names == "")) {
    centre_names <- paste0("centre_", seq_along(centres))
  }

  dfs <- lapply(seq_along(centres), function(i) {
    obj <- centres[[i]]
    df  <- if (inherits(obj, "stwd_join_result")) {
      obj$data
    } else if (is.list(obj) && "data" %in% names(obj) && is.data.frame(obj$data)) {
      obj$data
    } else if (is.data.frame(obj)) {
      obj
    } else {
      stop(sprintf("Element '%s' must be a data frame or a stwd_join_result.",
                   centre_names[[i]]))
    }
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!"center_name" %in% names(df)) {
      df$center_name <- centre_names[[i]]
    }
    df
  })

  common_cols <- Reduce(intersect, lapply(dfs, names))
  if (verbose) {
    message(sprintf("Common columns across %d centre(s): %d",
                    length(dfs), length(common_cols)))
  }

  if (keep_common_only) {
    dfs <- lapply(dfs, function(df) df[, common_cols, drop = FALSE])
  }

  dfs <- lapply(seq_along(dfs), function(i) {
    prep_coerce_dates(dfs[[i]], table_label = centre_names[[i]])
  })

  combined <- dplyr::bind_rows(dfs)

  centre_summary <- dplyr::bind_rows(lapply(seq_along(dfs), function(i) {
    df <- dfs[[i]]
    tibble::tibble(
      center_name              = centre_names[[i]],
      rows                     = nrow(df),
      unique_patients          = if (patient_id_col %in% names(df)) dplyr::n_distinct(df[[patient_id_col]], na.rm = TRUE) else NA_integer_,
      missing_organism         = if ("organism_name"    %in% names(df)) sum(is.na(df$organism_name))    else NA_integer_,
      missing_antibiotic       = if ("antibiotic_name"  %in% names(df)) sum(is.na(df$antibiotic_name))  else NA_integer_,
      missing_antibiotic_value = if ("antibiotic_value" %in% names(df)) sum(is.na(df$antibiotic_value)) else NA_integer_,
      missing_outcome          = if ("final_outcome"    %in% names(df)) sum(is.na(df$final_outcome))    else NA_integer_
    )
  }))

  duplicate_patients <- NULL
  if (patient_id_col %in% names(combined) && "center_name" %in% names(combined)) {
    duplicate_patients <- combined %>%
      dplyr::distinct(.data[[patient_id_col]], .data$center_name) %>%
      dplyr::count(.data[[patient_id_col]], name = "n_centres") %>%
      dplyr::filter(.data$n_centres > 1L)
    if (nrow(duplicate_patients) > 0L) {
      warning(sprintf(
        "%d patient ID(s) appear in more than one centre. Check `qc$duplicate_patients`.",
        nrow(duplicate_patients)
      ))
    }
  }

  if (verbose) {
    message(sprintf("Combined dataset: %d rows | %d columns.", nrow(combined), ncol(combined)))
    print(centre_summary)
  }

  list(
    data = combined,
    qc   = list(
      centre_summary             = centre_summary,
      common_columns             = common_cols,
      duplicate_patients_across_centres = duplicate_patients
    )
  )
}


# ---------------------------------------------------------------------------
# Section 6: Save helper
# ---------------------------------------------------------------------------

#' Validate and save a merged centre dataset
#'
#' Prints a summary then writes to CSV via \code{data.table::fwrite}.
#'
#' @param data Data frame to save.
#' @param path File path for the output CSV.
#' @param centre_name Character label used in the console summary.
#' @param patient_id_col Patient identifier column.
#'
#' @return Invisibly returns \code{data}.
#' @export
prep_save_centre <- function(data,
                             path,
                             centre_name    = "centre",
                             patient_id_col = "PatientInformation_id") {
  message(sprintf("\n=== %s - Save Summary ===", centre_name))
  message(sprintf("  Rows: %d  |  Cols: %d", nrow(data), ncol(data)))
  if (patient_id_col %in% names(data)) {
    message(sprintf("  Unique patients: %d",
                    data.table::uniqueN(data[[patient_id_col]])))
  }
  if ("organism_name"  %in% names(data)) message(sprintf("  NA organism_name:  %d", sum(is.na(data$organism_name))))
  if ("antibiotic_name" %in% names(data)) message(sprintf("  NA antibiotic_name: %d", sum(is.na(data$antibiotic_name))))

  t <- system.time(data.table::fwrite(data, path))
  file_mb <- round(file.size(path) / 1024^2, 1)
  message(sprintf("  Saved: %s  |  %.1f MB  |  %.2f sec", basename(path), file_mb, t[["elapsed"]]))
  message(sprintf("  Path: %s", path))

  invisible(data)
}


# ---------------------------------------------------------------------------
# Section 7: S3 print method
# ---------------------------------------------------------------------------

#' @export
print.stwd_join_result <- function(x, ...) {
  cat("Stewardship join result\n")
  cat(sprintf("  Centre     : %s\n", x$centre_name %||% "NA"))
  cat(sprintf("  Rows       : %d\n", nrow(x$data)))
  if ("PatientInformation_id" %in% names(x$data)) {
    cat(sprintf("  Unique pts : %d\n",
                dplyr::n_distinct(x$data$PatientInformation_id, na.rm = TRUE)))
  }
  if (!is.null(x$qc$merge_log) && nrow(x$qc$merge_log) > 0L) {
    cat("  Merge log  :\n")
    print(x$qc$merge_log[, c("step", "rows_before", "rows_after")], row.names = FALSE)
  }
  invisible(x)
}
