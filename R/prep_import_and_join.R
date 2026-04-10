# prep_import_and_join.R
# Data import, pivot, and wide-format functions for AMR data

#' Convert Wide Format to Long Format
#'
#' Converts wide format data (where each antibiotic is a column) to long format
#' (one row per organism-antibiotic combination). This is the first step before
#' normalization and analysis.
#'
#' @param data Data frame in wide format (antibiotics as columns)
#' @param antibiotic_cols Character vector. Names of antibiotic columns to pivot.
#'   If NULL, will auto-detect based on pattern. Default NULL.
#' @param pattern Character. Regex pattern to identify antibiotic columns if
#'   antibiotic_cols not provided. Default NULL (no auto-detection).
#' @param id_cols Character vector. Columns to keep as identifiers (not pivoted).
#'   Default c("patient_id", "event_id", "organism_name", "date_of_culture").
#' @param antibiotic_name_col Character. Name for the new column containing
#'   antibiotic names. Default "antibiotic_name".
#' @param antibiotic_value_col Character. Name for the new column containing
#'   susceptibility results. Default "antibiotic_value".
#' @param remove_missing Logical. Remove rows where antibiotic_value is NA, empty,
#'   or "-". Default TRUE.
#' @param create_event_id Logical. Create event_id column if it doesn't exist
#'   (uses row numbers). Default FALSE.
#'
#' @return Data frame in long format
#' @export
#'
#' @examples
#' \dontrun{
#' # Specify antibiotic columns explicitly
#' long_data <- prep_pivot_ast_wide_to_long(
#'   data = raw_data,
#'   antibiotic_cols = c("AMIKACIN", "GENTAMICIN", "CIPROFLOXACIN")
#' )
#'
#' # Auto-detect columns by pattern (columns 12-53)
#' long_data <- prep_pivot_ast_wide_to_long(
#'   data = raw_data,
#'   antibiotic_cols = names(raw_data)[12:53]
#' )
#'
#' # Auto-detect uppercase antibiotic names
#' long_data <- prep_pivot_ast_wide_to_long(
#'   data = raw_data,
#'   pattern = "^[A-Z]+$"
#' )
#' }
prep_pivot_ast_wide_to_long <- function(data,
                                        antibiotic_cols = NULL,
                                        pattern = NULL,
                                        id_cols = c("patient_id", "event_id", "organism_name", "date_of_culture"),
                                        antibiotic_name_col = "antibiotic_name",
                                        antibiotic_value_col = "antibiotic_value",
                                        remove_missing = TRUE,
                                        create_event_id = FALSE) {
  n_before <- nrow(data)

  # Create event_id if needed
  if (create_event_id && !"event_id" %in% names(data)) {
    message("Creating event_id column from row numbers...")
    data$event_id <- seq_len(nrow(data))
  }

  # Auto-detect antibiotic columns if not provided
  if (is.null(antibiotic_cols)) {
    if (!is.null(pattern)) {
      # Use pattern to detect
      antibiotic_cols <- names(data)[grepl(pattern, names(data))]
      message(sprintf(
        "Auto-detected %d antibiotic columns using pattern '%s'",
        length(antibiotic_cols), pattern
      ))
    } else {
      stop("Either antibiotic_cols or pattern must be provided")
    }
  }

  if (length(antibiotic_cols) == 0) {
    stop("No antibiotic columns found")
  }

  # Filter id_cols to only existing columns
  id_cols <- intersect(id_cols, names(data))

  message(sprintf(
    "Pivoting %d antibiotic columns to long format...",
    length(antibiotic_cols)
  ))

  # Check if all antibiotic columns exist
  missing_cols <- setdiff(antibiotic_cols, names(data))
  if (length(missing_cols) > 0) {
    warning(sprintf(
      "Some antibiotic columns not found: %s",
      paste(missing_cols, collapse = ", ")
    ))
    antibiotic_cols <- intersect(antibiotic_cols, names(data))
  }

  # Pivot to long format
  long_data <- data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(antibiotic_cols),
      names_to = antibiotic_name_col,
      values_to = antibiotic_value_col
    )

  n_after <- nrow(long_data)

  message(sprintf(
    "Pivoted: %d rows -> %d rows (wide -> long)",
    n_before, n_after
  ))

  # Remove missing values if requested
  if (remove_missing) {
    n_before_filter <- nrow(long_data)

    long_data <- long_data %>%
      dplyr::filter(
        !is.na(!!rlang::sym(antibiotic_value_col)),
        !!rlang::sym(antibiotic_value_col) != "",
        !!rlang::sym(antibiotic_value_col) != "-"
      )

    n_removed <- n_before_filter - nrow(long_data)

    if (n_removed > 0) {
      message(sprintf(
        "Removed %d rows with missing/empty values",
        n_removed
      ))
    }
  }

  message(sprintf(
    "Final long format: %d rows x %d columns",
    nrow(long_data), ncol(long_data)
  ))

  # Summary
  n_unique_abx <- dplyr::n_distinct(long_data[[antibiotic_name_col]])
  n_unique_events <- dplyr::n_distinct(long_data[["event_id"]], na.rm = TRUE)

  message(sprintf(
    "Contains: %d unique antibiotics across %d events",
    n_unique_abx, n_unique_events
  ))

  return(long_data)
}


#' Create Wide Format Dataset
#'
#' Converts long format (one row per organism-antibiotic) to wide format
#' (one row per event with antibiotic columns). Useful for analysis and
#' machine learning applications.
#'
#' @param data Data frame in long format
#' @param event_col Character. Event ID column. Default "event_id".
#' @param antibiotic_col Character. Antibiotic/class column to pivot.
#'   Default "antibiotic_normalized".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param prefix Character. Prefix for pivoted columns. Default "abx_".
#' @param keep_cols Character vector. Additional columns to keep from original data.
#'   Default c("patient_id", "organism_normalized", "date_of_culture").
#'
#' @return Wide format data frame (one row per event)
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic wide format
#' wide_data <- prep_create_wide_ast_matrix(data)
#'
#' # Class-level wide format
#' wide_data <- prep_create_wide_ast_matrix(
#'   data,
#'   antibiotic_col = "antibiotic_class",
#'   prefix = "class_"
#' )
#' }
prep_create_wide_ast_matrix <- function(data,
                                        event_col = "event_id",
                                        antibiotic_col = "antibiotic_normalized",
                                        susceptibility_col = "antibiotic_value",
                                        prefix = "abx_",
                                        keep_cols = c("patient_id", "organism_normalized", "date_of_culture")) {
  # Validate required columns
  if (!event_col %in% names(data)) {
    stop(sprintf("Column '%s' not found", event_col))
  }
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Column '%s' not found", antibiotic_col))
  }
  if (!susceptibility_col %in% names(data)) {
    stop(sprintf("Column '%s' not found", susceptibility_col))
  }

  # Filter keep_cols to only existing columns
  keep_cols <- intersect(keep_cols, names(data))
  keep_cols <- unique(c(event_col, keep_cols))

  message(sprintf(
    "Creating wide format: pivoting '%s' column...",
    antibiotic_col
  ))

  n_before <- nrow(data)
  n_events_before <- dplyr::n_distinct(data[[event_col]])

  # Select relevant columns
  data_subset <- data %>%
    dplyr::select(
      dplyr::all_of(keep_cols),
      !!rlang::sym(antibiotic_col),
      !!rlang::sym(susceptibility_col)
    )

  # Pivot to wide format
  wide_data <- data_subset %>%
    # Handle duplicates: take first value
    dplyr::group_by(
      !!rlang::sym(event_col),
      !!rlang::sym(antibiotic_col)
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    # Pivot
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(keep_cols),
      names_from = !!rlang::sym(antibiotic_col),
      values_from = !!rlang::sym(susceptibility_col),
      names_prefix = prefix
    )

  # Clean column names (replace spaces/special chars with underscores)
  names(wide_data) <- gsub("[^A-Za-z0-9_]", "_", names(wide_data))
  names(wide_data) <- gsub("_{2,}", "_", names(wide_data)) # Remove multiple underscores

  n_after <- nrow(wide_data)
  n_antibiotics <- ncol(wide_data) - length(keep_cols)

  message(sprintf(
    "Created wide format: %d events x %d antibiotics",
    n_after, n_antibiotics
  ))

  # Check for data loss
  if (n_after != n_events_before) {
    warning(sprintf(
      "[!] Row count mismatch: %d events in original data, %d rows in wide format",
      n_events_before, n_after
    ))
  }

  return(wide_data)
}
