# daly_pathogen_fraction.R
# Pathogen fraction and CFR calculations for DALY burden estimation

daly_calc_pathogen_fraction_nonfatal <- function(data,
                                 syndrome_col,
                                 syndrome_name,
                                 specimen_col = NULL,
                                 specimen_name = NULL,
                                 polymicrobial_col,
                                 patient_col,
                                 pathogen_col,
                                 outcome_col,
                                 discharged_value = "Discharged",
                                 glass_ref = NULL,
                                 facility_col = NULL,
                                 facility_name = NULL,
                                 pathogen_name = NULL) {
  # -- Input validation ------------------------------------------------------
  if (xor(is.null(specimen_col), is.null(specimen_name))) {
    stop("specimen_col and specimen_name must both be provided or both be NULL.")
  }

  required_cols <- c(
    syndrome_col, specimen_col, polymicrobial_col,
    patient_col, pathogen_col, outcome_col
  )
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  # -- Step 1: Filter syndrome + specimen (optional) + non-fatal -------------
  df <- data %>%
    dplyr::filter(
      .data[[syndrome_col]] == syndrome_name,
      .data[[outcome_col]] == discharged_value
    )
  if (!is.null(specimen_col)) {
    df <- df %>%
      dplyr::filter(.data[[specimen_col]] == specimen_name)
  }
  if (nrow(df) == 0) {
    stop("No non-fatal records remain after syndrome/specimen filtering.")
  }

  # -- Step 2: Optional single-facility restriction ---------------------------
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No records found for facility '%s'.", facility_name))
    }
  }

  # -- Step 3: Optional pathogen filter -------------------------------------
  if (!is.null(pathogen_name)) {
    df <- df %>% dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(df) == 0) {
      stop(sprintf(
        "No records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  # -- Step 4: GLASS filter -- polymicrobial patients only -------------------
  # Monomicrobial patients (polymicrobial_col == 0) are never filtered.
  if (!is.null(glass_ref)) {
    if (is.data.frame(glass_ref)) {
      valid_pathogens <- glass_ref %>%
        dplyr::filter(.data[["specimen"]] == specimen_name) %>%
        dplyr::pull(.data[["pathogen"]])
    } else {
      valid_pathogens <- glass_ref
    }

    df_mono <- df %>% dplyr::filter(.data[[polymicrobial_col]] == 0)
    df_poly <- df %>%
      dplyr::filter(.data[[polymicrobial_col]] == 1) %>%
      dplyr::filter(.data[[pathogen_col]] %in% valid_pathogens)

    n_removed <- sum(df[[polymicrobial_col]] == 1) - nrow(df_poly)
    if (n_removed > 0) {
      message(sprintf(
        "GLASS filter: removed %d row(s) from polymicrobial patients (pathogen not on GLASS list for '%s').",
        n_removed, specimen_name
      ))
    }
    df <- dplyr::bind_rows(df_mono, df_poly)
    if (nrow(df) == 0) {
      stop("No records remain after GLASS reference filtering.")
    }
  }

  # -- Step 5: Deduplicate to one row per patient x pathogen ----------------
  # The raw data has one row per antibiotic tested. Collapse to patient-pathogen
  # level first so that each patient-pathogen pair contributes exactly once.
  group_cols <- if (!is.null(facility_col)) c(facility_col, patient_col) else patient_col

  patient_pathogen <- df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group_cols, pathogen_col))))

  # -- Step 5b: Fractional weight 1/m_r -------------------------------------
  # m_r = distinct valid pathogens for patient r (within facility if relevant).
  patient_pathogen <- patient_pathogen %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      m_r    = dplyr::n_distinct(.data[[pathogen_col]]),
      weight = 1 / m_r
    ) %>%
    dplyr::ungroup()

  # -- Step 6: N^NF_LK = weighted sum per (facility, pathogen) --------------
  agg_cols <- if (!is.null(facility_col)) c(facility_col, pathogen_col) else pathogen_col

  N_LK <- patient_pathogen %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) %>%
    dplyr::summarise(N_NF_LK = sum(weight, na.rm = TRUE), .groups = "drop")

  # -- Step 7: N^NF_L = unique non-fatal patients per facility --------------
  fac_grp <- if (!is.null(facility_col)) facility_col else character(0)

  N_L <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(fac_grp))) %>%
    dplyr::summarise(
      N_NF_L = dplyr::n_distinct(.data[[patient_col]]),
      .groups = "drop"
    )

  # -- Step 8: Compute P'_LK and return -------------------------------------
  if (!is.null(facility_col)) {
    facility_level <- dplyr::left_join(N_LK, N_L, by = facility_col) %>%
      dplyr::mutate(P_Lk_prime = N_NF_LK / N_NF_L)

    # Pooled: sum numerators / sum denominators (not average of ratios)
    pooled <- N_LK %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(N_NF_LK = sum(N_NF_LK), .groups = "drop") %>%
      dplyr::mutate(
        N_NF_L     = sum(N_L$N_NF_L),
        P_Lk_prime = N_NF_LK / N_NF_L
      )

    message(sprintf(
      "P'LK computed: %d pathogen(s) across %d facility/facilities; pooled P'LK also returned.",
      dplyr::n_distinct(df[[pathogen_col]]),
      dplyr::n_distinct(df[[facility_col]])
    ))
    return(list(P_Lk_prime = pooled, facility_level = facility_level))
  } else {
    N_total <- N_L$N_NF_L
    result <- N_LK %>%
      dplyr::mutate(
        N_NF_L     = N_total,
        P_Lk_prime = N_NF_LK / N_NF_L
      )
    message(sprintf(
      "P'LK computed: %d pathogen(s), N^NF_L = %d patient(s).",
      nrow(result), N_total
    ))
    return(list(P_Lk_prime = result))
  }
}


# -- P_LK : Fatal pathogen distribution ----------------------------------------

#' Calculate fatal pathogen distribution (P_\{Lk\})
#'
#' Computes the fatal pathogen distribution for a given infectious syndrome
#' (L) using facility-level microbiology data. This quantity represents the
#' fractional contribution of each pathogen (k) to \strong{fatal} infection
#' cases, and is used in YLL / mortality burden estimation.
#'
#' The unit of analysis is the \strong{patient}. Each patient contributes
#' total weight 1, distributed equally across their valid pathogens.
#' For a patient with \eqn{m_r} valid pathogens, each pathogen receives
#' weight \eqn{1/m_r}.
#'
#' For \strong{polymicrobial patients} (\code{polymicrobial_col == 1}),
#' only pathogens listed in the GLASS reference (\code{glass_ref}) for the
#' given specimen type are retained before weighting. Monomicrobial patients
#' (\code{polymicrobial_col == 0}) are never filtered.
#'
#' The pooled formula across facilities is:
#' \deqn{P_{LK}^{\text{pooled}} =
#'   \frac{\sum_f N^{F}_{f,L,K}}{\sum_f N^{F}_{f,L}}}
#'
#' @param data            Data frame of facility-level microbiology records.
#' @param syndrome_col    Character. Column containing infectious syndrome labels (L).
#' @param syndrome_name   Character. Syndrome to analyse.
#' @param specimen_col    Character or NULL. Column containing specimen type.
#'   If \code{NULL}, no specimen filter is applied.
#' @param specimen_name   Character or NULL. Specimen to restrict to (e.g.,
#'   \code{"Blood"}). Required when \code{specimen_col} is provided.
#' @param polymicrobial_col Character. Column flagging polymicrobial patients
#'   (1 = polymicrobial, 0 = monomicrobial).
#' @param patient_col     Character. Unique patient identifier column.
#' @param pathogen_col    Character. Pathogen (organism) column (k).
#' @param outcome_col     Character. Final patient outcome column.
#' @param death_value     Character. Value indicating a fatal outcome.
#'   Default \code{"Death"}.
#' @param glass_ref       Character vector of valid pathogen names, or a data
#'   frame with columns \code{specimen} and \code{pathogen}. Applied to
#'   polymicrobial patients only. When \code{specimen_name} is \code{NULL} and
#'   \code{glass_ref} is a data frame, all pathogens in the reference are used
#'   regardless of specimen. \code{NULL} skips GLASS filtering.
#' @param facility_col    Character or NULL. Facility identifier column. When
#'   provided without \code{facility_name}, returns both facility-level and
#'   pooled P_LK.
#' @param facility_name   Character or NULL. Restricts to a single facility.
#' @param pathogen_name   Character vector or NULL. Filter to specific pathogen(s).
#'
#' @return A list:
#'   \describe{
#'     \item{P_Lk_fatal}{Pooled P_LK data frame: \code{pathogen_col},
#'       \code{N_F_LK}, \code{N_F_L}, \code{P_Lk_fatal}.}
#'     \item{facility_level}{Per-facility P_LK (only when \code{facility_col}
#'       is supplied and \code{facility_name} is NULL).}
#'   }
#'
#' @export
daly_calc_pathogen_fraction_fatal <- function(data,
                                 syndrome_col,
                                 syndrome_name,
                                 specimen_col = NULL,
                                 specimen_name = NULL,
                                 polymicrobial_col,
                                 patient_col,
                                 pathogen_col,
                                 outcome_col,
                                 death_value = "Death",
                                 glass_ref = NULL,
                                 facility_col = NULL,
                                 facility_name = NULL,
                                 pathogen_name = NULL) {
  # -- Input validation -------------------------------------------------------
  required_cols <- c(
    syndrome_col, polymicrobial_col,
    patient_col, pathogen_col, outcome_col
  )
  if (!is.null(specimen_col)) required_cols <- c(required_cols, specimen_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }
  if (!is.null(specimen_col) && is.null(specimen_name)) {
    stop("specimen_name must be provided when specimen_col is specified.")
  }

  # -- Step 1: Filter syndrome + (optional specimen) + fatal outcome ---------
  df <- data %>%
    dplyr::filter(
      .data[[syndrome_col]] == syndrome_name,
      .data[[outcome_col]] == death_value
    )
  if (!is.null(specimen_col)) {
    df <- df %>% dplyr::filter(.data[[specimen_col]] == specimen_name)
  }
  if (nrow(df) == 0) {
    stop(sprintf(
      "No fatal records remain after syndrome%s filtering (death_value='%s').",
      if (!is.null(specimen_col)) "/specimen" else "",
      death_value
    ))
  }

  # -- Step 2: Optional single-facility restriction ---------------------------
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No fatal records found for facility '%s'.", facility_name))
    }
  }

  # -- Step 3: Optional pathogen filter --------------------------------------
  if (!is.null(pathogen_name)) {
    df <- df %>% dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(df) == 0) {
      stop(sprintf(
        "No fatal records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  # -- Step 4: GLASS filter -- polymicrobial patients only --------------------
  # Monomicrobial patients (polymicrobial_col == 0) are never filtered.
  if (!is.null(glass_ref)) {
    if (is.data.frame(glass_ref)) {
      if (!is.null(specimen_name)) {
        valid_pathogens <- glass_ref %>%
          dplyr::filter(.data[["specimen"]] == specimen_name) %>%
          dplyr::pull(.data[["pathogen"]])
      } else {
        # No specimen filter: use all pathogens in the reference
        valid_pathogens <- unique(glass_ref[["pathogen"]])
      }
    } else {
      valid_pathogens <- glass_ref
    }

    df_mono <- df %>% dplyr::filter(.data[[polymicrobial_col]] == 0)
    df_poly <- df %>%
      dplyr::filter(.data[[polymicrobial_col]] == 1) %>%
      dplyr::filter(.data[[pathogen_col]] %in% valid_pathogens)

    n_removed <- sum(df[[polymicrobial_col]] == 1) - nrow(df_poly)
    if (n_removed > 0) {
      message(sprintf(
        "GLASS filter: removed %d row(s) from polymicrobial fatal patients (pathogen not on GLASS list%s).",
        n_removed,
        if (!is.null(specimen_name)) sprintf(" for '%s'", specimen_name) else ""
      ))
    }
    df <- dplyr::bind_rows(df_mono, df_poly)
    if (nrow(df) == 0) {
      stop("No fatal records remain after GLASS reference filtering.")
    }
  }

  # -- Step 5: Patient-level fractional weight 1/m_r -------------------------
  # m_r = distinct valid pathogens for fatal patient r.
  group_cols <- if (!is.null(facility_col)) c(facility_col, patient_col) else patient_col

  df <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      m_r    = dplyr::n_distinct(.data[[pathogen_col]]),
      weight = 1 / m_r
    ) %>%
    dplyr::ungroup()

  # -- Step 6: N^F_LK = weighted sum per (facility, pathogen) ---------------
  agg_cols <- if (!is.null(facility_col)) c(facility_col, pathogen_col) else pathogen_col

  N_LK <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) %>%
    dplyr::summarise(N_F_LK = sum(weight, na.rm = TRUE), .groups = "drop")

  # -- Step 7: N^F_L = unique fatal patients per facility --------------------
  fac_grp <- if (!is.null(facility_col)) facility_col else character(0)

  N_L <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(fac_grp))) %>%
    dplyr::summarise(
      N_F_L = dplyr::n_distinct(.data[[patient_col]]),
      .groups = "drop"
    )

  # -- Step 8: Compute P_LK and return ---------------------------------------
  if (!is.null(facility_col)) {
    facility_level <- dplyr::left_join(N_LK, N_L, by = facility_col) %>%
      dplyr::mutate(P_Lk_fatal = N_F_LK / N_F_L)

    # Pooled: sum numerators / sum denominators (not average of ratios)
    pooled <- N_LK %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(N_F_LK = sum(N_F_LK), .groups = "drop") %>%
      dplyr::mutate(
        N_F_L      = sum(N_L$N_F_L),
        P_Lk_fatal = N_F_LK / N_F_L
      )

    message(sprintf(
      "P_LK (fatal) computed: %d pathogen(s) across %d facility/facilities; pooled P_LK also returned.",
      dplyr::n_distinct(df[[pathogen_col]]),
      dplyr::n_distinct(df[[facility_col]])
    ))
    return(list(P_Lk_fatal = pooled, facility_level = facility_level))
  } else {
    N_total <- N_L$N_F_L
    result <- N_LK %>%
      dplyr::mutate(
        N_F_L      = N_total,
        P_Lk_fatal = N_F_LK / N_F_L
      )
    message(sprintf(
      "P_LK (fatal) computed: %d pathogen(s), N^F_L = %d fatal patient(s).",
      nrow(result), N_total
    ))
    return(list(P_Lk_fatal = result))
  }
}


# -- CFR_LK : Case fatality ratio by syndrome and pathogen ----------------------

#' Calculate case fatality ratio by syndrome and pathogen (CFR_\{Lk\})
#'
#' Computes the case fatality ratio (CFR) for each pathogen (k) within a
#' specified infectious syndrome (L) using facility-level microbiology data.
#'
#' The unit of analysis is the \strong{patient}. Each patient contributes
#' total weight 1, distributed equally across their valid pathogens via a
#' fractional weight \eqn{1/m_r} (where \eqn{m_r} is the number of distinct
#' valid pathogens for patient r). This ensures polymicrobial patients are
#' not double-counted.
#'
#' For \strong{polymicrobial patients} (\code{polymicrobial_col == 1}), only
#' pathogens listed in the GLASS reference (\code{glass_ref}) for the given
#' specimen type are retained before weighting. Monomicrobial patients
#' (\code{polymicrobial_col == 0}) are never filtered.
#'
#' When \code{facility_col} is supplied and \code{facility_name} is NULL,
#' both per-facility and pooled CFR are returned. The pooled CFR sums
#' weighted deaths and totals across facilities before dividing.
#'
#' @param data             Data frame of facility-level microbiology records.
#' @param syndrome_col     Character. Column containing infectious syndrome labels (L).
#' @param syndrome_name    Character. Syndrome to analyse.
#' @param specimen_col     Character. Column containing specimen type.
#' @param specimen_name    Character. Specimen to restrict to (e.g., \code{"Blood"}).
#' @param polymicrobial_col Character. Column flagging polymicrobial patients
#'   (1 = polymicrobial, 0 = monomicrobial).
#' @param patient_col      Character. Unique patient identifier column.
#' @param pathogen_col     Character. Pathogen (organism) column (k).
#' @param outcome_col      Character. Final patient outcome column.
#' @param death_value      Character. Value indicating death. Default \code{"Died"}.
#' @param glass_ref        Character vector of valid pathogen names, or a data
#'   frame with columns \code{specimen} and \code{pathogen}. Applied to
#'   polymicrobial patients only. \code{NULL} skips GLASS filtering.
#' @param facility_col     Character or NULL. Facility identifier column.
#' @param facility_name    Character or NULL. Restricts to a single facility.
#' @param pathogen_name    Character vector or NULL. Filter to specific pathogen(s).
#'
#' @return A list:
#'   \describe{
#'     \item{cfr_table}{Pooled CFR: \code{pathogen_col}, \code{weighted_deaths},
#'       \code{weighted_total}, \code{CFR_LK}.}
#'     \item{facility_level}{Per-facility CFR (only when \code{facility_col} is
#'       supplied and \code{facility_name} is NULL).}
#'   }
#'
#' @export

daly_calc_case_fatality <- function(data,
                             syndrome_col,
                             syndrome_name,
                             specimen_col,
                             specimen_name,
                             polymicrobial_col,
                             patient_col,
                             pathogen_col,
                             outcome_col,
                             death_value = "Died",
                             glass_ref = NULL,
                             facility_col = NULL,
                             facility_name = NULL,
                             pathogen_name = NULL) {
  # -- Input validation ------------------------------------------------------
  required_cols <- c(
    syndrome_col, specimen_col, polymicrobial_col,
    patient_col, pathogen_col, outcome_col
  )
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  # -- Step 1: Filter syndrome + specimen -----------------------------------
  df <- data %>%
    dplyr::filter(
      .data[[syndrome_col]] == syndrome_name,
      .data[[specimen_col]] == specimen_name
    )
  if (nrow(df) == 0) {
    stop("No records remain after syndrome/specimen filtering.")
  }

  # -- Step 2: Optional single-facility restriction --------------------------
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No records found for facility '%s'.", facility_name))
    }
  }

  # -- Step 3: Optional pathogen filter -------------------------------------
  if (!is.null(pathogen_name)) {
    df <- df %>% dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(df) == 0) {
      stop(sprintf(
        "No records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  # -- Step 4: GLASS filter -- polymicrobial patients only -------------------
  if (!is.null(glass_ref)) {
    if (is.data.frame(glass_ref)) {
      valid_pathogens <- glass_ref %>%
        dplyr::filter(.data[["specimen"]] == specimen_name) %>%
        dplyr::pull(.data[["pathogen"]])
    } else {
      valid_pathogens <- glass_ref
    }

    df_mono <- df %>% dplyr::filter(.data[[polymicrobial_col]] == 0)
    df_poly <- df %>%
      dplyr::filter(.data[[polymicrobial_col]] == 1) %>%
      dplyr::filter(.data[[pathogen_col]] %in% valid_pathogens)

    n_removed <- sum(df[[polymicrobial_col]] == 1) - nrow(df_poly)
    if (n_removed > 0) {
      message(sprintf(
        "GLASS filter: removed %d row(s) from polymicrobial patients (pathogen not on GLASS list for '%s').",
        n_removed, specimen_name
      ))
    }
    df <- dplyr::bind_rows(df_mono, df_poly)
    if (nrow(df) == 0) {
      stop("No records remain after GLASS reference filtering.")
    }
  }

  # -- Step 5: Patient-level fractional weight 1/m_r ------------------------
  group_cols <- if (!is.null(facility_col)) c(facility_col, patient_col) else patient_col

  df <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      m_r = dplyr::n_distinct(.data[[pathogen_col]]),
      weight = 1 / m_r,
      fatal = .data[[outcome_col]] == death_value
    ) %>%
    dplyr::ungroup()

  # -- Step 6: Weighted deaths + totals per (facility, pathogen) -------------
  agg_cols <- if (!is.null(facility_col)) c(facility_col, pathogen_col) else pathogen_col

  cfr_raw <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) %>%
    dplyr::summarise(
      weighted_deaths = sum(weight * fatal, na.rm = TRUE),
      weighted_total = sum(weight, na.rm = TRUE),
      CFR_LK = dplyr::if_else(
        weighted_total > 0,
        weighted_deaths / weighted_total,
        NA_real_
      ),
      .groups = "drop"
    )

  # -- Step 7: Pooled across facilities (if applicable) ----------------------
  if (!is.null(facility_col) && is.null(facility_name)) {
    pooled <- cfr_raw %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(
        weighted_deaths = sum(weighted_deaths),
        weighted_total = sum(weighted_total),
        CFR_LK = dplyr::if_else(
          weighted_total > 0,
          weighted_deaths / weighted_total,
          NA_real_
        ),
        .groups = "drop"
      )

    message(sprintf(
      "CFR_LK computed: %d pathogen(s) across %d facility/facilities; pooled CFR also returned.",
      dplyr::n_distinct(df[[pathogen_col]]),
      dplyr::n_distinct(df[[facility_col]])
    ))
    return(list(cfr_table = pooled, facility_level = cfr_raw))
  }

  message(sprintf(
    "CFR_LK computed: %d pathogen(s).",
    dplyr::n_distinct(df[[pathogen_col]])
  ))
  return(list(cfr_table = cfr_raw))
}


# -- Top N pathogens -----------------------------------------------------------

#' Identify top N pathogens by occurrence
#'
#' Ranks pathogens by the number of records (rows) in the dataset, optionally
#' filtered to a specific syndrome, specimen, facility, or outcome. Returns
#' the top N pathogens overall or broken down per facility.
#'
#' Use this to decide which pathogens to focus on before calling
#' \code{calculate_P_Lk_prime_BSI()}, \code{calculate_cfr_lk()}, or
#' \code{calculate_YLD()}.
#'
#' @param data          Data frame of facility-level records.
#' @param pathogen_col  Character. Pathogen column.
#' @param n             Integer. Number of top pathogens to return. Default 5.
#' @param syndrome_col  Character or NULL. Filter to this syndrome column.
#' @param syndrome_name Character or NULL. Syndrome value to filter to.
#' @param specimen_col  Character or NULL. Filter to this specimen column.
#' @param specimen_name Character or NULL. Specimen value to filter to.
#' @param outcome_col   Character or NULL. Filter to this outcome column.
#' @param outcome_name  Character or NULL. Outcome value to filter to
#'   (e.g., \code{"Died"} or \code{"Discharged"}).
#' @param facility_col  Character or NULL. Facility identifier column.
#'   When provided without \code{facility_name}, returns top N per facility.
#' @param facility_name Character or NULL. If provided, restricts to that
#'   facility only before ranking.
#'
#' @return Data frame with columns: \code{pathogen_col}, \code{n_records},
#'   \code{rank} (1 = most common), and \code{facility_col} if supplied.
#' @export

daly_get_top_pathogens <- function(data,
                              pathogen_col,
                              n = 5L,
                              syndrome_col = NULL,
                              syndrome_name = NULL,
                              specimen_col = NULL,
                              specimen_name = NULL,
                              outcome_col = NULL,
                              outcome_name = NULL,
                              facility_col = NULL,
                              facility_name = NULL) {
  # -- Input validation ------------------------------------------------------
  if (!pathogen_col %in% names(data)) {
    stop(sprintf("pathogen_col '%s' not found in data.", pathogen_col))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  df <- data

  # -- Optional filters ------------------------------------------------------
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    df <- df %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
  }
  if (!is.null(specimen_col) && !is.null(specimen_name)) {
    df <- df %>% dplyr::filter(.data[[specimen_col]] == specimen_name)
  }
  if (!is.null(outcome_col) && !is.null(outcome_name)) {
    df <- df %>% dplyr::filter(.data[[outcome_col]] == outcome_name)
  }
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
  }

  if (nrow(df) == 0) {
    warning("No records remain after applying filters.")
    return(data.frame())
  }

  # -- Count and rank per facility (if requested) or overall -----------------
  group_vars <- if (!is.null(facility_col)) {
    c(facility_col, pathogen_col)
  } else {
    pathogen_col
  }

  ranked <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(n_records = dplyr::n(), .groups = "drop")

  if (!is.null(facility_col)) {
    ranked <- ranked %>%
      dplyr::group_by(.data[[facility_col]]) %>%
      dplyr::mutate(rank = dplyr::min_rank(dplyr::desc(n_records))) %>%
      dplyr::filter(rank <= n) %>%
      dplyr::arrange(.data[[facility_col]], rank) %>%
      dplyr::ungroup()
  } else {
    ranked <- ranked %>%
      dplyr::mutate(rank = dplyr::min_rank(dplyr::desc(n_records))) %>%
      dplyr::filter(rank <= n) %>%
      dplyr::arrange(rank)
  }

  message(sprintf(
    "Top %d pathogen(s)%s returned.",
    n,
    if (!is.null(facility_col) && is.null(facility_name)) {
      sprintf(" per facility (%d facilities)", dplyr::n_distinct(df[[facility_col]]))
    } else {
      ""
    }
  ))

  return(ranked)
}
# -- YLD per pathogen ----------------------------------------------------------

#' Calculate YLD per pathogen
#'
#' Computes Years Lived with Disability (YLD) attributable to each pathogen K
#' within syndrome L:
#'
#' \deqn{YLD_K = \text{Incidence}_L \times P'_{LK} \times \text{DW\_sepsis}}
#'
#' The YLD weight per incident case is drawn from a GBD-derived reference table
#' stratified by Indian state (\code{yld_ref}) when \code{DW_sepsis} is not
#' provided. If \code{DW_sepsis} is provided directly, that scalar is used for
#' all rows and \code{yld_ref} is ignored.
#'
#' \strong{Facility-level mode} (when \code{facility_col} is supplied and
#' \code{facility_name} is NULL): \code{P_Lk_prime_tbl} must contain a
#' \code{facility_col} column (i.e., the \code{facility_level} element from
#' \code{calculate_P_Lk_prime()}). \code{incidence_data} must be a data frame
#' with \code{facility_col} and an incidence count column. A
#' \code{facility_state_map} data frame linking each facility to its state is
#' required only when \code{DW_sepsis} is not supplied.
#'
#' \strong{Pooled / no-facility mode}: \code{incidence_data} is a single
#' numeric scalar and \code{P_Lk_prime_tbl} has no facility column.
#' \code{state_name} selects the YLD weight; defaults to "India" if NULL.
#'
#' @param incidence_data  Numeric scalar (pooled) or data frame with
#'   \code{facility_col} + \code{incidence_col} columns (facility-level).
#' @param P_Lk_prime_tbl  Data frame from \code{calculate_P_Lk_prime()}.
#'   Use the \code{P_Lk_prime} element for pooled mode or the
#'   \code{facility_level} element for facility-level mode.
#' @param yld_ref         Data frame with columns \code{location_name} and
#'   \code{DW_sepsis}. Loaded from
#'   \file{inst/extdata/Proxy_YLD_per_case.xlsx}. Optional if
#'   \code{DW_sepsis} is provided.
#' @param DW_sepsis Numeric scalar or NULL. If provided, this value is used
#'   directly for all rows and \code{yld_ref} is ignored.
#' @param avg_los_years Numeric scalar or NULL. Overall average length of stay
#'   in years across all patients (resistant and susceptible combined), used to
#'   convert DW to a duration-weighted value:
#'   \code{effective_DW = DW_sepsis * avg_los_years}.
#'   If NULL, DW_sepsis is used as-is (caller is responsible for duration weighting).
#' @param state_name      Character or NULL. State for YLD weight in pooled
#'   mode. NULL uses the India row. Ignored in facility-level mode when
#'   \code{DW_sepsis} is provided.
#' @param facility_col    Character or NULL. Facility identifier column.
#' @param facility_name   Character or NULL. Restrict to one facility only.
#' @param facility_state_map Data frame with \code{facility_col} and
#'   \code{state_col} columns, mapping each facility to its state.
#'   Required in facility-level mode only when \code{DW_sepsis} is not provided.
#' @param state_col       Character. Column in \code{facility_state_map}
#'   containing state names that match \code{yld_ref$location_name}.
#'   Default \code{"state"}.
#' @param pathogen_col    Character. Pathogen column in \code{P_Lk_prime_tbl}.
#'   Default \code{"pathogen"}.
#' @param pathogen_name   Character vector or NULL. If provided, restricts
#'   output to those pathogen(s) only.
#' @param plk_col         Character. P'LK column in \code{P_Lk_prime_tbl}.
#'   Default \code{"P_Lk_prime"}.
#' @param incidence_col   Character. Incidence column when
#'   \code{incidence_data} is a data frame. Default \code{"n_cases"}.
#'
#' @return Data frame with columns: \code{pathogen_col}, \code{P_Lk_prime},
#'   \code{incidence_L}, \code{DW_sepsis}, \code{YLD},
#'   and \code{facility_col} / \code{state_col} when in facility-level mode.
#' @export
