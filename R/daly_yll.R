# daly_yll.R
# YLL (Years of Life Lost) calculation functions for AMR burden estimation

daly_load_life_expectancy <- function(le_path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to load the life expectancy table.")
  }
  if (!file.exists(le_path)) {
    stop(sprintf("Life expectancy file not found: %s", le_path))
  }

  raw <- readxl::read_excel(le_path, col_names = FALSE, .name_repair = "minimal")

  # Table layout (1-indexed rows) for life_expectancy_all.xlsx:
  #   Row  1       : title  ("LIFE EXPECTANCY 2019-2023 INDIA")
  #   Row  2       : header (Age category, India, state names ...)
  #   Rows  3 - 21 : Combined, 19 age bins
  #   Row  22      : "Male" section header
  #   Rows 23 - 41 : Male, 19 age bins
  #   Row  42      : "Female" section header
  #   Rows 43 - 61 : Female, 19 age bins
  # India life expectancy is always in column 2.

  age_bins <- c(
    "0-1", "1-5", "5-10", "10-15", "15-20", "20-25", "25-30",
    "30-35", "35-40", "40-45", "45-50", "50-55", "55-60",
    "60-65", "65-70", "70-75", "75-80", "80-85", "85+"
  )

  india_col <- 2L
  parse_le <- function(rows) suppressWarnings(as.numeric(raw[[india_col]][rows]))

  rbind(
    data.frame(
      sex = "Combined", age_bin = age_bins,
      life_expectancy = parse_le(3:21), stringsAsFactors = FALSE
    ),
    data.frame(
      sex = "Male", age_bin = age_bins,
      life_expectancy = parse_le(23:41), stringsAsFactors = FALSE
    ),
    data.frame(
      sex = "Female", age_bin = age_bins,
      life_expectancy = parse_le(43:61), stringsAsFactors = FALSE
    )
  )
}


#' Compute YLL Associated with AMR (Patient-Level, Facility-Direct)
#'
#' Computes years of life lost (YLL) associated with AMR directly from
#' patient-level facility records, without requiring population-level
#' P_LK (syndrome fractions) or R_kd (resistance profile scalars).
#'
#' For every fatal patient with pathogen k (and optionally syndrome L) the
#' individual YLL contribution is:
#' \deqn{\text{YLL}_{r,k} = \text{LE}(\text{age\_bin}_r, \text{sex}_r)
#'   \times w_{r,k}}
#' where \eqn{w_{r,k}} is the polymicrobial death weight for patient r and
#' pathogen k (= 1 for monomicrobial, 0-1 for polymicrobial episodes).
#' Total YLL associated:
#' \deqn{\text{YLL}_{\text{associated}} = \sum_{r,k} \text{YLL}_{r,k}}
#'
#' \strong{Polymicrobial weights} are computed via
#' \code{flag_polymicrobial()} + \code{compute_polymicrobial_weight()}
#' from \code{weight.R}.  When \code{facility_col} is provided the weights
#' are derived per facility (reflecting local organism distributions), with
#' automatic fallback to globally-pooled proportions for facilities whose
#' monomicrobial reference pool is smaller than \code{min_mono_per_facility}.
#' If \code{date_culture_col} is \code{NULL} polymicrobial flagging is
#' skipped and all weights default to 1.
#'
#' @param data                  Data frame of patient-level facility records.
#' @param outcome_col           Character.  Final outcome column.
#' @param death_value           Character.  Value(s) indicating a fatal outcome.
#'   Default \code{"Death"}.  Pass a vector to match multiple labels
#'   (e.g. \code{c("Death","Died")}).
#' @param pathogen_col          Character.  Pathogen / organism column (k).
#' @param patient_col           Character.  Unique patient identifier column.
#' @param age_bin_col           Character.  Column containing GBD-standard age
#'   bin labels (e.g. \code{"0-1"}, \code{"1-5"}, ..., \code{"85+"}).  Use
#'   \code{age_bin_map} to recode non-standard labels.
#' @param sex_col               Character.  Column containing patient sex.
#' @param facility_col          Character or \code{NULL}.  Facility identifier
#'   column.  When supplied, polymicrobial weights are computed per facility
#'   and results include a \code{by_facility} breakdown.  Default \code{NULL}.
#' @param syndrome_col          Character or \code{NULL}.  Syndrome column.
#'   When \code{NULL} all syndromes are pooled; when supplied the
#'   \code{by_syndrome} and \code{by_syndrome_pathogen} outputs are populated.
#' @param syndrome_name         Character or \code{NULL}.  If supplied, data
#'   are filtered to this syndrome before computation.  \code{NULL} = all.
#' @param date_culture_col      Character or \code{NULL}.  Culture date column
#'   used for polymicrobial episode detection.  When \code{NULL} polymicrobial
#'   flagging is skipped and all weights default to 1.
#' @param specimen_col          Character or \code{NULL}.  Specimen type column
#'   used alongside \code{date_culture_col} for polymicrobial detection.
#' @param le_path               Character.  Path to the India life expectancy
#'   xlsx file.  Defaults to the bundled \code{inst/extdata} copy.
#' @param male_value            Character.  Value in \code{sex_col} for males.
#'   Default \code{"Male"}.
#' @param female_value          Character.  Value in \code{sex_col} for females.
#'   Default \code{"Female"}.  All other values use the combined LE.
#' @param age_bin_map           Named character vector remapping non-standard
#'   age bin labels to LE-table labels.  Default \code{c("<1" = "0-1")}.
#' @param poly_weight_method    Character.  Method for computing polymicrobial
#'   death weights.  Default \code{"monomicrobial_proportion"}.
#' @param min_mono_per_facility Integer.  Minimum monomicrobial records per
#'   facility to use facility-specific weights; smaller facilities fall back
#'   to global weights.  Default \code{30L}.
#' @param gap_days              Integer.  Window (days) used for polymicrobial
#'   episode detection.  Default \code{14}.
#' @param stratify_by           Character vector or \code{NULL}.  Additional
#'   columns to aggregate results by.  Default \code{NULL}.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{total}}{Scalar: total YLL associated across all pathogens and
#'     facilities.}
#'   \item{\code{per_pathogen}}{Data frame: YLL summed per pathogen k, pooled
#'     across facilities.  Columns: \code{pathogen_col}, \code{n_patients},
#'     \code{YLL_associated_k}.}
#'   \item{\code{by_age_sex}}{Data frame: YLL by \code{age_bin_col} x sex.}
#'   \item{\code{by_pathogen_age_sex}}{Data frame: YLL by pathogen x
#'     \code{age_bin_col} x sex.}
#'   \item{\code{by_facility}}{Data frame (only when \code{facility_col} is
#'     supplied): per-facility YLL, one row per facility x pathogen.}
#'   \item{\code{by_syndrome}}{Data frame (only when \code{syndrome_col} is
#'     supplied and \code{syndrome_name} is \code{NULL}): YLL by syndrome.}
#'   \item{\code{by_syndrome_pathogen}}{Data frame (only when
#'     \code{syndrome_col} supplied): YLL by syndrome x pathogen.}
#'   \item{\code{stratified}}{Data frame (only when \code{stratify_by} is
#'     supplied): YLL aggregated by the requested columns.}
#'   \item{\code{patient_data}}{The death-cohort data frame used for
#'     computation, with \code{polymicrobial_weight}, \code{life_expectancy},
#'     and \code{yll_contribution} columns attached.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' yll <- daly_calc_yll_associated(
#'   data             = cohort_df,
#'   outcome_col      = "final_outcome",
#'   death_value      = "Death",
#'   pathogen_col     = "organism_name",
#'   patient_col      = "PatientInformation_id",
#'   age_bin_col      = "Age_bin",
#'   sex_col          = "gender",
#'   facility_col     = "center_name",
#'   syndrome_col     = "infectious_syndrome",
#'   date_culture_col = "culture_date",
#'   specimen_col     = "sample_type",
#'   stratify_by      = c("location", "infectious_syndrome")
#' )
#' yll$total
#' yll$per_pathogen
#' yll$by_age_sex
#' yll$stratified
#' }
daly_calc_yll_associated <- function(
  data,
  outcome_col,
  death_value = "Death",
  pathogen_col,
  patient_col,
  age_bin_col,
  sex_col,
  facility_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  date_culture_col = NULL,
  specimen_col = NULL,
  poly_weight_method = "monomicrobial_proportion",
  min_mono_per_facility = 30L,
  gap_days = 14,
  le_path = system.file("extdata", "life_expectancy_all.xlsx", package = "anumaan"),
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c("<1" = "0-1"),
  stratify_by = NULL
) {
  # -- Input validation -------------------------------------------------------
  required_cols <- c(
    outcome_col, pathogen_col, patient_col,
    age_bin_col, sex_col
  )
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)
  if (!is.null(date_culture_col)) required_cols <- c(required_cols, date_culture_col)
  if (!is.null(specimen_col)) required_cols <- c(required_cols, specimen_col)
  if (!is.null(stratify_by)) required_cols <- c(required_cols, stratify_by)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col)) {
    stop("syndrome_col must be supplied when syndrome_name is specified.")
  }

  # -- Step 1: Load life expectancy lookup -----------------------------------
  le_lookup <- daly_load_life_expectancy(le_path)

  # -- Step 2: Polymicrobial weights (per facility; global fallback) ---------
  # Computed on the FULL data before death filter so the monomicrobial
  # reference pool is not restricted to fatal cases only.


  .compute_poly_weights <- function(df_sub, label) {
    if (is.null(date_culture_col)) {
      df_sub$polymicrobial_weight <- 1.0
      df_sub$weight_method <- "none"
      df_sub$weight_confidence <- "n/a"
      return(df_sub)
    }

    df_flagged <- tryCatch(
      flag_polymicrobial(
        df_sub,
        patient_col  = patient_col,
        organism_col = pathogen_col
      ),
      error = function(e) {
        message(sprintf(
          "  [%s] flag_polymicrobial error: %s -- weights set to 1.",
          label, conditionMessage(e)
        ))
        df_sub$polymicrobial_weight <- 1
        df_sub$is_polymicrobial <- 0L
        df_sub$episode_id <- seq_len(nrow(df_sub))
        return(df_sub)
      }
    )

    # Ensure required columns exist
    if (!"is_polymicrobial" %in% names(df_flagged)) {
      message(sprintf("[%s] 'is_polymicrobial' missing -- assuming monomicrobial.", label))
      df_flagged$is_polymicrobial <- 0L
    }

    if (!"episode_id" %in% names(df_flagged)) {
      df_flagged$episode_id <- seq_len(nrow(df_flagged))
    }

    n_mono <- sum(df_flagged$is_polymicrobial == 0L, na.rm = TRUE)

    use_method <- if (n_mono >= min_mono_per_facility) {
      poly_weight_method
    } else {
      "equal"
    }

    if (n_mono < min_mono_per_facility) {
      message(sprintf(
        "  [%s] n_mono=%d < threshold %d -- falling back to equal weights.",
        label, n_mono, min_mono_per_facility
      ))
    }

    df_weighted <- compute_polymicrobial_weight(
      df_flagged,
      episode_col       = "episode_id",
      organism_col      = pathogen_col,
      polymicrobial_col = "is_polymicrobial",
      method            = use_method
    )

    return(df_weighted)
  }

  if (!is.null(facility_col)) {
    facilities <- sort(unique(data[[facility_col]]))
    message(sprintf(
      "Computing polymicrobial weights for %d facility/facilities...",
      length(facilities)
    ))
    data_weighted <- dplyr::bind_rows(lapply(facilities, function(fac) {
      df_fac <- data[data[[facility_col]] == fac, , drop = FALSE]
      .compute_poly_weights(df_fac, label = fac)
    }))
  } else {
    message("Computing polymicrobial weights (global)...")
    data_weighted <- .compute_poly_weights(data, label = "global")
  }

  # -- Step 3: Filter to death cohort ----------------------------------------
  df <- data_weighted %>%
    dplyr::filter(
      .data[[outcome_col]] %in% death_value,
      !is.na(.data[[pathogen_col]])
    )

  if (!is.null(syndrome_name)) {
    df <- df %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
  }

  if (nrow(df) == 0L) {
    stop("No fatal records found after applying filters. Check outcome_col, death_value, and syndrome_name.")
  }

  message(sprintf(
    "Death cohort: %d rows | %d patients | %d pathogens%s",
    nrow(df),
    dplyr::n_distinct(df[[patient_col]]),
    dplyr::n_distinct(df[[pathogen_col]]),
    if (!is.null(facility_col)) {
      sprintf(" | %d facilities", dplyr::n_distinct(df[[facility_col]]))
    } else {
      ""
    }
  ))

  # -- Step 4: Recode age bins; normalise sex --------------------------------
  if (length(age_bin_map) > 0L) {
    df[[age_bin_col]] <- dplyr::recode(
      as.character(df[[age_bin_col]]),
      !!!age_bin_map
    )
  }

  df <- df %>%
    dplyr::mutate(
      .sex_norm = dplyr::case_when(
        .data[[sex_col]] == male_value ~ "Male",
        .data[[sex_col]] == female_value ~ "Female",
        TRUE ~ "Combined"
      )
    )

  # -- Step 5: Join life expectancy (age bin x sex) --------------------------
  join_by <- stats::setNames(c("age_bin", "sex"), c(age_bin_col, ".sex_norm"))
  df <- dplyr::left_join(df, le_lookup, by = join_by)

  n_missing_le <- sum(is.na(df$life_expectancy))
  if (n_missing_le > 0L) {
    warning(sprintf(
      "%d row(s) had no life expectancy match (unrecognised age_bin or sex); YLL set to NA.",
      n_missing_le
    ))
  }

  # -- Step 6: YLL per patient-pathogen row ----------------------------------
  # death_weight = polymicrobial_weight from weight.R
  #   mono patient  -> weight = 1.0  (full LE attributed to this pathogen)
  #   poly patient  -> weight = 0-1  (fractional LE per pathogen)
  df <- df %>%
    dplyr::mutate(
      death_weight     = dplyr::coalesce(polymicrobial_weight, 1.0),
      yll_contribution = life_expectancy * death_weight
    )

  # -- Step 7: Aggregate -----------------------------------------------------

  .agg <- function(df, grp_cols) {
    df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) %>%
      dplyr::summarise(
        n_patients       = dplyr::n_distinct(.data[[patient_col]]),
        YLL_associated_k = sum(yll_contribution, na.rm = TRUE),
        .groups          = "drop"
      )
  }

  # Per pathogen (+ facility if supplied)
  path_grp <- c(if (!is.null(facility_col)) facility_col, pathogen_col)
  by_path_fac <- .agg(df, path_grp)

  # Pool to per-pathogen (across facilities)
  per_pathogen <- if (!is.null(facility_col)) {
    by_path_fac %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(
        n_patients       = sum(n_patients, na.rm = TRUE),
        YLL_associated_k = sum(YLL_associated_k, na.rm = TRUE),
        .groups          = "drop"
      )
  } else {
    by_path_fac
  }

  # By age_bin x sex
  by_age_sex <- .agg(df, c(age_bin_col, ".sex_norm")) %>%
    dplyr::rename(sex = ".sex_norm")

  # By pathogen x age_bin x sex
  by_pathogen_age_sex <- .agg(df, c(pathogen_col, age_bin_col, ".sex_norm")) %>%
    dplyr::rename(sex = ".sex_norm")

  # By facility (if supplied)
  by_facility <- if (!is.null(facility_col)) {
    .agg(df, c(facility_col, pathogen_col))
  } else {
    NULL
  }

  # By syndrome (if syndrome_col supplied and not filtered to one syndrome)
  by_syndrome <- if (!is.null(syndrome_col) && is.null(syndrome_name)) {
    .agg(df, syndrome_col)
  } else {
    NULL
  }

  by_syndrome_pathogen <- if (!is.null(syndrome_col)) {
    .agg(df, c(syndrome_col, pathogen_col))
  } else {
    NULL
  }

  # User-defined stratification
  stratified <- if (!is.null(stratify_by)) {
    .agg(df, unique(c(stratify_by, pathogen_col)))
  } else {
    NULL
  }

  # -- Step 8: Total and summary message -------------------------------------
  YLL_total <- sum(per_pathogen$YLL_associated_k, na.rm = TRUE)

  message(sprintf(
    "YLL associated: total = %.2f years | %d pathogens | %d patients",
    YLL_total,
    nrow(per_pathogen),
    sum(per_pathogen$n_patients, na.rm = TRUE)
  ))

  # -- Return -----------------------------------------------------------------
  out <- list(
    total               = YLL_total,
    per_pathogen        = per_pathogen,
    by_age_sex          = by_age_sex,
    by_pathogen_age_sex = by_pathogen_age_sex,
    patient_data        = df
  )
  if (!is.null(facility_col)) out$by_facility <- by_facility
  if (!is.null(syndrome_col)) out$by_syndrome <- by_syndrome
  if (!is.null(syndrome_col)) out$by_syndrome_pathogen <- by_syndrome_pathogen
  if (!is.null(stratify_by)) out$stratified <- stratified

  out
}


# -- YLL Attributable ----------------------------------------------------------

#' Compute YLL Attributable to AMR
#'
#' Takes the per-patient YLL data produced by \code{daly_calc_yll_associated()}
#' (which already contains \code{life_expectancy}, \code{death_weight}, and
#' \code{yll_contribution}) and multiplies by the mortality PAF from
#' \code{daly_calc_paf_mortality()} to produce AMR-attributable YLL.
#'
#' \strong{Two PAF modes:}
#' \describe{
#'   \item{PAF_k scalar mode (default)}{When \code{resistance_profile_col} is
#'     \code{NULL}, the overall mortality PAF per pathogen k (\code{PAF_k_mort})
#'     from \code{paf_mort} is used as a scalar multiplier.}
#'   \item{Per-profile mode}{When \code{resistance_profile_col} is supplied,
#'     each patient row is matched to its resistance profile delta and the
#'     profile-specific \code{PAF_mortality} is applied.  Patients whose
#'     profile does not appear in \code{paf_mort} receive \code{NA} and a
#'     warning is issued.}
#' }
#'
#' \deqn{\text{YLL}^{\text{attr}}_{i,k} =
#'   \text{yll\_contribution}_{i,k} \times \text{PAF}_{k(,\delta)}}
#'
#' @param yll_patient_data       Data frame: the \code{patient_data} element
#'   from \code{daly_calc_yll_associated()}, with \code{life_expectancy},
#'   \code{death_weight}, and \code{yll_contribution} columns.
#' @param paf_mort              Named list returned by
#'   \code{daly_calc_paf_mortality()}.
#' @param pathogen_col          Character.  Pathogen column.
#' @param patient_col           Character.  Patient identifier column.
#' @param age_bin_col           Character.  Age bin column.
#' @param sex_col               Character.  Normalised sex column.
#'   Default \code{".sex_norm"}.
#' @param resistance_profile_col Character or \code{NULL}.  Column holding the
#'   resistance profile identifier per patient row.  When supplied, per-profile
#'   PAF is applied instead of a scalar PAF_k.  Default \code{NULL}.
#' @param profile_col           Character.  Profile identifier column in the
#'   \code{paf_mort} per-profile data frames.  Default \code{"profile"}.
#' @param facility_col          Character or \code{NULL}.  Facility identifier
#'   column.  When supplied results include a \code{by_facility} breakdown.
#'   Default \code{NULL}.
#' @param syndrome_col          Character or \code{NULL}.  Syndrome column.
#' @param stratify_by           Character vector or \code{NULL}.  Additional
#'   columns to aggregate results by.  Default \code{NULL}.
#'
#' @return A list:
#'   \describe{
#'     \item{total}{Scalar: total AMR-attributable YLL.}
#'     \item{per_pathogen}{One row per pathogen: \code{pathogen_col},
#'       \code{n_patients}, \code{YLL_associated_k}, \code{PAF_k_mort},
#'       \code{YLL_attributable_k}.}
#'     \item{by_age_sex}{YLL attributable stratified by age_bin x sex.}
#'     \item{by_pathogen_age_sex}{YLL attributable by pathogen x age_bin x sex.}
#'     \item{by_facility}{Per-facility breakdown (when \code{facility_col}
#'       supplied).}
#'     \item{by_syndrome}{Per-syndrome breakdown (when \code{syndrome_col}
#'       supplied).}
#'     \item{by_syndrome_pathogen}{Per-syndrome x pathogen (when
#'       \code{syndrome_col} supplied).}
#'     \item{stratified}{User-defined stratification (when \code{stratify_by}
#'       supplied).}
#'     \item{patient_data}{Patient-level data augmented with \code{PAF_kd}
#'       and \code{YLL_attributable_contribution}.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' yll_assoc <- daly_calc_yll_associated(data = cohort, ...)
#' mort_or <- fit_mortality_rr_logistic(data = rr_data, ...)
#' profiles_or <- assign_rr_to_profiles(profiles_out,
#'   rr_table = mort_or,
#'   rr_col = "OR_death"
#' )
#' paf_mort <- daly_calc_paf_mortality(profiles_or)
#'
#' yll_attr <- daly_calc_yll_attributable(
#'   yll_patient_data = yll_assoc$patient_data,
#'   paf_mort         = paf_mort,
#'   pathogen_col     = "organism_name",
#'   patient_col      = "PatientInformation_id",
#'   age_bin_col      = "Age_bin",
#'   sex_col          = ".sex_norm",
#'   facility_col     = "center_name",
#'   syndrome_col     = "infectious_syndrome",
#'   stratify_by      = c("location", "infectious_syndrome")
#' )
#' }
daly_calc_yll_attributable <- function(
  yll_patient_data,
  paf_mort,
  pathogen_col,
  patient_col,
  age_bin_col,
  sex_col = ".sex_norm",
  resistance_profile_col = NULL,
  profile_col = "profile",
  facility_col = NULL,
  syndrome_col = NULL,
  stratify_by = NULL
) {
  # -- Input validation -------------------------------------------------------
  required_cols <- c(
    "yll_contribution", "life_expectancy", "death_weight",
    pathogen_col, patient_col, age_bin_col, sex_col
  )
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)
  if (!is.null(stratify_by)) required_cols <- c(required_cols, stratify_by)
  if (!is.null(resistance_profile_col)) required_cols <- c(required_cols, resistance_profile_col)

  missing_cols <- setdiff(required_cols, names(yll_patient_data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Missing column(s) in yll_patient_data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!is.list(paf_mort) || length(paf_mort) == 0L) {
    stop("paf_mort must be the non-empty named list from daly_calc_paf_mortality().")
  }

  use_profile_mode <- !is.null(resistance_profile_col)
  message(sprintf(
    "PAF mode: %s",
    if (use_profile_mode) "per-profile (PAF_kd)" else "PAF_k scalar"
  ))

  df <- yll_patient_data

  # -- Step 1: Build PAF lookup and join -------------------------------------

  # Always build scalar PAF_k lookup (used for per_pathogen output table).
  paf_scalar_df <- do.call(rbind, lapply(names(paf_mort), function(k) {
    data.frame(
      .pathogen = k,
      PAF_k_mort = paf_mort[[k]]$PAF_k_mort,
      stringsAsFactors = FALSE
    )
  }))
  names(paf_scalar_df)[1L] <- pathogen_col

  if (!use_profile_mode) {
    # Warn on pathogen mismatches
    pats_data <- unique(df[[pathogen_col]])
    pats_paf <- names(paf_mort)
    missing_paf <- setdiff(pats_data, pats_paf)
    missing_data <- setdiff(pats_paf, pats_data)
    if (length(missing_paf) > 0L) {
      warning(sprintf(
        "Pathogen(s) in yll_patient_data have no PAF in paf_mort: %s. YLL_attributable set to NA.",
        paste(missing_paf, collapse = ", ")
      ))
    }
    if (length(missing_data) > 0L) {
      warning(sprintf(
        "Pathogen(s) in paf_mort have no records in yll_patient_data: %s.",
        paste(missing_data, collapse = ", ")
      ))
    }

    df <- dplyr::left_join(df, paf_scalar_df, by = pathogen_col)
    df <- df %>%
      dplyr::mutate(
        PAF_kd                        = PAF_k_mort,
        YLL_attributable_contribution = yll_contribution * PAF_kd
      )
  } else {
    # Per-profile: build (pathogen, profile_col) -> PAF_mortality lookup
    paf_profile_df <- do.call(rbind, lapply(names(paf_mort), function(k) {
      pp <- paf_mort[[k]]$per_profile
      if (is.null(pp) || !is.data.frame(pp) || !profile_col %in% names(pp)) {
        return(NULL)
      }
      df_k <- pp[, c(profile_col, "PAF_mortality"), drop = FALSE]
      df_k[[pathogen_col]] <- k
      df_k
    }))
    if (is.null(paf_profile_df) || nrow(paf_profile_df) == 0L) {
      stop(paste0(
        "No per-profile PAF data found in paf_mort. ",
        "Check paf_mort$per_profile and the profile_col argument."
      ))
    }

    # Join: left(pathogen_col, resistance_profile_col) -> right(pathogen_col, profile_col)
    join_keys <- stats::setNames(
      c(pathogen_col, profile_col),
      c(pathogen_col, resistance_profile_col)
    )
    df <- dplyr::left_join(df, paf_profile_df, by = join_keys)

    n_missing_paf <- sum(is.na(df$PAF_mortality))
    if (n_missing_paf > 0L) {
      warning(sprintf(
        "%d row(s) had no PAF_mortality match (unrecognised pathogen/profile combination); YLL_attributable set to NA.",
        n_missing_paf
      ))
    }

    df <- df %>%
      dplyr::mutate(
        PAF_kd                        = PAF_mortality,
        YLL_attributable_contribution = yll_contribution * PAF_kd
      )
  }

  # -- Step 2: Aggregate -----------------------------------------------------

  .agg_attr <- function(df, grp_cols) {
    df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) %>%
      dplyr::summarise(
        n_patients         = dplyr::n_distinct(.data[[patient_col]]),
        YLL_associated_k   = sum(yll_contribution, na.rm = TRUE),
        YLL_attributable_k = sum(YLL_attributable_contribution, na.rm = TRUE),
        .groups            = "drop"
      )
  }

  # Rename sex_col -> "sex" in output tables
  .rename_sex <- function(tbl) {
    if (sex_col != "sex" && sex_col %in% names(tbl)) {
      names(tbl)[names(tbl) == sex_col] <- "sex"
    }
    tbl
  }

  # Per pathogen (+ facility if supplied)
  path_grp <- c(if (!is.null(facility_col)) facility_col, pathogen_col)
  by_path_fac <- .agg_attr(df, path_grp)

  # Pool to per-pathogen across facilities; attach PAF_k_mort
  per_pathogen <- if (!is.null(facility_col)) {
    by_path_fac %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(
        n_patients         = sum(n_patients, na.rm = TRUE),
        YLL_associated_k   = sum(YLL_associated_k, na.rm = TRUE),
        YLL_attributable_k = sum(YLL_attributable_k, na.rm = TRUE),
        .groups            = "drop"
      )
  } else {
    by_path_fac
  }
  per_pathogen <- dplyr::left_join(per_pathogen, paf_scalar_df, by = pathogen_col)

  # By age_bin x sex
  by_age_sex <- .agg_attr(df, c(age_bin_col, sex_col)) %>%
    .rename_sex()

  # By pathogen x age_bin x sex
  by_pathogen_age_sex <- .agg_attr(df, c(pathogen_col, age_bin_col, sex_col)) %>%
    .rename_sex()

  # By facility
  by_facility <- if (!is.null(facility_col)) {
    .agg_attr(df, c(facility_col, pathogen_col))
  } else {
    NULL
  }

  # By syndrome
  by_syndrome <- if (!is.null(syndrome_col)) {
    .agg_attr(df, syndrome_col)
  } else {
    NULL
  }

  by_syndrome_pathogen <- if (!is.null(syndrome_col)) {
    .agg_attr(df, c(syndrome_col, pathogen_col))
  } else {
    NULL
  }

  # User-defined stratification
  stratified <- if (!is.null(stratify_by)) {
    .agg_attr(df, unique(c(stratify_by, pathogen_col)))
  } else {
    NULL
  }

  # -- Step 3: Total and summary ----------------------------------------------
  YLL_attributable_total <- sum(per_pathogen$YLL_attributable_k, na.rm = TRUE)
  YLL_associated_total <- sum(per_pathogen$YLL_associated_k, na.rm = TRUE)

  message(sprintf(
    "YLL attributable: total = %.2f years | associated = %.2f years | %d pathogens | %d patients",
    YLL_attributable_total,
    YLL_associated_total,
    nrow(per_pathogen),
    sum(per_pathogen$n_patients, na.rm = TRUE)
  ))

  # -- Return -----------------------------------------------------------------
  out <- list(
    total               = YLL_attributable_total,
    per_pathogen        = per_pathogen,
    by_age_sex          = by_age_sex,
    by_pathogen_age_sex = by_pathogen_age_sex,
    patient_data        = df
  )
  if (!is.null(facility_col)) out$by_facility <- by_facility
  if (!is.null(syndrome_col)) out$by_syndrome <- by_syndrome
  if (!is.null(syndrome_col)) out$by_syndrome_pathogen <- by_syndrome_pathogen
  if (!is.null(stratify_by)) out$stratified <- stratified

  out
}


# -- PAF Mortality -------------------------------------------------------------

#' Compute Mortality Population Attributable Fraction per Resistance Profile
#'
#' Computes PAF_kd_mortality for each pathogen k and resistance profile delta
#' using the GBD multi-exposure Levin formula, substituting the mortality
#' odds ratio (OR_death) from \code{fit_mortality_rr_logistic()} in place of
#' the LOS relative risk used by \code{compute_paf_los()}:
#'
#' \deqn{\text{PAF}_{kd,\text{mort}} =
#'   \frac{R'_{K\delta}\,(\text{OR}_{K\delta} - 1)}
#'        {1 + \sum_\delta R'_{K\delta}\,(\text{OR}_{K\delta} - 1)}}
#'
#' The all-susceptible profile carries OR = 1 and contributes 0.
#' The denominator equals
#' \eqn{E[\text{OR}_k] = \sum_\delta R'_{K\delta} \cdot \text{OR}_{K\delta}},
#' numerically identical to the denominator produced by \code{compute_paf_los()}.
#'
#' Overall mortality PAF for pathogen k:
#'
#' \deqn{\text{PAF}_{k,\text{mort}} = \sum_\delta \text{PAF}_{kd,\text{mort}}
#'   = \frac{\sum_\delta R'_{K\delta}(\text{OR}_{K\delta}-1)}
#'          {1 + \sum_\delta R'_{K\delta}(\text{OR}_{K\delta}-1)}}
#'
#' \strong{Usage pipeline:}
#' \preformatted{
#'   # 1. Fit mortality OR per class
#'   mort_or <- fit_mortality_rr_logistic(data, ...)
#'
#'   # 2. Assign OR to profiles via max rule (rr_col = "OR_death")
#'   profiles_with_or <- assign_rr_to_profiles(
#'       profiles_output,
#'       rr_table = mort_or,
#'       rr_col   = "OR_death"
#'   )
#'
#'   # 3. Compute per-profile and overall mortality PAF
#'   paf_mort <- daly_calc_paf_mortality(profiles_with_or)
#' }
#'
#' @param profiles_with_rr Named list from \code{assign_rr_to_profiles()}
#'   called with \code{rr_col = "OR_death"}.  Each element is a profile
#'   data frame for one pathogen.
#' @param probability_col  Character.  Profile probability column name.
#'   Default \code{"probability"}.
#' @param rr_profile_col   Character.  Profile-level OR column as produced by
#'   \code{assign_rr_to_profiles()}.  Default \code{"RR_LOS_profile"}.
#' @param profile_col      Character.  Profile label column.
#'   Default \code{"profile"}.
#'
#' @return Named list (one entry per pathogen) containing:
#'   \itemize{
#'     \item \code{per_profile}: profile data frame augmented with
#'       \code{numerator_mort} (= \eqn{R'_{K\delta}(\text{OR}_{K\delta}-1)}),
#'       \code{PAF_mortality} (= numerator / denominator), and
#'       \code{denominator_mort} (= \eqn{1 + \sum_\delta} numerator).
#'     \item \code{PAF_k_mort}: overall mortality PAF for pathogen k.
#'     \item \code{denominator_mort}: shared denominator \eqn{E[\text{OR}_k]}.
#'   }
#' @export
daly_calc_paf_mortality <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile"
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col, profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found in profiles for '%s'.",
          col, path
        ))
      }
    }

    p <- df[[probability_col]] # R'_kd  (sums to 1)
    or <- df[[rr_profile_col]] # OR_death_kd
    numerator_vec <- p * (or - 1.0) # R'_kd * (OR_kd - 1)
    denom <- 1.0 + sum(numerator_vec, na.rm = TRUE)

    if (!is.finite(denom) || denom <= 0) {
      warning(sprintf(
        "'%s': PAF_mortality denominator = %.6g (must be > 0) -- all mortality ORs may be <= 1 or NA; skipping.",
        path, denom
      ))
      next
    }

    paf_vec <- numerator_vec / denom

    df$numerator_mort <- round(numerator_vec, 6L)
    df$PAF_mortality <- round(paf_vec, 6L)
    df$denominator_mort <- round(denom, 6L)

    paf_k_mort <- sum(paf_vec, na.rm = TRUE)

    out[[path]] <- list(
      per_profile      = df,
      PAF_k_mort       = round(paf_k_mort, 6L),
      denominator_mort = round(denom, 6L)
    )

    message(sprintf(
      "'%s': PAF_k_mortality = %.4f | E[OR_death] (denominator) = %.4f | %d profiles.",
      path, paf_k_mort, denom, nrow(df)
    ))
  }

  return(out)
}
