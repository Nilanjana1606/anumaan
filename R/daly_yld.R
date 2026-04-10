# daly_yld.R
# YLD (Years Lived with Disability) calculation functions for AMR burden estimation

daly_calc_yld_baseline <- function(incidence_data,
                          P_Lk_prime_tbl,
                          yld_ref = NULL,
                          DW_sepsis = NULL,
                          avg_los_years = NULL,
                          state_name = NULL,
                          facility_col = NULL,
                          facility_name = NULL,
                          facility_state_map = NULL,
                          state_col = "state",
                          pathogen_col = "pathogen",
                          pathogen_name = NULL,
                          plk_col = "P_Lk_prime",
                          incidence_col = "n_cases") {
  # -- Input validation ------------------------------------------------------
  use_scalar_proxy <- !is.null(DW_sepsis)

  if (use_scalar_proxy) {
    if (!is.numeric(DW_sepsis) || length(DW_sepsis) != 1 ||
      is.na(DW_sepsis)) {
      stop("DW_sepsis must be a single non-missing numeric value.")
    }
  } else {
    if (is.null(yld_ref) ||
      !all(c("location_name", "DW_sepsis") %in% names(yld_ref))) {
      stop(
        "Provide either DW_sepsis as a numeric scalar, or yld_ref ",
        "with columns: 'location_name', 'DW_sepsis'."
      )
    }
  }

  # Validate avg_los_years
  if (!is.null(avg_los_years)) {
    if (!is.numeric(avg_los_years) || length(avg_los_years) != 1 || is.na(avg_los_years) || avg_los_years <= 0) {
      stop("avg_los_years must be a single positive numeric value (mean LOS in years).")
    }
  }

  if (!pathogen_col %in% names(P_Lk_prime_tbl)) {
    stop(sprintf("pathogen_col '%s' not found in P_Lk_prime_tbl.", pathogen_col))
  }
  if (!plk_col %in% names(P_Lk_prime_tbl)) {
    stop(sprintf("plk_col '%s' not found in P_Lk_prime_tbl.", plk_col))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  pool_by_facility <- !is.null(facility_col) && is.null(facility_name) &&
    is.data.frame(incidence_data)

  # -- Pathogen filter -------------------------------------------------------
  if (!is.null(pathogen_name)) {
    P_Lk_prime_tbl <- P_Lk_prime_tbl %>%
      dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(P_Lk_prime_tbl) == 0) {
      stop(sprintf(
        "No rows in P_Lk_prime_tbl for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
  }

  # -- Single facility restriction -------------------------------------------
  if (!is.null(facility_name) && !is.null(facility_col)) {
    if (facility_col %in% names(P_Lk_prime_tbl)) {
      P_Lk_prime_tbl <- P_Lk_prime_tbl %>%
        dplyr::filter(.data[[facility_col]] == facility_name)
    }
    if (is.data.frame(incidence_data) &&
      facility_col %in% names(incidence_data)) {
      incidence_data <- incidence_data %>%
        dplyr::filter(.data[[facility_col]] == facility_name)
    }
  }

  # =========================================================================
  # POOLED / NO-FACILITY MODE
  # =========================================================================
  if (!pool_by_facility) {
    # Resolve scalar incidence
    if (is.data.frame(incidence_data)) {
      if (!incidence_col %in% names(incidence_data)) {
        stop(sprintf(
          "incidence_col '%s' not found in incidence_data.",
          incidence_col
        ))
      }
      incidence_L <- sum(incidence_data[[incidence_col]], na.rm = TRUE)
    } else {
      incidence_L <- as.numeric(incidence_data)
    }

    # Resolve YLD weight
    if (use_scalar_proxy) {
      loc <- "user_input"
      yld_weight <- as.numeric(DW_sepsis)
    } else {
      loc <- if (is.null(state_name)) "India" else state_name
      yld_weight <- yld_ref %>%
        dplyr::filter(location_name == loc) %>%
        dplyr::pull(DW_sepsis)

      if (length(yld_weight) == 0) {
        stop(sprintf(
          "Location '%s' not found in yld_ref. Available: %s",
          loc,
          paste(head(yld_ref$location_name, 10), collapse = ", ")
        ))
      }
    }

    effective_dw <- if (!is.null(avg_los_years)) yld_weight * avg_los_years else yld_weight

    result <- P_Lk_prime_tbl %>%
      dplyr::select(dplyr::all_of(c(pathogen_col, plk_col))) %>%
      dplyr::mutate(
        incidence_L = incidence_L,
        DW_sepsis = yld_weight,
        avg_los_years = if (!is.null(avg_los_years)) avg_los_years else NA_real_,
        effective_DW = effective_dw,
        YLD = incidence_L * .data[[plk_col]] * effective_dw
      )

    message(sprintf(
      "YLD computed (pooled, location='%s'): %d pathogen(s), effective_DW=%.6f, total YLD = %.2f.",
      loc, nrow(result), effective_dw, sum(result$YLD, na.rm = TRUE)
    ))

    return(result)
  }

  # =========================================================================
  # FACILITY-LEVEL MODE
  # =========================================================================

  if (!use_scalar_proxy) {
    if (is.null(facility_state_map)) {
      stop("facility_state_map is required in facility-level mode.")
    }
    if (!all(c(facility_col, state_col) %in% names(facility_state_map))) {
      stop(sprintf(
        "facility_state_map must have columns '%s' and '%s'.",
        facility_col, state_col
      ))
    }
  }

  if (!facility_col %in% names(P_Lk_prime_tbl)) {
    stop(sprintf(
      "facility_col '%s' not found in P_Lk_prime_tbl. ",
      facility_col,
      "Pass the 'facility_level' element from calculate_P_Lk_prime()."
    ))
  }
  if (!incidence_col %in% names(incidence_data)) {
    stop(sprintf(
      "incidence_col '%s' not found in incidence_data.",
      incidence_col
    ))
  }

  # Join P'LK with incidence
  result <- P_Lk_prime_tbl %>%
    dplyr::select(dplyr::all_of(c(facility_col, pathogen_col, plk_col))) %>%
    dplyr::left_join(
      incidence_data %>%
        dplyr::select(dplyr::all_of(c(facility_col, incidence_col))),
      by = facility_col
    )

  if (use_scalar_proxy) {
    result <- result %>%
      dplyr::mutate(DW_sepsis = as.numeric(.env$DW_sepsis))
  } else {
    result <- result %>%
      dplyr::left_join(
        facility_state_map %>%
          dplyr::select(dplyr::all_of(c(facility_col, state_col))) %>%
          dplyr::distinct(),
        by = facility_col
      ) %>%
      dplyr::left_join(
        yld_ref %>%
          dplyr::select(location_name, DW_sepsis) %>%
          dplyr::rename(!!state_col := location_name),
        by = state_col
      )

    # Warn if any facilities couldn't be matched to a state YLD weight
    missing_yld <- result %>%
      dplyr::filter(is.na(DW_sepsis)) %>%
      dplyr::pull(.data[[facility_col]]) %>%
      unique()

    if (length(missing_yld) > 0) {
      warning(sprintf(
        "No YLD weight found for facility/state: %s. Using India-wide fallback.",
        paste(missing_yld, collapse = ", ")
      ))
      india_yld <- yld_ref %>%
        dplyr::filter(location_name == "India") %>%
        dplyr::pull(DW_sepsis)

      result <- result %>%
        dplyr::mutate(
          DW_sepsis = dplyr::if_else(
            is.na(DW_sepsis), india_yld, DW_sepsis
          )
        )
    }
  }

  result <- result %>%
    dplyr::mutate(
      avg_los_years = if (!is.null(avg_los_years)) avg_los_years else NA_real_,
      effective_DW = if (!is.null(avg_los_years)) DW_sepsis * avg_los_years else DW_sepsis,
      YLD = .data[[incidence_col]] * .data[[plk_col]] * effective_DW
    )

  message(sprintf(
    "YLD computed (facility-level): %d facility/facilities, %d pathogen(s), total YLD = %.2f.",
    dplyr::n_distinct(result[[facility_col]]),
    dplyr::n_distinct(result[[pathogen_col]]),
    sum(result$YLD, na.rm = TRUE)
  ))

  return(result)
}


# -- PAF_LOS : population attributable fraction for length of stay ---------------

#' Compute PAF for length of stay per resistance profile
#'
#' Computes the population attributable fraction (PAF) for length of stay
#' from a named list of profile data frames, each containing resistance-profile
#' probabilities and profile-level relative risk of LOS.
#'
#' For each pathogen the formula is:
#'   PAF = sum_d R'_kd * (RR_kd - 1) / (1 + sum_d R'_kd * (RR_kd - 1))
#'
#' where R'_kd is the probability of resistance profile d and RR_kd is the
#' corresponding relative LOS multiplier.
#'
#' @param profiles_with_rr Named list returned by \code{assign_rr_to_profiles()}
#'   or \code{filter_profiles_to_rr_classes()}. Each entry is a data frame of
#'   resistance profiles for one pathogen.
#' @param probability_col Character. Column name for profile probability
#'   (must sum to 1 within each pathogen). Default \code{"probability"}.
#' @param rr_profile_col Character. Column name for profile-level LOS relative
#'   risk. Default \code{"RR_LOS_profile"}.
#' @param profile_col Character. Column name for the profile identifier.
#'   Default \code{"profile"}.
#'
#' @return Named list (one entry per pathogen) with columns \code{profile},
#'   \code{probability}, \code{rr_profile_col}, \code{numerator}, \code{PAF_LOS},
#'   and \code{denominator} added to the input profile data frame.
#' @export

daly_calc_paf_los <- function(
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
    rr <- df[[rr_profile_col]] # RR_LOS_kd
    numerator_vec <- p * (rr - 1.0) # R'_kd * (RR_kd - 1)
    denom <- 1.0 + sum(numerator_vec, na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) {
      warning(sprintf(
        "'%s': PAF denominator = %.6g (must be > 0) -- all relative LOS may be < 1 or NA; skipping.",
        path, denom
      ))
      next
    }
    paf_vec <- numerator_vec / denom

    df$numerator <- round(numerator_vec, 6L)
    df$PAF_LOS <- round(paf_vec, 6L)
    df$denominator <- round(denom, 6L)

    paf_k <- sum(paf_vec)

    out[[path]] <- list(
      per_profile = df,
      PAF_k       = round(paf_k, 6L),
      denominator = round(denom, 6L)
    )
    message(sprintf(
      "'%s': PAF_k = %.4f | E[relative LOS] (denominator) = %.4f | %d profiles.",
      path, paf_k, denom, nrow(df)
    ))
  }

  return(out)
}


# -- Step 5 --------------------------------------------------------------------

#' Compute Associated-Burden Fractions per Resistance Profile
#'
#' Computes the fraction of total expected YLD burden that *occurs in*
#' infections with each resistance profile delta (YLDs associated with
#' resistance).  This is a **burden partition**, not a counterfactual.
#'
#' For a single drug-class d the formula simplifies to:
#'
#'   Fraction_assoc_Kd = R'_Kd * RR_Kd / [(1 - R'_Kd) + R'_Kd * RR_Kd]
#'
#' With resistance profiles the denominator becomes the expected RR across
#' ALL profiles (including the all-susceptible profile with RR = 1):
#'
#'   E_RR_k = sum_delta  R'_K_delta * RR_K_delta
#'          = 1 + sum_delta  R'_K_delta * (RR_K_delta - 1)   [equivalent]
#'
#' Per-profile associated fraction:
#'   fraction_K_delta = R'_K_delta * RR_K_delta / E_RR_k
#'
#' Overall associated fraction (all resistant profiles combined):
#'   Fraction_k = sum_\{delta != 0\}  fraction_K_delta
#'
#' where delta != 0 denotes profiles with at least one resistant class.
#'
#' Note: E_RR_k is numerically identical to the `denominator` produced by
#' daly_calc_paf_los() -- both equal 1 + sum_d R'_kd*(RR_kd - 1).
#'
#' @param profiles_with_rr Named list from assign_rr_to_profiles() or
#'   filter_profiles_to_rr_classes().  Each entry is a profile data frame.
#' @param probability_col Character.  Profile probability column.
#'   Default \code{"probability"}.
#' @param rr_profile_col Character.  Profile-level RR column.
#'   Default \code{"RR_LOS_profile"}.
#'
#' @return Named list (one entry per pathogen) containing:
#'   \itemize{
#'     \item \code{per_profile}: profile data frame augmented with
#'       \code{numerator_assoc} (= p * rr) and \code{fraction_assoc}
#'       (= p * rr / E_RR_k).
#'     \item \code{Fraction_k}: overall associated fraction for the pathogen
#'       (sum of \code{fraction_assoc} over all resistant profiles).
#'     \item \code{E_RR_k}: expected RR = sum_delta R'_K_delta * RR_K_delta.
#'   }
#' @export
daly_calc_fraction_associated_yld <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  # Columns that are metadata / computed -- NOT class indicator columns
  non_class_cols <- c(
    "profile", probability_col, rr_profile_col,
    "dominant_class", "CI_lower_profile", "CI_upper_profile",
    "numerator", "PAF_LOS", "denominator",
    "numerator_assoc", "fraction_assoc"
  )

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found in profiles for '%s'.",
          col, path
        ))
      }
    }

    p <- df[[probability_col]] # R'_K_delta  (sums to 1 after normalisation)
    rr <- df[[rr_profile_col]] # RR_LOS_K_delta (1.0 for all-susceptible profile)

    # E[RR_k] = sum_delta  R'_K_delta * RR_K_delta
    # Equivalent: 1 + sum_delta R'_K_delta*(RR-1) = PAF denominator
    E_RR_k <- sum(p * rr)

    if (E_RR_k <= 0) {
      warning(sprintf("'%s': E[relative LOS] <= 0 -- skipping.", path))
      next
    }

    # Per-profile numerator and fraction
    numerator_assoc <- p * rr
    fraction_assoc_vec <- numerator_assoc / E_RR_k

    df$numerator_assoc <- round(numerator_assoc, 6L)
    df$fraction_assoc <- round(fraction_assoc_vec, 6L)

    # Resistant profiles: at least one binary class indicator column == 1.
    # The all-susceptible profile has every class column = 0 and RR = 1.
    class_cols <- setdiff(names(df), non_class_cols)
    if (length(class_cols) > 0L) {
      is_resistant <- rowSums(
        df[, class_cols, drop = FALSE] == 1L,
        na.rm = TRUE
      ) > 0L
    } else {
      # Fallback when class columns absent: use RR > 1 as proxy
      is_resistant <- rr > 1.0
    }

    Fraction_k <- sum(fraction_assoc_vec[is_resistant])

    out[[path]] <- list(
      per_profile = df,
      Fraction_k  = round(Fraction_k, 6L),
      E_RR_k      = round(E_RR_k, 6L)
    )

    message(sprintf(
      "'%s': Fraction_k (associated) = %.4f | E[relative LOS] = %.4f | %d profiles (%d resistant).",
      path, Fraction_k, E_RR_k, nrow(df), sum(is_resistant)
    ))
  }

  return(out)
}

daly_calc_yld_associated <- function(
  yld_k_tbl,
  fraction_assoc_list,
  pathogen_col = "pathogen",
  yld_col = "YLD",
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
) {
  if (!is.data.frame(yld_k_tbl)) {
    stop("yld_k_tbl must be a data frame (output of calculate_YLD()).")
  }
  if (!is.list(fraction_assoc_list)) {
    stop("fraction_assoc_list must be the list returned by daly_calc_fraction_associated_yld().")
  }
  if (!pathogen_col %in% names(yld_k_tbl)) {
    stop(sprintf("pathogen_col '%s' not found in yld_k_tbl.", pathogen_col))
  }
  if (!yld_col %in% names(yld_k_tbl)) {
    stop(sprintf("yld_col '%s' not found in yld_k_tbl.", yld_col))
  }

  # Build flat lookup: one row per pathogen
  frac_rows <- lapply(names(fraction_assoc_list), function(k) {
    res <- fraction_assoc_list[[k]]
    if (is.null(res)) {
      return(NULL)
    }

    # Defaults if per_profile unavailable
    R_k_delta_txt <- NA_character_
    LOS_k_delta_txt <- NA_character_

    if (!is.null(res$per_profile) && is.data.frame(res$per_profile)) {
      df <- res$per_profile

      # Classify resistant profiles the same way as Step 5
      non_class_cols <- c(
        "profile", probability_col, rr_profile_col,
        "dominant_class", "CI_lower_profile", "CI_upper_profile",
        "numerator", "PAF_LOS", "denominator",
        "numerator_assoc", "fraction_assoc"
      )

      class_cols <- setdiff(names(df), non_class_cols)
      rr <- if (rr_profile_col %in% names(df)) df[[rr_profile_col]] else rep(NA_real_, nrow(df))

      if (length(class_cols) > 0L) {
        is_resistant <- rowSums(df[, class_cols, drop = FALSE] == 1L, na.rm = TRUE) > 0L
      } else {
        is_resistant <- rr > 1.0
      }

      if (probability_col %in% names(df) && rr_profile_col %in% names(df)) {
        p_res <- df[[probability_col]][is_resistant]
        rr_res <- df[[rr_profile_col]][is_resistant]

        if (length(p_res) > 0L) {
          R_k_delta_txt <- paste(round(p_res, 6L), collapse = ",")
        }
        if (length(rr_res) > 0L) {
          LOS_k_delta_txt <- paste(round(rr_res, 6L), collapse = ",")
        }
      }
    }

    data.frame(
      .pathogen = k,
      Fraction_k = res$Fraction_k,
      R_k_delta = R_k_delta_txt,
      LOS_k_delta = LOS_k_delta_txt,
      stringsAsFactors = FALSE
    )
  })

  frac_df <- do.call(rbind, frac_rows)

  if (is.null(frac_df) || nrow(frac_df) == 0L) {
    stop("fraction_assoc_list is empty -- no fractions to join.")
  }

  names(frac_df)[names(frac_df) == ".pathogen"] <- pathogen_col

  out <- merge(yld_k_tbl, frac_df, by = pathogen_col, all.x = TRUE)

  out$YLD_associated <- out[[yld_col]] * out$Fraction_k

  n_unmatched <- sum(is.na(out$Fraction_k))
  if (n_unmatched > 0L) {
    warning(sprintf(
      "%d row(s) in yld_k_tbl had no matching Fraction_k; YLD_associated set to NA.",
      n_unmatched
    ))
  }

  message(sprintf(
    "YLD_associated computed: %d pathogen(s), total YLD_associated = %.4f.",
    sum(!is.na(out$YLD_associated)),
    sum(out$YLD_associated, na.rm = TRUE)
  ))

  return(out)
}


# -- Step 7 --------------------------------------------------------------------

#' Compute YLDs Attributable to Resistance
#'
#' Multiplies \code{YLD_k} (from \code{calculate_YLD()}) by the LOS-based
#' PAF (from \code{daly_calc_paf_los()}).
#'
#'   YLD_attributable_k = YLD_k * PAF_k
#'
#' Answers: "How much disability burden exists *only because* infections were
#' resistant instead of susceptible?"  This is a counterfactual -- it measures
#' the excess burden driven purely by resistance.
#'
#' Note: YLD_attributable_k < YLD_associated_k always, because
#'   PAF_k = Fraction_k * (1 - 1/E_RR_k)  <  Fraction_k.
#'
#' @param yld_k_tbl Data frame from \code{calculate_YLD()} containing at
#'   least a pathogen column and a YLD column.
#' @param paf_los_list Named list from \code{daly_calc_paf_los()}.
#' @param pathogen_col Character.  Pathogen column in \code{yld_k_tbl}.
#'   Default \code{"pathogen"}.
#' @param yld_col Character.  YLD column in \code{yld_k_tbl}.
#'   Default \code{"YLD"}.
#'
#' @return \code{yld_k_tbl} augmented with columns \code{PAF_k},
#'   \code{denominator}, and \code{YLD_attributable}.
#' @export
daly_calc_yld_attributable <- function(
  yld_k_tbl,
  paf_los_list,
  pathogen_col = "pathogen",
  yld_col = "YLD"
) {
  if (!is.data.frame(yld_k_tbl)) {
    stop("yld_k_tbl must be a data frame (output of calculate_YLD()).")
  }
  if (!is.list(paf_los_list)) {
    stop("paf_los_list must be the list returned by daly_calc_paf_los().")
  }
  if (!pathogen_col %in% names(yld_k_tbl)) {
    stop(sprintf("pathogen_col '%s' not found in yld_k_tbl.", pathogen_col))
  }
  if (!yld_col %in% names(yld_k_tbl)) {
    stop(sprintf("yld_col '%s' not found in yld_k_tbl.", yld_col))
  }

  # Build flat lookup: one row per pathogen
  paf_rows <- lapply(names(paf_los_list), function(k) {
    res <- paf_los_list[[k]]
    if (is.null(res)) {
      return(NULL)
    }
    data.frame(
      .pathogen = k,
      PAF_k = res$PAF_k,
      denominator = res$denominator,
      stringsAsFactors = FALSE
    )
  })
  paf_df <- do.call(rbind, paf_rows)

  if (is.null(paf_df) || nrow(paf_df) == 0L) {
    stop("paf_los_list is empty -- no PAFs to join.")
  }

  names(paf_df)[names(paf_df) == ".pathogen"] <- pathogen_col

  out <- merge(yld_k_tbl, paf_df, by = pathogen_col, all.x = TRUE)

  out$YLD_attributable <- out[[yld_col]] * out$PAF_k

  n_unmatched <- sum(is.na(out$PAF_k))
  if (n_unmatched > 0L) {
    warning(sprintf(
      "%d row(s) in yld_k_tbl had no matching PAF_k; YLD_attributable set to NA.",
      n_unmatched
    ))
  }

  message(sprintf(
    "YLD_attributable computed: %d pathogen(s), total YLD_attributable = %.4f.",
    sum(!is.na(out$YLD_attributable)),
    sum(out$YLD_attributable, na.rm = TRUE)
  ))

  return(out)
}


# ==============================================================================
# MORTALITY-BASED RELATIVE RISK FOR AMR BURDEN
# ==============================================================================
#
# Mixed-effects logistic regression per antibiotic class per pathogen:
#
#   logit(P(death_i)) = beta_0 + beta_1*Resistance_ci + beta_2*Age_i + beta_3*Sex_i
#                     + beta_4*HAI_i + beta_5*ICU_i + beta_6*Comorbidity_i
#                     + u_facility   [u_facility ~ N(0, sigma^2)]
#
#   OR_death_kc = exp(beta_1) : odds ratio for death, resistant vs susceptible,
#                controlling for age, sex, infection acquisition route,
#                disease severity (ICU), comorbidity, and facility clustering.
#
# Pipeline:
#   M1. derive_infection_type_for_mortality()  -- HAI/CAI with death-specific
#                                                 data-quality reporting
#   M2. .derive_icu_binary()                  -- "ever in ICU" binary per patient
#   M3. .encode_comorbidity_mortality()       -- standardise comorbidity column
#   M4. .check_hai_icu_collinearity()         -- collinearity guard (phi coeff)
#   M5. .build_class_resistance_wide()        -- binary class matrix (reused)
#   fit_mortality_rr_logistic()               -- main glmer loop
#
# References:
#   Antimicrobial Resistance Collaborators. Lancet. 2022.


# -- Step M1 -------------------------------------------------------------------

#' Derive Infection Type (HAI / CAI) for Mortality RR Model
#'
#' Variant of \code{derive_infection_type()} tailored for the mortality
#' relative-risk model. Key differences from the LOS version:
#' \enumerate{
#'   \item Processes \strong{all} patients (not just discharged), because the
#'         mortality model requires both dead and surviving patients.
#'   \item Recognises explicit HAI/CAI labels first (and common synonyms such
#'         as "hospital-acquired infection" / "community acquired infection").
#'         Date-gap derivation is applied only when the label is absent or
#'         ambiguous.
#'   \item Performs a dedicated death-date data-quality check: patients whose
#'         \code{final_outcome} equals \code{death_value} but whose outcome
#'         date is missing -- yet admission or culture date is present -- are
#'         listed. HAI/CAI derivation is unaffected (it uses admission vs
#'         culture date), but the missing death date is surfaced for review.
#' }
#'
#' @param data Data frame.
#' @param infection_type_col Character. Raw infection type column.
#'   Default \code{"type_of_infection"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_culture_col Character. Default
#'   \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_date_col Character. Date of final outcome (death date
#'   for deceased patients). Default \code{"final_outcome_date"}.
#' @param death_value Character. Value in \code{final_outcome_col} that
#'   indicates death. Default \code{"Death"}.
#' @param hai_threshold_hours Numeric. Gap threshold in hours between
#'   admission and culture date. Default \code{48}.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#'
#' @return \code{data} with column \code{infection_type_derived}
#'   (\code{"HAI"} / \code{"CAI"} / \code{"Not Known"}).
#'   Attribute \code{"missing_death_date_patients"} is a data frame of rows
#'   where death is confirmed but the outcome date is absent.
#' @export
