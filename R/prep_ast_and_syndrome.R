# prep_ast_and_syndrome.R
# AST value standardization, contaminant flagging, MDR/XDR classification,
# syndrome/infection-type derivation, and antibiotic-level collapse functions


#' Standardize Susceptibility Values
#'
#' Maps various susceptibility/resistance representations to standard
#' "S" (Susceptible), "R" (Resistant), "I" (Intermediate) values.
#'
#' @param data Data frame containing susceptibility column
#' @param col Character. Name of susceptibility column. Default "antibiotic_value".
#'
#' @return Data frame with standardized susceptibility values and added column
#'   \code{antibiotic_value_std}
#'
#' @export
prep_standardize_ast_values <- function(data, col = "antibiotic_value") {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col))
  }

  data$antibiotic_value_std <- dplyr::case_when(
    toupper(data[[col]]) %in% c("R", "RESISTANT", "RES") ~ "R",
    toupper(data[[col]]) %in% c("S", "SENSITIVE", "SUSCEPTIBLE", "SUS") ~ "S",
    toupper(data[[col]]) %in% c("I", "INTERMEDIATE", "INT") ~ "I",
    TRUE ~ NA_character_
  )

  na_count <- sum(is.na(data$antibiotic_value_std))
  if (na_count > 0) {
    warning(sprintf(
      "Susceptibility column: %d values could not be standardized to S/R/I",
      na_count
    ))
  }

  # Summary
  susc_counts <- table(data$antibiotic_value_std, useNA = "ifany")
  message("Susceptibility distribution:")
  print(susc_counts)

  return(data)
}


#' Clean Antibiotic Susceptibility Values
#'
#' Extracts clean S/I/R values from messy antibiotic result columns.
#' Handles common patterns like "S (HIGH LEVEL)", "R   Escherichia coli", etc.
#'
#' @param data Data frame with antibiotic susceptibility data
#' @param value_col Character. Column name containing susceptibility values.
#'   Default "antibiotic_value".
#' @param strict Logical. If TRUE, only accept S/I/R values. If FALSE, attempt
#'   to parse from messy strings. Default FALSE.
#'
#' @return Data frame with cleaned antibiotic_value column
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   antibiotic_value = c("S", "R   E. coli", "S (HIGH LEVEL)", "I")
#' )
#' prep_clean_ast_values(data)
#' }
prep_clean_ast_values <- function(data,
                                  value_col = "antibiotic_value",
                                  strict = FALSE) {
  if (!value_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", value_col))
  }

  n_before <- nrow(data)
  n_unique_before <- dplyr::n_distinct(data[[value_col]], na.rm = TRUE)

  message(sprintf(
    "Cleaning antibiotic values: %d unique values found",
    n_unique_before
  ))

  # Create a cleaned column
  data$temp_value_clean <- trimws(as.character(data[[value_col]]))

  if (strict) {
    # Strict mode: only accept exact S/I/R
    data[[value_col]] <- ifelse(
      data$temp_value_clean %in% c("S", "I", "R"),
      data$temp_value_clean,
      NA_character_
    )
  } else {
    # Lenient mode: extract S/I/R from messy strings
    data[[value_col]] <- sapply(data$temp_value_clean, function(val) {
      if (is.na(val) || val == "") {
        return(NA_character_)
      }

      # Convert to uppercase for matching
      val_upper <- toupper(val)

      # Extract first character that is S, I, or R
      # This handles cases like:
      # - "R   Escherichia coli" -> "R"
      # - "S (HIGH LEVEL)" -> "S"
      # - "I (HIGH LEVEL SYNERGY)" -> "I"

      # Check if first character is S/I/R
      first_char <- substr(val_upper, 1, 1)
      if (first_char %in% c("S", "I", "R")) {
        return(first_char)
      }

      # Check if string contains S, I, or R anywhere
      if (grepl("^R\\s", val_upper) || grepl("^R\\(", val_upper) || grepl("^R$", val_upper)) {
        return("R")
      }
      if (grepl("^S\\s", val_upper) || grepl("^S\\(", val_upper) || grepl("^S$", val_upper)) {
        return("S")
      }
      if (grepl("^I\\s", val_upper) || grepl("^I\\(", val_upper) || grepl("^I$", val_upper)) {
        return("I")
      }

      # Check for lowercase s, i, r
      if (val == "s") {
        return("S")
      }
      if (val == "i") {
        return("I")
      }
      if (val == "r") {
        return("R")
      }

      # Special cases
      if (grepl("resistant", val_upper)) {
        return("R")
      }
      if (grepl("susceptible", val_upper) || grepl("sensitive", val_upper)) {
        return("S")
      }
      if (grepl("intermediate", val_upper)) {
        return("I")
      }

      # Cannot parse - return NA
      return(NA_character_)
    }, USE.NAMES = FALSE)
  }

  # Clean up temp column
  data$temp_value_clean <- NULL

  n_unique_after <- dplyr::n_distinct(data[[value_col]], na.rm = TRUE)
  n_na <- sum(is.na(data[[value_col]]))

  message(sprintf(
    "Cleaned: %d unique values -> %d (S/I/R)",
    n_unique_before, n_unique_after
  ))

  if (n_na > 0) {
    message(sprintf(
      "[!] %d values could not be parsed (%.1f%%)",
      n_na, 100 * n_na / n_before
    ))
  }

  # Show distribution
  value_dist <- table(data[[value_col]], useNA = "ifany")
  message("\nValue distribution:")
  print(value_dist)

  return(data)
}


#' Recode Intermediate (I) Susceptibility Values
#'
#' Converts "I" (Intermediate) values to either "S" or "R" based on antibiotic type.
#' For Colistin: I -> S (following clinical guidelines)
#' For all other antibiotics: I -> R (conservative approach for surveillance)
#'
#' @param data Data frame with antibiotic susceptibility data
#' @param antibiotic_col Character. Column name with antibiotic names.
#'   Default "antibiotic_name".
#' @param value_col Character. Column name with S/I/R values.
#'   Default "antibiotic_value".
#' @param colistin_to_s Logical. Convert Colistin I to S. Default TRUE.
#' @param others_to_r Logical. Convert other antibiotics' I to R. Default TRUE.
#'
#' @return Data frame with recoded intermediate values
#' @export
#'
#' @examples
#' \dontrun{
#' # Recode all I values according to standard rules
#' data <- prep_recode_intermediate_ast(data)
#'
#' # Only recode Colistin
#' data <- prep_recode_intermediate_ast(data, others_to_r = FALSE)
#' }
prep_recode_intermediate_ast <- function(data,
                                         antibiotic_col = "antibiotic_name",
                                         value_col = "antibiotic_value",
                                         colistin_to_s = TRUE,
                                         others_to_r = TRUE) {
  # Input validation
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", antibiotic_col))
  }
  if (!value_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", value_col))
  }

  # Count I values before recoding
  n_intermediate_before <- sum(data[[value_col]] == "I", na.rm = TRUE)

  if (n_intermediate_before == 0) {
    message("No intermediate (I) values found. Returning data unchanged.")
    return(data)
  }

  # Create working copy
  data$temp_value <- data[[value_col]]
  data$temp_antibiotic <- tolower(trimws(data[[antibiotic_col]]))

  # Recode Colistin I -> S
  if (colistin_to_s) {
    colistin_i_idx <- which(
      data$temp_value == "I" &
        grepl("colistin", data$temp_antibiotic, ignore.case = TRUE)
    )

    if (length(colistin_i_idx) > 0) {
      data$temp_value[colistin_i_idx] <- "S"
      message(sprintf("Recoded %d Colistin I -> S", length(colistin_i_idx)))
    }
  }

  # Recode all other antibiotics I -> R
  if (others_to_r) {
    other_i_idx <- which(
      data$temp_value == "I" &
        !grepl("colistin", data$temp_antibiotic, ignore.case = TRUE)
    )

    if (length(other_i_idx) > 0) {
      data$temp_value[other_i_idx] <- "R"
      message(sprintf("Recoded %d non-Colistin I -> R", length(other_i_idx)))
    }
  }

  # Update original column
  data[[value_col]] <- data$temp_value

  # Clean up
  data$temp_value <- NULL
  data$temp_antibiotic <- NULL

  # Summary
  n_intermediate_after <- sum(data[[value_col]] == "I", na.rm = TRUE)
  message(sprintf(
    "Intermediate values: %d -> %d",
    n_intermediate_before, n_intermediate_after
  ))

  return(data)
}


#' Get Contaminant List from Reference File
#'
#' Loads contaminant organisms from common_commensals.csv and filters
#' by syndrome and/or specimen type.
#'
#' @param syndrome Character. Optional syndrome name to filter by
#'   (e.g., "Bloodstream infections", "Urinary tract infections").
#' @param specimen_type Character. Optional specimen type to filter by
#'   (e.g., "Blood culture", "Urine culture").
#' @param return_all Logical. If TRUE, returns all contaminants from all
#'   syndromes. Default FALSE.
#'
#' @return Character vector of contaminant organism names (lowercase)
#'
#' @export
#'
#' @examples
#' # Get all blood culture contaminants
#' prep_get_contaminant_list(syndrome = "Bloodstream infections")
#'
#' # Get all contaminants for a specific specimen type
#' prep_get_contaminant_list(specimen_type = "Blood culture")
#'
#' # Get all contaminants across all syndromes
#' prep_get_contaminant_list(return_all = TRUE)
prep_get_contaminant_list <- function(syndrome = NULL,
                                      specimen_type = NULL,
                                      return_all = FALSE) {
  # Load common commensals reference file
  commensals_path <- find_extdata_file("common_commensals.csv")

  if (commensals_path == "" || !file.exists(commensals_path)) {
    warning("common_commensals.csv not found in inst/extdata/. Returning empty list.")
    return(character(0))
  }

  # Read the CSV file
  commensals_data <- readr::read_csv(
    commensals_path,
    show_col_types = FALSE
  )

  # Filter by syndrome if provided
  if (!is.null(syndrome) && !return_all) {
    commensals_data <- commensals_data %>%
      dplyr::filter(Syndrome == syndrome)
  }

  # Filter by specimen type if provided
  if (!is.null(specimen_type) && !return_all) {
    commensals_data <- commensals_data %>%
      dplyr::filter(`Type of culture/specimen` == specimen_type)
  }

  # Extract and parse the contaminant lists
  contaminants <- commensals_data %>%
    dplyr::pull(`Common commensals`) %>%
    stringr::str_split(";\\s*") %>%
    unlist() %>%
    stringr::str_trim() %>%
    unique()

  # Remove any empty strings
  contaminants <- contaminants[contaminants != ""]

  # Generate flexible patterns for each contaminant
  # This allows matching variations like "Staph" -> "Staphylococcus"
  contaminant_patterns <- lapply(contaminants, function(name) {
    name_lower <- tolower(name)

    # Split into genus and species
    parts <- strsplit(name_lower, "\\s+")[[1]]

    patterns <- c(name_lower) # Original full name

    if (length(parts) >= 1) {
      genus <- parts[1]

      # Only add genus as a standalone pattern for genus-level entries
      # (e.g. "Micrococcus species", "Aerococcus spp.") -- not for specific species
      # like "Staphylococcus epidermidis", which would incorrectly match S. aureus
      is_genus_level <- length(parts) == 1 ||
        grepl("^sp(p\\.?|ecies)$", parts[2], ignore.case = TRUE)

      # Abbreviation for combining with species: "^s\." matches "s. epidermidis"
      genus_abbrev_base <- paste0("^", substr(genus, 1, 1), "\\.")
      # Standalone abbreviation requires dot + whitespace to avoid matching
      # full genus names: "^a\.\s" matches "A. spp." but NOT "Acinetobacter spp."
      genus_abbrev_standalone <- paste0("^", substr(genus, 1, 1), "\\.\\s")

      if (is_genus_level) {
        patterns <- c(patterns, genus, genus_abbrev_standalone)
      } else {
        patterns <- c(patterns, genus_abbrev_standalone)
      }

      # If there's a species, add genus+species combinations
      if (length(parts) >= 2) {
        species <- parts[2]

        # Only add species as standalone if it is a real binomial species name:
        # - Skip "spp." / "species" -- these match any organism with that suffix
        # - Skip when genus contains a hyphen (descriptor format like "coagulase-negative")
        #   because the second word is another genus name (e.g., "Staphylococcus"),
        #   not a true species, and would match unrelated organisms
        is_real_species <- !grepl("^sp(p\\.?|ecies)$", species, ignore.case = TRUE) &&
          !grepl("-", genus)

        if (is_real_species) {
          patterns <- c(patterns, species)
        }
        patterns <- c(
          patterns,
          paste0(genus_abbrev_base, "\\s*", species), # e.g., "^s\.\\s*epidermidis"
          paste(genus, species) # Full name
        )
      }

      # Handle common abbreviations
      if (grepl("staphylococcus", genus)) {
        if (is_genus_level) {
          # Genus-level entry: "staph" alone is safe (e.g., "Staphylococcus spp.")
          patterns <- c(patterns, "staph", "coag.*neg", "coagulase.*negative")
        } else {
          # Species-specific entry: bind "staph" to the species to avoid matching S. aureus
          patterns <- c(patterns, paste0("staph.*", parts[2]), "coag.*neg", "coagulase.*negative")
        }
      }
      if (grepl("streptococcus", genus)) {
        patterns <- c(patterns, "strep")
      }
      if (grepl("escherichia", genus)) {
        patterns <- c(patterns, "e\\.?\\s*coli")
      }
      if (grepl("pseudomonas", genus)) {
        patterns <- c(patterns, "pseudo")
      }
      if (grepl("klebsiella", genus)) {
        patterns <- c(patterns, "kleb")
      }
      if (grepl("acinetobacter", genus)) {
        patterns <- c(patterns, "acin")
      }
      if (grepl("enterococcus", genus)) {
        patterns <- c(patterns, "entero")
      }
      if (grepl("corynebacterium", genus)) {
        patterns <- c(patterns, "coryno", "diphtheroids?")
      }
      if (grepl("bacillus", genus)) {
        patterns <- c(patterns, "bacil")
      }
      if (grepl("micrococcus", genus)) {
        patterns <- c(patterns, "micro")
      }
      if (grepl("cutibacterium|propionibacterium", genus)) {
        patterns <- c(patterns, "propioni", "cuti", "p\\.?\\s*acnes")
      }
      if (grepl("lactobacillus", genus)) {
        patterns <- c(patterns, "lacto")
      }
    }

    # Return unique patterns
    list(
      original = name,
      patterns = unique(patterns)
    )
  })

  # Return structure with both original names and matching patterns
  result <- list(
    names = contaminants,
    patterns = contaminant_patterns
  )

  return(result)
}


#' Check if Organism is a Contaminant
#'
#' Checks if a given organism name matches any known contaminant for a
#' specific syndrome or specimen type using flexible pattern matching.
#'
#' @param organism_name Character vector. Organism name(s) to check.
#' @param syndrome Character. Optional syndrome name to filter contaminants.
#' @param specimen_type Character. Optional specimen type to filter contaminants.
#'
#' @return Logical vector indicating if each organism is a contaminant
#'
#' @export
#'
#' @examples
#' # Check if organism is a blood culture contaminant
#' prep_is_contaminant("Staph epidermidis", syndrome = "Bloodstream infections")
#' prep_is_contaminant("E. coli", syndrome = "Bloodstream infections")
#'
#' # Check multiple organisms
#' prep_is_contaminant(c("Staph", "Klebsiella"), specimen_type = "Blood culture")
prep_is_contaminant <- function(organism_name,
                                syndrome = NULL,
                                specimen_type = NULL) {
  # Get contaminant list with patterns
  contaminant_data <- prep_get_contaminant_list(
    syndrome = syndrome,
    specimen_type = specimen_type
  )

  # Check each organism name against all contaminant patterns
  sapply(organism_name, function(org) {
    if (is.na(org) || org == "") {
      return(FALSE)
    }

    org_lower <- tolower(trimws(org))

    # Check against each contaminant's patterns
    for (contam_info in contaminant_data$patterns) {
      for (pattern in contam_info$patterns) {
        if (grepl(pattern, org_lower, perl = TRUE)) {
          return(TRUE)
        }
      }
    }

    return(FALSE)
  }, USE.NAMES = FALSE)
}


#' Flag Contaminant Organisms
#'
#' Identifies likely contaminant organisms using multi-path logic.
#'
#' @param data Data frame
#' @param method Character. "auto" (try all methods), "device_based",
#'   "heuristic", "provided". Default "auto".
#' @param organism_col Character. Normalized organism column.
#' @param specimen_col Character. Specimen type column.
#'
#' @return Data frame with is_contaminant, contaminant_confidence,
#'   contaminant_method columns
#' @export
prep_flag_contaminants <- function(data,
                                   method = "auto",
                                   organism_col = "organism_normalized",
                                   specimen_col = "specimen_type") {
  # Check if required columns exist
  if (!organism_col %in% names(data)) {
    message(sprintf("[!] Column '%s' not found. Skipping contaminant flagging.", organism_col))
    data$is_contaminant <- FALSE
    data$contaminant_confidence <- "insufficient_data"
    data$contaminant_method <- "skipped"
    return(data)
  }

  if (!specimen_col %in% names(data)) {
    data[[specimen_col]] <- NA_character_
  }

  # Initialize
  data$is_contaminant <- FALSE
  data$contaminant_confidence <- "unknown"
  data$contaminant_method <- NA_character_

  # Path C: Use provided flag if exists
  if ("pathogen_contaminant" %in% names(data)) {
    data$is_contaminant <- data$pathogen_contaminant == 1
    data$contaminant_method <- "provided"
    data$contaminant_confidence <- "high"
    message("Using provided contaminant flags")
    return(data)
  }

  # Path A: Device-based (if device data available)
  if (all(c("device_inserted", "device_insertion_date", "date_of_culture") %in% names(data))) {
    data <- data %>%
      dplyr::mutate(
        days_since_device = as.numeric(
          difftime(date_of_culture, device_insertion_date, units = "days")
        )
      )

    # CoNS within 48h of line insertion = likely contaminant
    data <- data %>%
      dplyr::mutate(
        is_contaminant = dplyr::case_when(
          !!rlang::sym(organism_col) %in% c(
            "staphylococcus epidermidis",
            "staphylococcus haemolyticus"
          ) &
            !!rlang::sym(specimen_col) == "blood" &
            days_since_device <= 2 ~ TRUE,
          TRUE ~ is_contaminant
        ),
        contaminant_confidence = dplyr::if_else(
          !!rlang::sym(organism_col) %in% c(
            "staphylococcus epidermidis",
            "staphylococcus haemolyticus"
          ) &
            !!rlang::sym(specimen_col) == "blood",
          "high",
          contaminant_confidence
        ),
        contaminant_method = dplyr::if_else(
          is.na(contaminant_method),
          "device_based",
          contaminant_method
        )
      )

    n_device_based <- sum(data$contaminant_method == "device_based", na.rm = TRUE)
    message(sprintf("Device-based contaminants: %d", n_device_based))
  }

  # Path B: Heuristic (organism + specimen type)
  contam_blood <- prep_get_contaminant_list(syndrome = "Bloodstream infections")
  contam_urine <- prep_get_contaminant_list(syndrome = "Urinary tract infections / pyelonephritis")

  # Store previous contaminant flags
  prev_contaminant <- data$is_contaminant
  prev_method <- data$contaminant_method

  data <- data %>%
    dplyr::mutate(
      # Create temporary safe specimen column
      .specimen_lower = tolower(!!rlang::sym(specimen_col)),
      is_contaminant = dplyr::case_when(
        # Blood contaminants -- grepl handles "Blood-peripheral", "Blood-central catheter" etc.
        grepl("blood", .specimen_lower) &
          !!rlang::sym(organism_col) %in% contam_blood$names ~ TRUE,

        # Urine contaminants -- grepl handles "Urine culture", "Urine" etc.
        grepl("urine", .specimen_lower) &
          !!rlang::sym(organism_col) %in% contam_urine$names ~ TRUE,

        # Keep existing TRUE
        prev_contaminant == TRUE ~ TRUE,
        TRUE ~ FALSE
      ),
      contaminant_confidence = dplyr::case_when(
        is_contaminant & !is.na(prev_method) & prev_method == "device_based" ~ "high",
        is_contaminant ~ "medium",
        TRUE ~ "low"
      ),
      contaminant_method = dplyr::if_else(
        is.na(prev_method) | prev_method == "",
        "heuristic",
        prev_method
      ),
      # Remove temporary column
      .specimen_lower = NULL
    )

  # Summary
  n_contam <- sum(data$is_contaminant)
  n_high <- sum(data$contaminant_confidence == "high", na.rm = TRUE)
  n_medium <- sum(data$contaminant_confidence == "medium", na.rm = TRUE)

  message(sprintf(
    "Contaminants flagged: %d total (high: %d, medium: %d)",
    n_contam, n_high, n_medium
  ))

  return(data)
}


#' Classify Infection-Related Mortality
#'
#' Determines if death was related to infection using date window logic.
#' Implements proxy logic when dates are missing.
#'
#' @param data Data frame
#' @param outcome_col Character. Outcome column. Default "final_outcome".
#' @param event_date_col Character. Event/culture date. Default "date_of_culture".
#' @param outcome_date_col Character. Outcome date. Default "date_of_final_outcome".
#' @param window Numeric. Days after event to classify death as infection-related.
#'   Default 14.
#'
#' @return Data frame with mortality_infection, mortality_method,
#'   mortality_confidence columns
#' @export
prep_classify_mortality <- function(data,
                                    outcome_col = "final_outcome",
                                    event_date_col = "date_of_culture",
                                    outcome_date_col = "date_of_final_outcome",
                                    window = 14) {
  # Check columns
  if (!outcome_col %in% names(data)) {
    stop(sprintf("Outcome column '%s' not found", outcome_col))
  }

  # Initialize
  data$mortality_infection <- "No"
  data$mortality_method <- NA_character_
  data$mortality_confidence <- NA_character_

  # Path A: Full date information available (HIGH CONFIDENCE)
  if (all(c(event_date_col, outcome_date_col) %in% names(data))) {
    data <- data %>%
      dplyr::mutate(
        gap_days = as.numeric(
          difftime(!!rlang::sym(outcome_date_col),
            !!rlang::sym(event_date_col),
            units = "days"
          )
        ),
        within_window = !is.na(gap_days) & gap_days >= 0 & gap_days <= window,
        mortality_infection = dplyr::case_when(
          !!rlang::sym(outcome_col) == "Died" & within_window ~ "Yes",
          !!rlang::sym(outcome_col) == "Died" & !within_window ~ "No",
          TRUE ~ "No"
        ),
        mortality_method = dplyr::if_else(
          !!rlang::sym(outcome_col) == "Died" & !is.na(gap_days),
          "date_calculated",
          NA_character_
        ),
        mortality_confidence = dplyr::case_when(
          mortality_method == "date_calculated" ~ "high",
          TRUE ~ NA_character_
        )
      )

    n_high_conf <- sum(data$mortality_confidence == "high", na.rm = TRUE)
    message(sprintf(
      "Classified mortality using dates (%d-day window): %d high confidence",
      window, n_high_conf
    ))
  } else {
    message("Date columns not available for mortality classification")
  }

  # Path B: PROXY - Only outcome available (LOW CONFIDENCE)
  data <- data %>%
    dplyr::mutate(
      mortality_infection = dplyr::case_when(
        # Keep high confidence classifications
        !is.na(mortality_confidence) ~ mortality_infection,

        # Proxy: Died but no dates -> mark as "Possible"
        !!rlang::sym(outcome_col) == "Died" ~ "Possible",
        TRUE ~ "No"
      ),
      mortality_method = dplyr::case_when(
        !is.na(mortality_method) ~ mortality_method,
        !!rlang::sym(outcome_col) == "Died" ~ "proxy_outcome_only",
        TRUE ~ NA_character_
      ),
      mortality_confidence = dplyr::case_when(
        !is.na(mortality_confidence) ~ mortality_confidence,
        mortality_method == "proxy_outcome_only" ~ "low",
        TRUE ~ NA_character_
      )
    )

  # Summary
  mortality_summary <- table(
    Method = data$mortality_method,
    Result = data$mortality_infection,
    useNA = "ifany"
  )

  message("\nMortality classification summary:")
  print(mortality_summary)

  n_proxy <- sum(data$mortality_method == "proxy_outcome_only", na.rm = TRUE)
  if (n_proxy > 0) {
    message(sprintf(
      "\n[!] Warning: %d deaths classified using PROXY (dates missing). Low confidence.",
      n_proxy
    ))
  }

  return(data)
}


#' Classify MDR (Multidrug Resistant)
#'
#' Classifies isolates as MDR using Magiorakos 2012 criteria.
#'
#' @param data Data frame
#' @param definition Character. "Magiorakos" or "WHO". Default "Magiorakos".
#' @param organism_group_col Character. Organism group column for
#'   pathogen-specific thresholds.
#'
#' @return Data frame with mdr, mdr_confidence, mdr_method,
#'   n_resistant_categories columns
#' @export
#' @references
#' Magiorakos AP et al. Clin Microbiol Infect. 2012;18(3):268-281.
prep_classify_mdr <- function(data,
                              definition = "Magiorakos",
                              organism_group_col = "org_group") {
  # Requires class-level resistance data
  if (!"class_result_event" %in% names(data)) {
    stop("Must run prep_collapse_class_level() before MDR classification")
  }

  # Get thresholds
  thresholds <- get_magiorakos_thresholds()

  # Count resistant categories per event
  resistant_counts <- data %>%
    dplyr::filter(class_result_event == "R") %>%
    dplyr::group_by(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::summarise(
      n_resistant_categories = dplyr::n_distinct(antibiotic_class),
      resistant_categories = paste(unique(antibiotic_class), collapse = "; "),
      .groups = "drop"
    )

  # Total categories tested
  total_tested <- data %>%
    dplyr::group_by(event_id) %>%
    dplyr::summarise(
      n_total_categories = dplyr::n_distinct(antibiotic_class),
      .groups = "drop"
    )

  # Merge
  mdr_data <- resistant_counts %>%
    dplyr::left_join(total_tested, by = "event_id") %>%
    dplyr::left_join(thresholds, by = stats::setNames("organism_group", organism_group_col))

  # Apply MDR threshold
  mdr_data <- mdr_data %>%
    dplyr::mutate(
      mdr_threshold = dplyr::coalesce(mdr_threshold, 3), # Default to 3 if no match
      mdr = n_resistant_categories >= mdr_threshold,
      mdr_confidence = dplyr::case_when(
        n_total_categories >= 8 ~ "high",
        n_total_categories >= 5 ~ "medium",
        n_total_categories >= 3 ~ "low",
        TRUE ~ "insufficient_data"
      ),
      mdr_method = definition
    )

  # Join back to main data
  data <- data %>%
    dplyr::left_join(
      mdr_data %>% dplyr::select(
        event_id, mdr, mdr_confidence, mdr_method,
        n_resistant_categories, resistant_categories
      ),
      by = "event_id"
    ) %>%
    dplyr::mutate(mdr = tidyr::replace_na(mdr, FALSE))

  # Summary
  n_mdr <- sum(data$mdr & data$mdr_confidence != "insufficient_data", na.rm = TRUE)
  n_total <- dplyr::n_distinct(data$event_id)

  message(sprintf(
    "MDR classification (%s): %d/%d events (%.1f%%)",
    definition, n_mdr, n_total, 100 * n_mdr / n_total
  ))

  return(data)
}


#' Classify XDR (Extensively Drug Resistant)
#'
#' Classifies isolates as XDR using Magiorakos 2012 criteria.
#'
#' @param data Data frame
#' @param definition Character. "Magiorakos" or "WHO". Default "Magiorakos".
#' @param organism_group_col Character. Organism group column.
#'
#' @return Data frame with xdr, xdr_confidence, xdr_method columns
#' @export
#' @references
#' Magiorakos AP et al. Clin Microbiol Infect. 2012;18(3):268-281.
prep_classify_xdr <- function(data,
                              definition = "Magiorakos",
                              organism_group_col = "org_group") {
  # Requires MDR to be run first
  if (!"mdr" %in% names(data)) {
    message("Running MDR classification first...")
    data <- prep_classify_mdr(data, definition, organism_group_col)
  }

  # Get thresholds
  thresholds <- get_magiorakos_thresholds()

  # Count susceptible categories
  susceptible_counts <- data %>%
    dplyr::filter(class_result_event == "S") %>%
    dplyr::group_by(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::summarise(
      n_susceptible_categories = dplyr::n_distinct(antibiotic_class),
      .groups = "drop"
    )

  # Get total categories per organism
  xdr_data <- data %>%
    dplyr::distinct(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::left_join(susceptible_counts, by = c("event_id", organism_group_col)) %>%
    dplyr::left_join(thresholds, by = stats::setNames("organism_group", organism_group_col)) %>%
    dplyr::mutate(
      n_susceptible_categories = tidyr::replace_na(n_susceptible_categories, 0),
      # XDR: susceptible to <=2 categories
      xdr = n_susceptible_categories <= 2,
      xdr_confidence = dplyr::case_when(
        !is.na(total_categories) ~ "high", # Known pathogen
        TRUE ~ "medium"
      ),
      xdr_method = definition
    )

  # Join back
  data <- data %>%
    dplyr::left_join(
      xdr_data %>% dplyr::select(event_id, xdr, xdr_confidence, xdr_method),
      by = "event_id"
    ) %>%
    dplyr::mutate(xdr = tidyr::replace_na(xdr, FALSE))

  # Summary
  n_xdr <- sum(data$xdr, na.rm = TRUE)
  n_total <- dplyr::n_distinct(data$event_id)

  message(sprintf(
    "XDR classification (%s): %d/%d events (%.1f%%)",
    definition, n_xdr, n_total, 100 * n_xdr / n_total
  ))

  return(data)
}


#' Derive HAI/CAI Infection Type
#'
#' Infers Community-Acquired (CAI) vs Hospital-Acquired (HAI) infection
#' using admission-culture date gap.
#'
#' @param data Data frame
#' @param infection_type_col Character. Infection type column.
#'   Default "infection_type".
#' @param admission_col Character. Admission date. Default "date_of_admission".
#' @param culture_col Character. Culture date. Default "date_of_culture".
#' @param hai_cutoff Numeric. Days after admission to classify as HAI.
#'   Default 2 (48 hours).
#' @param overwrite Logical. Recalculate even if present. Default FALSE.
#'
#' @return Data frame with infection_type enriched
#' @export
#'
#' @examples
#' \dontrun{
#' # Default 2-day cutoff
#' data_enriched <- prep_derive_hai_cai(data)
#'
#' # 3-day cutoff
#' data_enriched <- prep_derive_hai_cai(data, hai_cutoff = 3)
#' }
prep_derive_hai_cai <- function(data,
                                infection_type_col = "infection_type",
                                admission_col = "date_of_admission",
                                culture_col = "date_of_culture",
                                hai_cutoff = 2,
                                overwrite = FALSE) {
  # Check columns
  has_infection_type <- infection_type_col %in% names(data)
  has_admission <- admission_col %in% names(data)
  has_culture <- culture_col %in% names(data)

  if (!has_infection_type) {
    data[[infection_type_col]] <- NA_character_
  }

  if (!has_admission || !has_culture) {
    message(sprintf(
      "[!] Cannot infer infection type: missing '%s' or '%s'",
      admission_col, culture_col
    ))
    return(data)
  }

  n_before_missing <- sum(is.na(data[[infection_type_col]]))

  message(sprintf(
    "Inferring infection type using %d-day HAI cutoff...",
    hai_cutoff
  ))

  # Infer infection type
  data <- data %>%
    dplyr::mutate(
      days_to_culture = dplyr::case_when(
        !is.na(!!rlang::sym(admission_col)) & !is.na(!!rlang::sym(culture_col)) ~
          as.numeric(difftime(!!rlang::sym(culture_col), !!rlang::sym(admission_col), units = "days")),
        TRUE ~ NA_real_
      ),
      inferred_type = dplyr::case_when(
        !is.na(days_to_culture) & days_to_culture >= hai_cutoff ~ "HAI",
        !is.na(days_to_culture) & days_to_culture < hai_cutoff ~ "CAI",
        TRUE ~ NA_character_
      ),
      !!infection_type_col := dplyr::case_when(
        overwrite & !is.na(inferred_type) ~ inferred_type,
        is.na(!!rlang::sym(infection_type_col)) & !is.na(inferred_type) ~ inferred_type,
        TRUE ~ !!rlang::sym(infection_type_col)
      ),
      infection_type_method = dplyr::case_when(
        !is.na(inferred_type) ~ sprintf("inferred_%dday_cutoff", hai_cutoff),
        !is.na(!!rlang::sym(infection_type_col)) ~ "provided",
        TRUE ~ NA_character_
      ),
      infection_type_confidence = dplyr::case_when(
        infection_type_method == "provided" ~ "high",
        !is.na(inferred_type) ~ "medium",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-inferred_type, -days_to_culture)

  n_after_missing <- sum(is.na(data[[infection_type_col]]))
  n_enriched <- n_before_missing - n_after_missing

  if (n_enriched > 0) {
    message(sprintf(
      "Enriched infection_type: %d rows filled",
      n_enriched
    ))
  }

  # Summary
  type_summary <- data %>%
    dplyr::count(!!rlang::sym(infection_type_col), infection_type_method) %>%
    dplyr::arrange(dplyr::desc(n))

  message("\nInfection type distribution:")
  print(type_summary)

  return(data)
}


#' Collapse to Antibiotic Level (OPTIONAL - Run When YOU Decide)
#'
#' **IMPORTANT**: This function removes duplicate tests. Only run this when
#' you have reviewed your data and decided to collapse duplicates.
#'
#' Aggregates multiple test results for the same organism-antibiotic combination
#' within an event. Uses "any R -> R" logic where resistance in any test
#' results in resistant classification.
#'
#' **When to use**: After you've cleaned and normalized data, if you have
#' multiple tests for the same patient-organism-antibiotic and want one
#' result per combination.
#'
#' **What it removes**: Duplicate rows based on event_id + organism + antibiotic
#'
#' @param data Data frame with susceptibility results
#' @param event_col Character. Event ID column. Default "event_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_normalized".
#' @param susceptibility_col Character. Susceptibility column (S/I/R).
#'   Default "antibiotic_value".
#' @param aggregation_rule Character. Rule for aggregation: "any_R" (default,
#'   any R -> R), "most_resistant" (R > I > S), "most_common" (mode).
#'
#' @return Aggregated data frame (one row per event-organism-antibiotic)
#' @export
#'
#' @examples
#' \dontrun{
#' # Check for duplicates first
#' data %>%
#'   group_by(event_id, organism_normalized, antibiotic_normalized) %>%
#'   filter(n() > 1)
#'
#' # Then decide to collapse using any R -> R
#' collapsed <- prep_collapse_antibiotic_level(data)
#'
#' # Or use most common result
#' collapsed <- prep_collapse_antibiotic_level(
#'   data,
#'   aggregation_rule = "most_common"
#' )
#' }
prep_collapse_antibiotic_level <- function(data,
                                           event_col = "event_id",
                                           organism_col = "organism_normalized",
                                           antibiotic_col = "antibiotic_normalized",
                                           susceptibility_col = "antibiotic_value",
                                           aggregation_rule = "any_R") {
  # Validate columns
  required_cols <- c(event_col, organism_col, antibiotic_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Validate rule
  valid_rules <- c("any_R", "most_resistant", "most_common")
  if (!aggregation_rule %in% valid_rules) {
    stop(sprintf("aggregation_rule must be one of: %s", paste(valid_rules, collapse = ", ")))
  }

  n_before <- nrow(data)

  message(sprintf(
    "Collapsing to antibiotic level using rule: %s",
    aggregation_rule
  ))

  # Identify duplicates
  duplicates <- data %>%
    dplyr::group_by(
      !!rlang::sym(event_col),
      !!rlang::sym(organism_col),
      !!rlang::sym(antibiotic_col)
    ) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()

  n_duplicates <- dplyr::n_distinct(
    duplicates[[event_col]],
    duplicates[[organism_col]],
    duplicates[[antibiotic_col]]
  )

  if (n_duplicates > 0) {
    message(sprintf(
      "Found %d event-organism-antibiotic combinations with multiple tests",
      n_duplicates
    ))
  }

  # Apply aggregation rule
  if (aggregation_rule == "any_R") {
    # Any R -> R
    collapsed <- data %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::arrange(dplyr::desc(!!rlang::sym(susceptibility_col))) %>%
      dplyr::mutate(
        final_result = dplyr::case_when(
          any(!!rlang::sym(susceptibility_col) == "R", na.rm = TRUE) ~ "R",
          any(!!rlang::sym(susceptibility_col) == "I", na.rm = TRUE) ~ "I",
          any(!!rlang::sym(susceptibility_col) == "S", na.rm = TRUE) ~ "S",
          TRUE ~ NA_character_
        ),
        aggregation_method = "any_R",
        n_tests_aggregated = dplyr::n()
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  } else if (aggregation_rule == "most_resistant") {
    # R > I > S hierarchy
    collapsed <- data %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::mutate(
        resistance_rank = dplyr::case_when(
          !!rlang::sym(susceptibility_col) == "R" ~ 3,
          !!rlang::sym(susceptibility_col) == "I" ~ 2,
          !!rlang::sym(susceptibility_col) == "S" ~ 1,
          TRUE ~ 0
        ),
        n_tests_aggregated = dplyr::n()
      ) %>%
      dplyr::arrange(dplyr::desc(resistance_rank)) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        final_result = !!rlang::sym(susceptibility_col),
        aggregation_method = "most_resistant"
      ) %>%
      dplyr::select(-resistance_rank) %>%
      dplyr::ungroup()
  } else if (aggregation_rule == "most_common") {
    # Mode (most frequent result)
    collapsed <- data %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::mutate(n_tests_aggregated = dplyr::n()) %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col),
        !!rlang::sym(susceptibility_col),
        .add = FALSE
      ) %>%
      dplyr::mutate(result_count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::arrange(dplyr::desc(result_count)) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        final_result = !!rlang::sym(susceptibility_col),
        aggregation_method = "most_common"
      ) %>%
      dplyr::select(-result_count) %>%
      dplyr::ungroup()
  }

  # Update susceptibility column with final result
  collapsed <- collapsed %>%
    dplyr::mutate(!!susceptibility_col := final_result) %>%
    dplyr::select(-final_result)

  n_after <- nrow(collapsed)
  n_removed <- n_before - n_after

  message(sprintf(
    "Collapsed: %d -> %d rows (%d duplicates removed)",
    n_before, n_after, n_removed
  ))

  return(collapsed)
}


#' Collapse to Class Level
#'
#' Aggregates resistance at antibiotic class level instead of individual drugs.
#' Uses "any R in class -> class R" logic.
#'
#' @param data Data frame with antibiotic class information
#' @param event_col Character. Event ID column. Default "event_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param class_col Character. Antibiotic class column. Default "antibiotic_class".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param extra_cols Character vector or NULL. Additional columns to carry
#'   through the aggregation. Default NULL.
#'
#' @return Aggregated data frame (one row per event-organism-class)
#' @export
#'
#' @examples
#' \dontrun{
#' class_level <- prep_collapse_class_level(data)
#' }
prep_collapse_class_level <- function(data,
                                      event_col = "event_id",
                                      organism_col = "organism_normalized",
                                      class_col = "antibiotic_class",
                                      susceptibility_col = "antibiotic_value",
                                      extra_cols = NULL) {
  # ---------------------------------------------------------------------------
  # Validate required columns
  # ---------------------------------------------------------------------------
  required_cols <- c(event_col, organism_col, class_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Validate optional extra columns
  if (!is.null(extra_cols)) {
    missing_extra <- setdiff(extra_cols, names(data))
    if (length(missing_extra) > 0) {
      stop(sprintf(
        "Missing extra columns requested: %s",
        paste(missing_extra, collapse = ", ")
      ))
    }
  }

  n_before <- nrow(data)
  message("Collapsing to antibiotic class level...")

  # ---------------------------------------------------------------------------
  # Build grouping variables dynamically
  # ---------------------------------------------------------------------------
  group_vars <- c(event_col, organism_col, class_col)

  # ---------------------------------------------------------------------------
  # Core aggregation
  # ---------------------------------------------------------------------------
  collapsed <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      class_resistance = dplyr::case_when(
        any(.data[[susceptibility_col]] == "R", na.rm = TRUE) ~ "R",
        any(.data[[susceptibility_col]] == "I", na.rm = TRUE) ~ "I",
        any(.data[[susceptibility_col]] == "S", na.rm = TRUE) ~ "S",
        TRUE ~ NA_character_
      ),
      n_drugs_in_class = dplyr::n(),
      n_resistant = sum(.data[[susceptibility_col]] == "R", na.rm = TRUE),
      pct_resistant_in_class = 100 * n_resistant / n_drugs_in_class,
      drugs_tested = paste(
        sort(unique(.data[["antibiotic_normalized"]])),
        collapse = "; "
      ),
      .groups = "drop"
    )

  # ---------------------------------------------------------------------------
  # Attach optional extra columns (collapsed safely)
  # ---------------------------------------------------------------------------
  if (!is.null(extra_cols)) {
    extra_summary <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(extra_cols),
          ~ paste(sort(unique(.x)), collapse = "; "),
          .names = "{.col}"
        ),
        .groups = "drop"
      )

    collapsed <- collapsed %>%
      dplyr::left_join(extra_summary, by = group_vars)
  }

  # ---------------------------------------------------------------------------
  # Metadata
  # ---------------------------------------------------------------------------
  collapsed <- collapsed %>%
    dplyr::mutate(collapse_method = "class_any_R")

  n_after <- nrow(collapsed)

  message(sprintf(
    "Collapsed: %d rows -> %d class-level rows",
    n_before, n_after
  ))

  # ---------------------------------------------------------------------------
  # Class-level resistance summary (informational)
  # ---------------------------------------------------------------------------
  class_summary <- collapsed %>%
    dplyr::count(.data[[class_col]], class_resistance) %>%
    tidyr::pivot_wider(
      names_from = class_resistance,
      values_from = n,
      values_fill = 0
    )

  message("\nClass-level resistance distribution:")
  print(class_summary)

  return(collapsed)
}


#' Create Resistance Profile
#'
#' Generates a resistance profile string for each event summarizing
#' resistance patterns. Useful for identifying common resistance phenotypes.
#'
#' @param data Data frame with resistance data
#' @param event_col Character. Event ID column. Default "event_id".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_normalized".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param format Character. Output format:
#'   - "resistant_list": List resistant drugs only (default)
#'   - "full_pattern": Full S/R pattern string
#'   - "class_summary": Resistant classes only
#' @param class_col Character. Class column (required for format = "class_summary").
#'   Default "antibiotic_class".
#'
#' @return Data frame with resistance_profile column added
#' @export
#'
#' @examples
#' \dontrun{
#' # List resistant drugs
#' data_with_profile <- prep_create_resistance_profile(data)
#'
#' # Full S/R pattern
#' data_with_profile <- prep_create_resistance_profile(
#'   data,
#'   format = "full_pattern"
#' )
#'
#' # Resistant classes only
#' data_with_profile <- prep_create_resistance_profile(
#'   data,
#'   format = "class_summary"
#' )
#' }
prep_create_resistance_profile <- function(data,
                                           event_col = "event_id",
                                           antibiotic_col = "antibiotic_normalized",
                                           susceptibility_col = "antibiotic_value",
                                           format = "resistant_list",
                                           class_col = "antibiotic_class") {
  # Validate columns
  required_cols <- c(event_col, antibiotic_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (format == "class_summary" && !class_col %in% names(data)) {
    stop(sprintf("Column '%s' required for format = 'class_summary'", class_col))
  }

  # Validate format
  valid_formats <- c("resistant_list", "full_pattern", "class_summary")
  if (!format %in% valid_formats) {
    stop(sprintf("format must be one of: %s", paste(valid_formats, collapse = ", ")))
  }

  message(sprintf("Creating resistance profiles (format: %s)...", format))

  # Format A: List resistant drugs only
  if (format == "resistant_list") {
    profiles <- data %>%
      dplyr::filter(!!rlang::sym(susceptibility_col) == "R") %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::summarise(
        resistance_profile = paste(sort(unique(!!rlang::sym(antibiotic_col))), collapse = "; "),
        n_resistant = dplyr::n_distinct(!!rlang::sym(antibiotic_col)),
        .groups = "drop"
      )

    # Add profile to original data
    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(
        resistance_profile = dplyr::coalesce(resistance_profile, "None"),
        n_resistant = dplyr::coalesce(n_resistant, 0L),
        profile_format = "resistant_list"
      )
  }

  # Format B: Full S/R pattern
  else if (format == "full_pattern") {
    profiles <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(!!rlang::sym(antibiotic_col)) %>%
      dplyr::summarise(
        resistance_profile = paste(
          sprintf(
            "%s:%s",
            !!rlang::sym(antibiotic_col),
            !!rlang::sym(susceptibility_col)
          ),
          collapse = "; "
        ),
        n_resistant = sum(!!rlang::sym(susceptibility_col) == "R", na.rm = TRUE),
        n_tested = dplyr::n(),
        .groups = "drop"
      )

    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(profile_format = "full_pattern")
  }

  # Format C: Resistant classes summary
  else if (format == "class_summary") {
    profiles <- data %>%
      dplyr::filter(!!rlang::sym(susceptibility_col) == "R") %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::summarise(
        resistance_profile = paste(sort(unique(!!rlang::sym(class_col))), collapse = "; "),
        n_resistant_classes = dplyr::n_distinct(!!rlang::sym(class_col)),
        .groups = "drop"
      )

    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(
        resistance_profile = dplyr::coalesce(resistance_profile, "None"),
        n_resistant_classes = dplyr::coalesce(n_resistant_classes, 0L),
        profile_format = "class_summary"
      )
  }

  # Summary of common profiles
  profile_freq <- data %>%
    dplyr::filter(!duplicated(!!rlang::sym(event_col))) %>%
    dplyr::count(resistance_profile, sort = TRUE) %>%
    utils::head(10)

  message("\nTop 10 resistance profiles:")
  print(profile_freq)

  return(data)
}
# select.R
# Resistance class selection using beta-lactam hierarchy and RR ranking

#' Select Resistance Class
#'
#' Selects a single resistance class per event using beta-lactam hierarchy
#' and relative risk (RR) values. Prevents double-counting in burden estimation
#' by choosing the most clinically relevant resistant class.
#'
#' Selection logic:
#' 1. Filter to resistant classes only (R)
#' 2. Apply beta-lactam hierarchy (Carbapenems > 4GC > 3GC > ...)
#' 3. Within same hierarchy rank, prioritize by RR value (higher RR first)
#' 4. If tied, select alphabetically for reproducibility
#'
#' @param data Data frame with resistance and RR information
#' @param event_col Character. Event ID column. Default "event_id".
#' @param class_col Character. Antibiotic class column. Default "antibiotic_class".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param rr_col Character. RR value column. Default "rr_value".
#'   If missing, only hierarchy is used.
#' @param hierarchy Named numeric vector. Custom hierarchy (class name -> rank).
#'   If NULL, uses default from get_beta_lactam_hierarchy().
#' @param filter_resistant Logical. If TRUE, only consider resistant (R) classes.
#'   Default TRUE.
#'
#' @return Data frame filtered to one resistance class per event
#' @export
#'
#' @examples
#' \dontrun{
#' # Select single resistance class per event
#' selected <- select_resistance_class(data)
#'
#' # Include susceptible classes too
#' selected <- select_resistance_class(data, filter_resistant = FALSE)
#' }
select_resistance_class <- function(data,
                                    event_col = "event_id",
                                    class_col = "antibiotic_class",
                                    susceptibility_col = "antibiotic_value",
                                    rr_col = "rr_value",
                                    hierarchy = NULL,
                                    filter_resistant = TRUE) {
  # Validate columns
  required_cols <- c(event_col, class_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Get hierarchy
  if (is.null(hierarchy)) {
    hierarchy <- get_beta_lactam_hierarchy()
  }

  # Check if RR column exists
  use_rr <- rr_col %in% names(data)
  if (!use_rr) {
    message(sprintf(
      "RR column '%s' not found. Using hierarchy only for selection.",
      rr_col
    ))
  }

  n_before <- nrow(data)
  n_events_before <- dplyr::n_distinct(data[[event_col]])

  message(sprintf(
    "Selecting resistance classes using hierarchy%s...",
    ifelse(use_rr, " + RR", "")
  ))

  # Filter to resistant classes if requested
  if (filter_resistant) {
    data <- data %>%
      dplyr::filter(!!rlang::sym(susceptibility_col) == "R")

    message(sprintf(
      "Filtered to resistant classes: %d rows",
      nrow(data)
    ))
  }

  # Apply prioritization
  selected <- prioritize_resistance(
    data = data,
    event_col = event_col,
    class_col = class_col,
    rr_col = if (use_rr) rr_col else NULL,
    hierarchy = hierarchy
  )

  n_after <- nrow(selected)
  n_events_after <- dplyr::n_distinct(selected[[event_col]])

  message(sprintf(
    "Selected: %d rows from %d events (avg %.2f classes/event before -> 1.0 after)",
    n_after,
    n_events_after,
    n_before / n_events_before
  ))

  # Show selection summary
  selection_summary <- selected %>%
    dplyr::count(!!rlang::sym(class_col), name = "n_events") %>%
    dplyr::arrange(dplyr::desc(n_events)) %>%
    utils::head(10)

  message("\nTop 10 selected classes:")
  print(selection_summary)

  return(selected)
}


#' Prioritize Resistance
#'
#' Helper function that applies hierarchy + RR ranking to select
#' the most important resistance class per event.
#'
#' @param data Data frame
#' @param event_col Character. Event ID column.
#' @param class_col Character. Class column.
#' @param rr_col Character or NULL. RR column. If NULL, uses hierarchy only.
#' @param hierarchy Named numeric vector. Hierarchy mapping.
#'
#' @return Data frame with one row per event (highest priority class)
#' @export
#'
#' @keywords internal
prioritize_resistance <- function(data,
                                  event_col,
                                  class_col,
                                  rr_col = NULL,
                                  hierarchy) {
  # Map classes to hierarchy ranks
  hierarchy_df <- data.frame(
    class = names(hierarchy),
    hierarchy_rank = as.numeric(hierarchy),
    stringsAsFactors = FALSE
  )
  names(hierarchy_df)[1] <- class_col

  data <- data %>%
    dplyr::left_join(hierarchy_df, by = class_col)

  # Assign default rank for unmapped classes (lowest priority)
  max_rank <- max(hierarchy_df$hierarchy_rank, na.rm = TRUE)
  data <- data %>%
    dplyr::mutate(
      hierarchy_rank = dplyr::coalesce(hierarchy_rank, max_rank + 1)
    )

  # Priority logic
  if (!is.null(rr_col) && rr_col %in% names(data)) {
    # Priority: hierarchy rank (lower = better) -> RR (higher = better) -> alphabetical
    selected <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(
        hierarchy_rank, # Lower rank = higher priority
        dplyr::desc(!!rlang::sym(rr_col)), # Higher RR = higher priority
        !!rlang::sym(class_col) # Alphabetical tie-breaker
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        selection_method = "hierarchy_rr",
        selection_confidence = "high"
      )
  } else {
    # Priority: hierarchy rank only -> alphabetical
    selected <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(
        hierarchy_rank,
        !!rlang::sym(class_col)
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        selection_method = "hierarchy_only",
        selection_confidence = "medium"
      )
  }

  # Clean up temporary column
  selected <- selected %>%
    dplyr::select(-hierarchy_rank)

  # Count how many classes were dropped per event
  classes_per_event <- data %>%
    dplyr::group_by(!!rlang::sym(event_col)) %>%
    dplyr::summarise(n_classes = dplyr::n(), .groups = "drop")

  multi_class_events <- classes_per_event %>%
    dplyr::filter(n_classes > 1)

  if (nrow(multi_class_events) > 0) {
    message(sprintf(
      "Applied selection to %d events with multiple resistant classes",
      nrow(multi_class_events)
    ))

    # Show example of selection
    example_event <- multi_class_events[[event_col]][1]
    example_before <- data %>%
      dplyr::filter(!!rlang::sym(event_col) == example_event) %>%
      dplyr::select(
        !!rlang::sym(event_col),
        !!rlang::sym(class_col),
        dplyr::any_of(c(rr_col, "hierarchy_rank"))
      )

    example_after <- selected %>%
      dplyr::filter(!!rlang::sym(event_col) == example_event) %>%
      dplyr::select(
        !!rlang::sym(event_col),
        !!rlang::sym(class_col),
        selection_method
      )

    message(sprintf("\nExample selection for event '%s':", example_event))
    message("Before (all resistant classes):")
    print(example_before)
    message("After (selected class):")
    print(example_after)
  }

  return(selected)
}
#' Flag Polymicrobial Infections (No Specimen Type, No Episode ID)
#'
#' Flags polymicrobial status by counting distinct organisms per patient
#' (optionally within facility/syndrome scope). No specimen type and no
#' episode_id are used.
#'
#' @param data          Data frame with patient and organism columns.
#' @param patient_col   Character. Patient ID column. Default \code{"patient_id"}.
#' @param organism_col  Character. Organism column.
#'   Default \code{"organism_normalized"}.
#' @param facility_col  Character or \code{NULL}. Facility column. When supplied,
#'   polymicrobial status is scoped within each facility.
#' @param facility_name Character or \code{NULL}. When supplied with
#'   \code{facility_col}, data are first filtered to that facility.
#' @param syndrome_col  Character or \code{NULL}. Syndrome column. When supplied,
#'   polymicrobial status is scoped within each syndrome.
#' @param syndrome_name Character or \code{NULL}. When supplied with
#'   \code{syndrome_col}, data are first filtered to that syndrome.
#'
#' @return Data frame with \code{n_organisms} and \code{is_polymicrobial} (0/1).
#' @export
flag_polymicrobial <- function(data,
                               patient_col = "patient_id",
                               organism_col = "organism_normalized",
                               facility_col = NULL,
                               facility_name = NULL,
                               syndrome_col = NULL,
                               syndrome_name = NULL) {
  # -- Validate columns --------------------------------------------------------
  required_cols <- c(patient_col, organism_col)
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be supplied when facility_name is specified.")
  }
  if (!is.null(syndrome_name) && is.null(syndrome_col)) {
    stop("syndrome_col must be supplied when syndrome_name is specified.")
  }

  # -- Optional filters --------------------------------------------------------
  if (!is.null(facility_name)) {
    data <- data %>% dplyr::filter(.data[[facility_col]] == facility_name)
    message(sprintf("Filtered to facility: %s (%d rows)", facility_name, nrow(data)))
  }
  if (!is.null(syndrome_name)) {
    data <- data %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
    message(sprintf("Filtered to syndrome: %s (%d rows)", syndrome_name, nrow(data)))
  }

  message("Identifying polymicrobial infections (no specimen type, no episode_id)...")

  # -- Temp columns ------------------------------------------------------------
  data <- data %>%
    dplyr::mutate(
      .temp_patient  = as.character(.data[[patient_col]]),
      .temp_organism = as.character(.data[[organism_col]])
    )

  if (!is.null(facility_col)) data$.temp_facility <- as.character(data[[facility_col]])
  if (!is.null(syndrome_col)) data$.temp_syndrome <- as.character(data[[syndrome_col]])

  data$.temp_patient[is.na(data$.temp_patient)] <- "__NA__"

  # Group scope: patient (+ optional facility/syndrome)
  group_cols <- c(
    ".temp_patient",
    if (!is.null(facility_col)) ".temp_facility",
    if (!is.null(syndrome_col)) ".temp_syndrome"
  )

  # Distinct patient-context x organism
  distinct_cols <- c(group_cols, ".temp_organism")
  df_unique <- data %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(distinct_cols)))

  # Count organisms per patient-context
  poly_counts <- df_unique %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_organisms      = dplyr::n_distinct(.temp_organism),
      is_polymicrobial = as.integer(n_organisms > 1),
      .groups          = "drop"
    )

  # Join back to full data
  data <- data %>%
    dplyr::left_join(poly_counts, by = group_cols)

  # Summary
  n_groups <- poly_counts %>% nrow()
  n_poly <- sum(poly_counts$is_polymicrobial == 1, na.rm = TRUE)
  pct_poly <- if (n_groups > 0) 100 * n_poly / n_groups else 0

  message(sprintf(
    "\nPolymicrobial: %d/%d patient-context groups (%.1f%%)",
    n_poly, n_groups, pct_poly
  ))
  message("\nOrganism count distribution per patient-context:")
  org_dist <- poly_counts %>%
    dplyr::count(n_organisms, name = "n_groups") %>%
    dplyr::arrange(n_organisms)
  print(org_dist)

  # Cleanup temp columns
  data <- data %>%
    dplyr::select(-dplyr::starts_with(".temp_"))

  return(data)
}


#' Compute Polymicrobial Weights
#'
#' Calculates proportional weights for polymicrobial infections.  Monomicrobial
#' patients always receive weight = 1.0.  For polymicrobial patients, weight is
#' computed per organism within each episode using one of three methods.
#'
#' When \code{facility_col} or \code{syndrome_col} are supplied, the
#' monomicrobial reference pool (method \code{"monomicrobial_proportion"}) is
#' computed within each stratum, so the reference distribution is local to each
#' facility / syndrome rather than global.
#'
#' @param data             Data frame with \code{episode_id},
#'   \code{is_polymicrobial} (0/1), and organism columns (output of
#'   \code{flag_polymicrobial()}).
#' @param episode_col      Character. Episode ID column.
#'   Default \code{"episode_id"}.
#' @param organism_col     Character. Organism column.
#'   Default \code{"organism_normalized"}.
#' @param polymicrobial_col Character. Polymicrobial flag column (0/1).
#'   Default \code{"is_polymicrobial"}.
#' @param method           Character. Weighting method:
#'   \code{"monomicrobial_proportion"} (default), \code{"equal"}, or
#'   \code{"manual"}.
#' @param weight_map       Named numeric vector. Custom organism weights when
#'   \code{method = "manual"}.
#' @param facility_col     Character or \code{NULL}. Facility column. When
#'   supplied, monomicrobial proportions are computed per-facility so each
#'   facility's local organism distribution is used as reference.
#' @param facility_name    Character or \code{NULL}. When supplied together
#'   with \code{facility_col}, data are first filtered to that facility.
#' @param syndrome_col     Character or \code{NULL}. Syndrome column. When
#'   supplied, monomicrobial proportions are computed per-syndrome.
#' @param syndrome_name    Character or \code{NULL}. When supplied together
#'   with \code{syndrome_col}, data are first filtered to that syndrome.
#'
#' @return Data frame with \code{polymicrobial_weight} column (range 0-1),
#'   plus \code{weight_method} and \code{weight_confidence} audit columns.
#'   \code{episode_id} is removed (internal use only).
#' @export
#'
#' @examples
#' \dontrun{
#' # Global monomicrobial proportion weights
#' data_weighted <- compute_polymicrobial_weight(data_flagged)
#'
#' # Per-facility reference distribution
#' data_weighted <- compute_polymicrobial_weight(
#'   data_flagged,
#'   facility_col = "center_name"
#' )
#'
#' # Filter to one facility + one syndrome, then weight
#' data_weighted <- compute_polymicrobial_weight(
#'   data_flagged,
#'   facility_col  = "center_name",
#'   facility_name = "PGIMER",
#'   syndrome_col  = "infectious_syndrome",
#'   syndrome_name = "Bloodstream infection"
#' )
#' }
compute_polymicrobial_weight <- function(data,
                                         episode_col = "episode_id",
                                         organism_col = "organism_normalized",
                                         polymicrobial_col = "is_polymicrobial",
                                         method = "monomicrobial_proportion",
                                         weight_map = NULL,
                                         facility_col = NULL,
                                         facility_name = NULL,
                                         syndrome_col = NULL,
                                         syndrome_name = NULL) {
  # -- Validate columns --------------------------------------------------------
  required_cols <- c(episode_col, organism_col, polymicrobial_col)
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  valid_methods <- c("monomicrobial_proportion", "equal", "manual")
  if (!method %in% valid_methods) {
    stop(sprintf("method must be one of: %s", paste(valid_methods, collapse = ", ")))
  }

  if (method == "manual" && is.null(weight_map)) {
    stop("weight_map must be provided when method = 'manual'.")
  }

  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be supplied when facility_name is specified.")
  }
  if (!is.null(syndrome_name) && is.null(syndrome_col)) {
    stop("syndrome_col must be supplied when syndrome_name is specified.")
  }

  # -- Optional filters --------------------------------------------------------
  if (!is.null(facility_name)) {
    data <- data %>% dplyr::filter(.data[[facility_col]] == facility_name)
    message(sprintf("Filtered to facility: %s (%d rows)", facility_name, nrow(data)))
  }
  if (!is.null(syndrome_name)) {
    data <- data %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
    message(sprintf("Filtered to syndrome: %s (%d rows)", syndrome_name, nrow(data)))
  }

  message(sprintf("Computing polymicrobial weights using method: %s", method))

  # Stratification columns for per-facility / per-syndrome reference
  strat_cols <- c(
    if (!is.null(facility_col)) facility_col,
    if (!is.null(syndrome_col)) syndrome_col
  )

  # -- Initialise weight column ------------------------------------------------
  data <- data %>% dplyr::mutate(polymicrobial_weight = 1.0)

  # -- Method A: Monomicrobial proportion -------------------------------------
  if (method == "monomicrobial_proportion") {
    # Compute monomicrobial proportions -- per stratum if strat_cols provided,
    # globally otherwise.
    mono_base <- data %>% dplyr::filter(.data[[polymicrobial_col]] == 0)

    if (length(strat_cols) > 0) {
      mono_proportions <- mono_base %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(strat_cols, organism_col)))) %>%
        dplyr::summarise(n_mono = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(strat_cols))) %>%
        dplyr::mutate(
          total_mono = sum(n_mono),
          proportion = n_mono / total_mono
        ) %>%
        dplyr::ungroup()
    } else {
      mono_proportions <- mono_base %>%
        dplyr::count(!!rlang::sym(organism_col), name = "n_mono") %>%
        dplyr::mutate(
          total_mono = sum(n_mono),
          proportion = n_mono / total_mono
        )
    }

    message(sprintf(
      "Calculated monomicrobial proportions for %d organism%s%s",
      nrow(mono_proportions),
      if (nrow(mono_proportions) != 1) "s" else "",
      if (length(strat_cols) > 0) {
        sprintf(" across %s strata", paste(strat_cols, collapse = " x "))
      } else {
        ""
      }
    ))

    # Join proportions back (by organism + strat_cols if present)
    join_by <- c(strat_cols, organism_col)
    data <- data %>%
      dplyr::left_join(
        mono_proportions %>%
          dplyr::select(dplyr::all_of(c(join_by, "proportion"))),
        by = join_by
      ) %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        polymicrobial_weight = dplyr::case_when(
          .data[[polymicrobial_col]] == 0 ~ 1.0,
          is.na(proportion) ~ 1.0 / dplyr::n(), # equal fallback
          TRUE ~ proportion / sum(proportion, na.rm = TRUE)
        ),
        weight_method = dplyr::case_when(
          .data[[polymicrobial_col]] == 0 ~ "monomicrobial",
          is.na(proportion) ~ "equal_fallback",
          TRUE ~ "monomicrobial_proportion"
        ),
        weight_confidence = dplyr::case_when(
          .data[[polymicrobial_col]] == 0 ~ "high",
          weight_method == "equal_fallback" ~ "low",
          TRUE ~ "high"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-proportion)

    n_fallback <- sum(data$weight_method == "equal_fallback", na.rm = TRUE)
    if (n_fallback > 0) {
      message(sprintf(
        "Warning: %d isolate(s) used equal weighting (no monomicrobial reference in stratum).",
        n_fallback
      ))
    }
  }

  # -- Method B: Equal weighting -----------------------------------------------
  else if (method == "equal") {
    data <- data %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        polymicrobial_weight = 1.0 / dplyr::n(),
        weight_method        = "equal",
        weight_confidence    = "medium"
      ) %>%
      dplyr::ungroup()

    message("Applied equal weighting to all organisms within episodes.")
  }

  # -- Method C: Manual weights ------------------------------------------------
  else if (method == "manual") {
    manual_weights_df <- data.frame(
      organism = names(weight_map),
      manual_weight = as.numeric(weight_map),
      stringsAsFactors = FALSE
    )
    names(manual_weights_df)[1L] <- organism_col

    data <- data %>%
      dplyr::left_join(manual_weights_df, by = organism_col) %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        polymicrobial_weight = dplyr::case_when(
          !is.na(manual_weight) ~ manual_weight / sum(manual_weight, na.rm = TRUE),
          TRUE ~ 1.0 / dplyr::n() # fallback for unmapped organisms
        ),
        weight_method = dplyr::case_when(
          !is.na(manual_weight) ~ "manual",
          TRUE ~ "equal_fallback"
        ),
        weight_confidence = dplyr::case_when(
          weight_method == "manual" ~ "user_defined",
          TRUE ~ "low"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-manual_weight)

    message(sprintf("Applied manual weights for %d organisms.", length(weight_map)))
  }

  # -- Validation: weights should sum to 1.0 per episode ----------------------
  weight_sums <- data %>%
    dplyr::group_by(!!rlang::sym(episode_col)) %>%
    dplyr::summarise(
      total_weight = sum(polymicrobial_weight, na.rm = TRUE),
      .groups      = "drop"
    )

  invalid_sums <- weight_sums %>%
    dplyr::filter(abs(total_weight - 1.0) > 0.01)

  if (nrow(invalid_sums) > 0) {
    warning(sprintf(
      "%d episode(s) have weights not summing to 1.0 (mixed-method fallback).",
      nrow(invalid_sums)
    ))
  }

  # -- Summary -----------------------------------------------------------------
  message("\nWeight statistics (polymicrobial isolates only):")
  weight_summary <- data %>%
    dplyr::filter(.data[[polymicrobial_col]] == 1) %>%
    dplyr::summarise(
      mean_weight   = mean(polymicrobial_weight, na.rm = TRUE),
      median_weight = median(polymicrobial_weight, na.rm = TRUE),
      min_weight    = min(polymicrobial_weight, na.rm = TRUE),
      max_weight    = max(polymicrobial_weight, na.rm = TRUE)
    )
  print(weight_summary)

  # Remove episode_id (internal use only; avoid confusion with event_id)
  if (episode_col %in% names(data)) {
    data <- data %>% dplyr::select(-dplyr::all_of(episode_col))
  }

  return(data)
}
