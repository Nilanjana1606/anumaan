# prep_clean_and_standardize.R
# Cleaning and standardization functions for AMR data preprocessing


# Store package root when this file is sourced (set at source time)
.amrburden_pkg_root <- NULL

# Try to detect package root at source time
local({
  # Get the path of this source file
  source_path <- tryCatch(
    {
      # Check various ways to get the source file path
      if (exists("ofile", envir = parent.frame(2))) {
        get("ofile", envir = parent.frame(2))
      } else {
        # Try to get from sys.frames
        for (i in seq_len(sys.nframe())) {
          env <- sys.frame(i)
          if (exists("ofile", envir = env)) {
            return(get("ofile", envir = env))
          }
        }
        NULL
      }
    },
    error = function(e) NULL
  )

  if (!is.null(source_path) && file.exists(source_path)) {
    # This file is in R/, so package root is parent of R/
    r_dir <- dirname(normalizePath(source_path))
    pkg_root <- dirname(r_dir)

    # Verify it's the anumaan package
    desc_path <- file.path(pkg_root, "DESCRIPTION")
    if (file.exists(desc_path)) {
      .amrburden_pkg_root <<- pkg_root
    }
  }
})


#' Find Package Data File
#'
#' Helper function to locate data files in inst/extdata.
#' Works both when the package is installed and during development.
#'
#' @param filename Character. Name of the file to find (e.g., "organisms.csv").
#' @return Character. Full path to the file, or empty string if not found.
#' @keywords internal
find_extdata_file <- function(filename) {
  # First try system.file (works when package is installed)
  file_path <- system.file("extdata", filename, package = "anumaan")

  if (file_path != "" && file.exists(file_path)) {
    return(file_path)
  }

  # Strategy 1: Use package root detected at source time
  if (!is.null(.amrburden_pkg_root)) {
    extdata_path <- file.path(.amrburden_pkg_root, "inst", "extdata", filename)
    if (file.exists(extdata_path)) {
      return(normalizePath(extdata_path))
    }
  }

  # Strategy 2: Check relative to current working directory
  # Cover common directory structures (including deep nesting like analysis/X/Y/)
  possible_roots <- c(
    ".", # Working dir is package root
    "..", # Working dir is R/ or subdirectory
    "../..", # Working dir is 2 levels deep
    "../../..", # Working dir is 3 levels deep (e.g., analysis/X/Y/)
    "../../../..", # Working dir is 4 levels deep
    "anumaan", # Working dir is parent of package
    "../anumaan", # Working dir is sibling of package
    "../../anumaan", # Working dir is nested sibling
    "../../../anumaan" # Working dir is deeply nested
  )

  for (root in possible_roots) {
    extdata_path <- file.path(root, "inst", "extdata", filename)
    if (file.exists(extdata_path)) {
      return(normalizePath(extdata_path))
    }
  }

  # File not found
  return("")
}


#' Standardize Column Names to Package Convention
#'
#' Maps incoming dataset column names to standardized names used throughout
#' the package. Supports exact matching and optional fuzzy matching for
#' unmatched columns.
#'
#' @param data A data frame with raw column names
#' @param mapping Named list where names are standard column names and values
#'   are character vectors of acceptable aliases. Default uses
#'   \code{default_column_mappings}.
#' @param fuzzy_match Logical. If TRUE, attempts fuzzy matching for unmapped
#'   columns using string distance. Default TRUE.
#' @param fuzzy_threshold Numeric. Maximum string distance (0-1) for fuzzy
#'   matching. Lower values require closer matches. Default 0.3.
#' @param interactive Logical. If TRUE and fuzzy matches found, prompts user
#'   for confirmation. Default FALSE (auto-accept).
#'
#' @return A list with components:
#'   \itemize{
#'     \item data: Data frame with standardized column names
#'     \item mapping_log: List documenting which columns were mapped and how
#'     \item unmapped: Character vector of columns that couldn't be mapped
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' raw_data <- data.frame(
#'   PatientID = 1:10,
#'   Organism = rep("E. coli", 10),
#'   Drug = rep("Ampicillin", 10)
#' )
#' result <- prep_standardize_column_names(raw_data)
#' clean_data <- result$data
#' }
prep_standardize_column_names <- function(data,
                                     mapping = default_column_mappings,
                                     fuzzy_match = TRUE,
                                     fuzzy_threshold = 0.3,
                                     interactive = FALSE) {
  original_names <- names(data)
  new_names <- names(data)
  mapping_log <- list()

  # Phase 1: Exact matching
  for (std_name in names(mapping)) {
    aliases <- mapping[[std_name]]
    matches <- which(original_names %in% aliases)

    if (length(matches) > 1) {
      warning(sprintf(
        "Multiple columns match '%s': %s. Using first match.",
        std_name,
        paste(original_names[matches], collapse = ", ")
      ))
      matches <- matches[1]
    }

    if (length(matches) == 1) {
      new_names[matches] <- std_name
      mapping_log[[std_name]] <- list(
        original = original_names[matches],
        method = "exact_match"
      )
    }
  }

  # Phase 2: Fuzzy matching (optional)
  if (fuzzy_match) {
    unmapped_idx <- which(new_names == original_names)
    unmapped <- original_names[unmapped_idx]

    for (std_name in names(mapping)) {
      if (std_name %in% new_names) next # Already mapped

      # Calculate string distances
      distances <- stringdist::stringdist(
        tolower(std_name),
        tolower(unmapped),
        method = "jw" # Jaro-Winkler distance
      )

      best_match_idx <- which.min(distances)

      if (length(best_match_idx) > 0 && distances[best_match_idx] < fuzzy_threshold) {
        if (interactive) {
          message(sprintf(
            "Fuzzy match: '%s' -> '%s' (distance: %.2f)",
            unmapped[best_match_idx],
            std_name,
            distances[best_match_idx]
          ))
          confirm <- readline(prompt = "Accept this mapping? (y/n): ")
          accept <- tolower(confirm) == "y"
        } else {
          accept <- TRUE
          message(sprintf(
            "Auto-accepted fuzzy match: '%s' -> '%s' (distance: %.2f)",
            unmapped[best_match_idx],
            std_name,
            distances[best_match_idx]
          ))
        }

        if (accept) {
          new_names[unmapped_idx[best_match_idx]] <- std_name
          mapping_log[[std_name]] <- list(
            original = unmapped[best_match_idx],
            method = "fuzzy_match",
            distance = distances[best_match_idx]
          )
        }
      }
    }
  }

  # Apply new names
  names(data) <- new_names

  # Report unmapped columns
  unmapped_final <- setdiff(new_names, names(mapping))
  if (length(unmapped_final) > 0) {
    message(
      "Columns not mapped to standard names: ",
      paste(unmapped_final, collapse = ", ")
    )
  }

  return(list(
    data = data,
    mapping_log = mapping_log,
    unmapped = unmapped_final
  ))
}


#' Normalize Organism Names
#'
#' Normalizes organism names using organisms.csv reference file.
#' Automatically handles abbreviations (E. coli), case variations, and typos.
#'
#' @param data Data frame with organism column
#' @param organism_col Character. Organism column name. Default "organism_name".
#' @param add_organism_group Logical. Add organism_group column from CSV. Default TRUE.
#' @param add_resistance_flags Logical. Add resistance flag columns (0/1). Default TRUE.
#'   Creates: is_MRSA, is_MSSA, is_MRCONS, is_MSCONS
#'
#' @return Data frame with organism_normalized, organism_group, and optionally resistance flag columns
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(organism = c("E. coli", "S. aureus", "MRSA"))
#' result <- prep_standardize_organisms(data, organism_col = "organism")
#' }
prep_standardize_organisms <- function(data,
                               organism_col = "organism_name",
                               add_organism_group = TRUE,
                               add_resistance_flags = TRUE) {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", organism_col))
  }

  # Load organisms.csv reference
  csv_path <- find_extdata_file("organisms.csv")

  if (csv_path == "" || !file.exists(csv_path)) {
    warning("organisms.csv not found in inst/extdata/. Organism normalization skipped.")
    data$organism_normalized <- data[[organism_col]]
    if (add_organism_group) data$organism_group <- NA_character_
    return(data)
  }

  # Load reference organisms
  org_ref <- readr::read_csv(csv_path, show_col_types = FALSE)

  # Helper function: extract keywords from organism name
  extract_keywords <- function(name) {
    if (is.na(name) || name == "") {
      return(character(0))
    }

    # Remove special characters, lowercase, split by spaces
    cleaned <- tolower(gsub("[^a-z0-9\\s]", " ", name))
    words <- unlist(strsplit(cleaned, "\\s+"))

    # Remove common filler words
    filler <- c(
      "spp", "sp", "species", "not", "specified", "other", "resistant",
      "positive", "negative", "gen", "generation"
    )
    words <- words[!words %in% filler & nchar(words) > 1]

    return(unique(words))
  }

  # Create lookup: extract keywords from reference organisms
  org_ref$keywords <- lapply(org_ref$organism_name, extract_keywords)
  org_ref$ref_lower <- tolower(trimws(org_ref$organism_name))

  # Process each unique input organism
  data$temp_org_input <- trimws(as.character(data[[organism_col]]))

  # ============================================================================
  # STEP 1: Detect and flag methicillin resistance patterns (MRSA, MRCONS, etc.)
  # ============================================================================
  if (add_resistance_flags) {
    data$is_MRSA <- 0L
    data$is_MSSA <- 0L
    data$is_MRCONS <- 0L
    data$is_MSCONS <- 0L

    mrsa_pattern <- "\\bmrsa\\b|\\(mrsa\\)|methicillin[- ]*resist[a-z]*[- ]*s[a-z]*[- ]*aureus|s\\.?\\s*aureus.*mrsa|staphylococcus\\s+aureus.*mrsa"
    data$is_MRSA[grepl(mrsa_pattern, data$temp_org_input, ignore.case = TRUE)] <- 1L

    mssa_pattern <- "\\bmssa\\b|\\(mssa\\)|methicillin[- ]*sens[a-z]*[- ]*s[a-z]*[- ]*aureus|methicillin[- ]*suscept[a-z]*[- ]*s[a-z]*[- ]*aureus|s\\.?\\s*aureus.*mssa|staphylococcus\\s+aureus.*mssa"
    data$is_MSSA[grepl(mssa_pattern, data$temp_org_input, ignore.case = TRUE)] <- 1L

    mrcons_pattern <- "\\bmr[- ]?cons\\b|\\bmrcos\\b|\\(mr[- ]?cons\\)|methicillin[- ]*resist[a-z]*[- ]*coagulase[- ]*neg"
    data$is_MRCONS[grepl(mrcons_pattern, data$temp_org_input, ignore.case = TRUE)] <- 1L

    mscons_pattern <- "\\bms[- ]?cons\\b|\\(ms[- ]?cons\\)|methicillin[- ]*sens[a-z]*[- ]*coagulase[- ]*neg|methicillin[- ]*suscept[a-z]*[- ]*coagulase[- ]*neg"
    data$is_MSCONS[grepl(mscons_pattern, data$temp_org_input, ignore.case = TRUE)] <- 1L

    data$temp_org_input[data$is_MRSA == 1 | data$is_MSSA == 1] <- "Staphylococcus aureus"

    cons_pattern <- "\\bcons\\b|\\bcns\\b|\\(cons\\)|\\(cns\\)|coagulase[- ]*neg|\\bco[a]?g[- ]*neg"
    data$temp_org_input[data$is_MRCONS == 1 | data$is_MSCONS == 1 | grepl(cons_pattern, data$temp_org_input, ignore.case = TRUE)] <- "Coagulase-negative staphylococci"
  } else {
    cons_pattern <- "\\bcons\\b|\\bcns\\b|\\(cons\\)|\\(cns\\)|coagulase[- ]*neg|\\bco[a]?g[- ]*neg|\\bmr[- ]?cons\\b|\\bmrcos\\b|\\bms[- ]?cons\\b"
    data$temp_org_input[grepl(cons_pattern, data$temp_org_input, ignore.case = TRUE)] <- "Coagulase-negative staphylococci"

    mrsa_mssa_pattern <- "\\bmrsa\\b|\\(mrsa\\)|\\bmssa\\b|\\(mssa\\)"
    data$temp_org_input[grepl(mrsa_mssa_pattern, data$temp_org_input, ignore.case = TRUE)] <- "Staphylococcus aureus"
  }

  # ============================================================================
  # STEP 2: Standardize species abbreviations
  # ============================================================================
  data$temp_org_input <- gsub("\\bsp\\.\\s*$", "spp.", data$temp_org_input, ignore.case = TRUE)
  data$temp_org_input <- gsub("\\bsp\\s+", "spp. ", data$temp_org_input, ignore.case = TRUE)
  data$temp_org_input <- gsub("\\bsp\\b", "spp.", data$temp_org_input, ignore.case = TRUE)
  data$temp_org_input <- gsub("non[- ]?ferm[ea]nt[ia]ng", "Non-fermenting", data$temp_org_input, ignore.case = TRUE)

  data$organism_normalized <- NA_character_

  unique_inputs <- unique(data$temp_org_input[!is.na(data$temp_org_input) & data$temp_org_input != ""])
  organism_map <- setNames(rep(NA_character_, length(unique_inputs)), unique_inputs)

  for (input_org in unique_inputs) {
    input_keywords <- extract_keywords(input_org)

    if (length(input_keywords) == 0) {
      organism_map[input_org] <- tolower(input_org)
      next
    }

    scores <- sapply(1:nrow(org_ref), function(i) {
      ref_keywords <- org_ref$keywords[[i]]
      if (length(ref_keywords) == 0) return(0)
      overlap <- sum(input_keywords %in% ref_keywords)
      if (tolower(input_org) == org_ref$ref_lower[i]) return(1000)
      if (grepl(org_ref$ref_lower[i], tolower(input_org), fixed = TRUE) ||
        grepl(tolower(input_org), org_ref$ref_lower[i], fixed = TRUE)) {
        overlap <- overlap + 2
      }
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) return(0)
      return(overlap / total_unique)
    })

    max_score <- max(scores)
    if (max_score > 0.3) {
      top_matches <- which(scores == max_score)
      if (length(top_matches) == 1) {
        organism_map[input_org] <- org_ref$ref_lower[top_matches[1]]
      } else {
        input_lower <- tolower(input_org)
        tie_distances <- sapply(top_matches, function(idx) {
          adist(input_lower, org_ref$ref_lower[idx], ignore.case = TRUE)[1, 1]
        })
        best_tie_idx <- top_matches[which.min(tie_distances)]
        organism_map[input_org] <- org_ref$ref_lower[best_tie_idx]
      }
    } else {
      input_lower <- tolower(input_org)
      input_words <- strsplit(input_lower, "\\s+")[[1]]
      input_genus <- input_words[1]
      if (input_genus %in% c("non", "multi", "methicillin")) {
        input_genus <- if (length(input_words) > 1) input_words[2] else input_words[1]
      }
      genus_info <- lapply(org_ref$ref_lower, function(ref) {
        ref_words <- strsplit(ref, "\\s+")[[1]]
        ref_genus <- ref_words[1]
        if (ref_genus %in% c("non", "multi", "methicillin")) {
          ref_genus <- if (length(ref_words) > 1) ref_words[2] else ref_words[1]
        }
        genus_dist <- adist(input_genus, ref_genus, ignore.case = TRUE)[1, 1]
        list(ref_genus = ref_genus, genus_dist = genus_dist)
      })
      genus_distances <- sapply(genus_info, function(x) x$genus_dist)
      genus_threshold <- if (nchar(input_genus) >= 4) 2 else 1
      genus_matches <- which(genus_distances <= genus_threshold)
      if (length(genus_matches) > 0) {
        best_match_idx <- genus_matches[1]
        best_match_dist <- Inf
        for (idx in genus_matches) {
          full_dist <- adist(input_lower, org_ref$ref_lower[idx], ignore.case = TRUE)[1, 1]
          if (full_dist < best_match_dist) {
            best_match_dist <- full_dist
            best_match_idx <- idx
          }
        }
        organism_map[input_org] <- org_ref$ref_lower[best_match_idx]
      } else {
        organism_map[input_org] <- tolower(input_org)
      }
    }
  }

  data$organism_normalized <- organism_map[data$temp_org_input]
  data$organism_normalized[is.na(data$temp_org_input) | data$temp_org_input == ""] <- NA_character_

  if (add_organism_group) {
    data <- data %>%
      dplyr::left_join(
        org_ref %>% dplyr::select(ref_lower, organism_group),
        by = c("organism_normalized" = "ref_lower")
      )
  }

  if (add_organism_group && "organism_group" %in% names(data)) {
    norm_lower <- tolower(data$organism_normalized)
    data$organism_group[!is.na(norm_lower) & grepl("chryseobacterium", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("vibrio", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("ralstonia", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("delftia", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("elizabethkingia", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("ochrobacterium", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("shewanella", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("non fermenter", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("non lactose fermenting", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("pantoea", norm_lower)] <- "Enterobacterales"
    data$organism_group[!is.na(norm_lower) & grepl("kluyvera", norm_lower)] <- "Enterobacterales"
    data$organism_group[!is.na(norm_lower) & grepl("leclercia", norm_lower)] <- "Enterobacterales"
  }

  data$temp_org_input <- NULL

  n_matched <- sum(!is.na(data$organism_normalized))
  message(sprintf(
    "Normalized %d/%d organisms (%.1f%%)",
    n_matched, nrow(data), 100 * n_matched / nrow(data)
  ))

  n_unique <- dplyr::n_distinct(data$organism_normalized, na.rm = TRUE)
  message(sprintf("Result: %d unique organisms", n_unique))

  if (add_resistance_flags) {
    n_mrsa <- sum(data$is_MRSA == 1, na.rm = TRUE)
    n_mssa <- sum(data$is_MSSA == 1, na.rm = TRUE)
    n_mrcons <- sum(data$is_MRCONS == 1, na.rm = TRUE)
    n_mscons <- sum(data$is_MSCONS == 1, na.rm = TRUE)
    if (n_mrsa > 0 || n_mssa > 0 || n_mrcons > 0 || n_mscons > 0) {
      message(sprintf(
        "Resistance flags: MRSA=%d, MSSA=%d, MRCONS=%d, MSCONS=%d",
        n_mrsa, n_mssa, n_mrcons, n_mscons
      ))
    }
  }

  if (add_organism_group) {
    n_with_group <- sum(!is.na(data$organism_group))
    n_without_group <- sum(is.na(data$organism_group) & !is.na(data$organism_normalized))
    message(sprintf(
      "With organism_group: %d (%.1f%%)",
      n_with_group, 100 * n_with_group / nrow(data)
    ))
    if (n_without_group > 0) {
      message(sprintf("Warning: %d organisms matched but missing organism_group", n_without_group))
    }
  }

  return(data)
}


#' Parse Dates Safely
#'
#' Attempts to parse date columns using multiple common formats via lubridate.
#' Returns Date objects or NA for unparseable values.
#'
#' @param data Data frame containing date columns
#' @param date_columns Character vector of column names to parse as dates.
#'   Default includes common date fields.
#'
#' @return Data frame with date columns converted to Date class
#'
#' @export
prep_parse_dates <- function(data, date_columns = c(
                          "date_of_admission",
                          "date_of_culture",
                          "date_of_final_outcome",
                          "DOB"
                        )) {
  for (col in date_columns) {
    if (!col %in% names(data)) {
      message(sprintf("Date column '%s' not found, skipping", col))
      next
    }

    if (inherits(data[[col]], "Date")) {
      next
    }

    original_na <- sum(is.na(data[[col]]))

    data[[col]] <- suppressWarnings(
      lubridate::parse_date_time(
        data[[col]],
        orders = c(
          "Ymd", "ymd", "dmy", "mdy", "Y-m-d", "d-m-Y", "d/m/Y",
          "Ymd HMS", "ymd HMS", "dmy HMS", "mdy HMS"
        ),
        tz = "UTC"
      )
    ) %>% lubridate::as_date()

    new_na <- sum(is.na(data[[col]]))
    parsed_na <- new_na - original_na

    if (parsed_na > 0) {
      warning(sprintf(
        "Column '%s': %d values could not be parsed as dates",
        col, parsed_na
      ))
    }

    message(sprintf(
      "Parsed '%s': %d dates (%d NA)",
      col, sum(!is.na(data[[col]])), new_na
    ))
  }

  return(data)
}


#' Standardize Sex Values
#'
#' Maps various gender/sex representations to standard "M" or "F" values.
#'
#' @param data Data frame containing gender column
#' @param col Character. Name of gender column. Default "gender".
#'
#' @return Data frame with standardized gender values
#'
#' @export
prep_standardize_sex <- function(data, col = "gender") {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col))
  }

  original_values <- unique(data[[col]])

  data[[col]] <- dplyr::case_when(
    toupper(data[[col]]) %in% c("M", "MALE", "MAN") ~ "M",
    toupper(data[[col]]) %in% c("F", "FEMALE", "WOMAN") ~ "F",
    TRUE ~ NA_character_
  )

  na_count <- sum(is.na(data[[col]]))
  if (na_count > 0) {
    warning(sprintf(
      "Gender column: %d values could not be standardized to M/F",
      na_count
    ))
  }

  message(sprintf(
    "Standardized gender: M=%d, F=%d, NA=%d",
    sum(data[[col]] == "M", na.rm = TRUE),
    sum(data[[col]] == "F", na.rm = TRUE),
    na_count
  ))

  return(data)
}


#' Standardize Outcome Values
#'
#' Maps various outcome representations to standard values:
#' "Died", "Discharged", "LAMA", "Unknown"
#'
#' @param data Data frame containing outcome column
#' @param col Character. Name of outcome column. Default "final_outcome".
#'
#' @return Data frame with standardized outcome values
#'
#' @export
prep_standardize_outcome <- function(data, col = "final_outcome") {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col))
  }

  data[[col]] <- dplyr::case_when(
    toupper(data[[col]]) %in% c("DIED", "DEATH", "DECEASED", "EXPIRED") ~ "Died",
    toupper(data[[col]]) %in% c("DISCHARGED", "DISCHARGE", "ALIVE", "SURVIVED") ~ "Discharged",
    toupper(data[[col]]) %in% c("LAMA", "DAMA", "LEFT", "ABSCONDED") ~ "LAMA",
    TRUE ~ "Unknown"
  )

  outcome_counts <- table(data[[col]])
  message("Outcome distribution:")
  print(outcome_counts)

  return(data)
}


#' Normalize Specimen/Sample Type
#'
#' Normalizes specimen/sample type names and adds sample_category and
#' sterile_classification from the reference CSV file.
#'
#' @param data Data frame with specimen column
#' @param specimen_col Character. Specimen column name. Default "specimen_type".
#' @param add_categories Logical. Add sample_category and sterile_classification.
#'   Default TRUE.
#'
#' @return Data frame with specimen_normalized, sample_category, and
#'   sterile_classification columns
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(specimen = c("Blood", "Urine", "CSF"))
#' result <- prep_standardize_specimens(data, specimen_col = "specimen")
#' }
prep_standardize_specimens <- function(data,
                               specimen_col = "specimen_type",
                               add_categories = TRUE) {
  if (!specimen_col %in% names(data)) {
    warning(sprintf("Column '%s' not found. Skipping specimen normalization.", specimen_col))
    return(data)
  }

  csv_path <- find_extdata_file("sample_type_classification.csv")

  if (csv_path == "" || !file.exists(csv_path)) {
    warning("sample_type_classification.csv not found in inst/extdata/. Specimen normalization skipped.")
    data$specimen_normalized <- data[[specimen_col]]
    if (add_categories) {
      data$sample_category <- NA_character_
      data$sterile_classification <- NA_character_
    }
    return(data)
  }

  specimen_ref <- readr::read_csv(csv_path, show_col_types = FALSE)

  data$temp_spec_input <- tolower(trimws(data[[specimen_col]]))
  data$temp_spec_clean <- data$temp_spec_input

  data$temp_spec_clean <- gsub("\\s*culture\\s*/\\s*sensitivity.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*culture.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*c/s.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*sensitivity.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*\\(.*?\\)\\s*", " ", data$temp_spec_clean)
  data$temp_spec_clean <- gsub("\\s*\\[.*?\\]\\s*", " ", data$temp_spec_clean)
  data$temp_spec_clean <- gsub("\\s*from\\s+.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*(catheter|foley|clean catch|voided).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*-\\s*(catheter|foley|midstream|peripheral|central).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("^(midstream|peripheral|central|clean|voided|fresh)\\s+", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*-\\s*(specimen|sample|report|test|site).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*:\\s*(specimen|sample|report|test|site).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\b(specimen|sample|site|swab from|aspirate from|tip|fluid from)\\b", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s+", " ", data$temp_spec_clean)
  data$temp_spec_clean <- trimws(data$temp_spec_clean)
  data$temp_spec_clean <- gsub("\\bet\\b", "endotracheal", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\beta\\b", "endotracheal aspirate", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\bbal\\s+fluid\\b", "bal", data$temp_spec_clean, ignore.case = TRUE)

  core_specimen_types <- c(
    "blood", "urine", "csf", "sputum", "stool", "pus", "bile",
    "faeces", "wound", "throat", "nasal", "vaginal", "cervical",
    "bal", "eta", "endotracheal"
  )

  for (spec_type in core_specimen_types) {
    if (any(grepl(paste0("\\b", spec_type, "\\b"), data$temp_spec_clean))) {
      idx <- grepl(paste0("\\b", spec_type, "\\b"), data$temp_spec_clean)
      if (any(idx)) {
        data$temp_spec_clean[idx] <- gsub(paste0(".*\\b(", spec_type, ")\\b.*"), "\\1", data$temp_spec_clean[idx])
      }
    }
  }

  data$temp_spec_clean <- trimws(data$temp_spec_clean)

  specimen_ref$temp_ref_lower <- tolower(trimws(specimen_ref$specimen_name))

  extract_spec_keywords <- function(name) {
    if (is.na(name) || name == "") return(character(0))
    cleaned <- tolower(gsub("[^a-z0-9\\s/]", " ", name))
    words <- unlist(strsplit(cleaned, "[\\s/]+"))
    filler <- c("a", "an", "the", "of", "from", "to", "and", "or", "is", "in", "on")
    words <- words[!words %in% filler & nchar(words) >= 2]
    return(unique(words))
  }

  specimen_ref$keywords <- lapply(specimen_ref$temp_ref_lower, extract_spec_keywords)

  data$specimen_normalized <- NA_character_
  data$sample_category <- NA_character_
  data$sterile_classification <- NA_character_

  unique_inputs <- unique(data$temp_spec_clean[!is.na(data$temp_spec_clean) & data$temp_spec_clean != ""])
  specimen_map <- setNames(rep(NA_integer_, length(unique_inputs)), unique_inputs)

  for (input_spec in unique_inputs) {
    exact_match <- which(specimen_ref$temp_ref_lower == input_spec)
    if (length(exact_match) > 0) {
      specimen_map[input_spec] <- exact_match[1]
      next
    }

    input_keywords <- extract_spec_keywords(input_spec)
    if (length(input_keywords) == 0) next

    scores <- sapply(1:nrow(specimen_ref), function(i) {
      ref_keywords <- specimen_ref$keywords[[i]]
      if (length(ref_keywords) == 0) return(0)
      overlap <- sum(input_keywords %in% ref_keywords)
      if (grepl(specimen_ref$temp_ref_lower[i], input_spec, fixed = TRUE) ||
        grepl(input_spec, specimen_ref$temp_ref_lower[i], fixed = TRUE)) {
        overlap <- overlap + 2
      }
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) return(0)
      return(overlap / total_unique)
    })

    max_score <- max(scores)
    if (max_score > 0.3) {
      top_matches <- which(scores == max_score)
      specimen_map[input_spec] <- top_matches[1]
    } else {
      distances <- sapply(specimen_ref$temp_ref_lower, function(ref) {
        adist(input_spec, ref, ignore.case = TRUE)[1, 1]
      })
      min_dist_idx <- which.min(distances)
      min_dist <- distances[min_dist_idx]
      if (min_dist <= 2) {
        specimen_map[input_spec] <- min_dist_idx
      }
    }
  }

  valid_mappings <- specimen_map[!is.na(specimen_map)]
  if (length(valid_mappings) > 0) {
    lookup_df <- data.frame(
      temp_spec_clean = names(valid_mappings),
      ref_idx = as.integer(valid_mappings),
      stringsAsFactors = FALSE
    )
    lookup_df$specimen_normalized <- specimen_ref$specimen_name[lookup_df$ref_idx]
    lookup_df$sample_category <- specimen_ref$sample_category[lookup_df$ref_idx]
    lookup_df$sterile_classification <- specimen_ref$sterile_classification[lookup_df$ref_idx]
    lookup_df$ref_idx <- NULL

    data$specimen_normalized <- NULL
    data$sample_category <- NULL
    data$sterile_classification <- NULL

    data <- dplyr::left_join(data, lookup_df, by = "temp_spec_clean")
  }

  data$temp_spec_input <- NULL
  data$temp_spec_clean <- NULL

  n_matched <- sum(!is.na(data$specimen_normalized))
  message(sprintf(
    "Specimen normalization: %d/%d matched (%.1f%%)",
    n_matched, nrow(data), 100 * n_matched / nrow(data)
  ))

  if (add_categories) {
    n_sterile <- sum(data$sterile_classification == "Sterile site", na.rm = TRUE)
    n_nonsterile <- sum(data$sterile_classification == "Non-sterile site", na.rm = TRUE)
    message(sprintf(
      "Sterile classification: %d sterile, %d non-sterile",
      n_sterile, n_nonsterile
    ))
  }

  return(data)
}


#' Normalize Antibiotic Names
#'
#' Standardizes antibiotic names using fuzzy matching against WHO reference.
#' Similar approach to organism normalization.
#'
#' @param data Data frame with antibiotic data
#' @param antibiotic_col Character. Column name with antibiotic names.
#'   Default "antibiotic_name".
#' @param who_table Data frame. WHO AWaRe classification table. If NULL,
#'   loads from inst/extdata/WHO_aware_class.csv.
#' @param add_class Logical. Add antibiotic_class column. Default TRUE.
#' @param add_aware Logical. Add aware_category column. Default TRUE.
#'
#' @return Data frame with antibiotic_normalized, antibiotic_class, aware_category
#' @export
#'
#' @examples
#' \dontrun{
#' data <- prep_standardize_antibiotics(data, antibiotic_col = "antibiotic_name")
#' }
prep_standardize_antibiotics <- function(data,
                                 antibiotic_col = "antibiotic_name",
                                 who_table = NULL,
                                 add_class = TRUE,
                                 add_aware = TRUE) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", antibiotic_col))
  }

  if (is.null(who_table)) {
    who_path <- find_extdata_file("WHO_aware_class.csv")

    if (who_path == "" || !file.exists(who_path)) {
      warning("WHO_aware_class.csv not found in inst/extdata/. Antibiotic normalization skipped.")
      data$antibiotic_normalized <- data[[antibiotic_col]]
      if (add_class) data$antibiotic_class <- NA_character_
      if (add_aware) data$aware_category <- NA_character_
      return(data)
    }

    who_table <- readr::read_csv(who_path, show_col_types = FALSE)
  }

  if ("Antibiotic" %in% names(who_table) && !"antibiotic_name" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(antibiotic_name = Antibiotic)
  }
  if ("Class" %in% names(who_table) && !"antibiotic_class" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(antibiotic_class = Class)
  }
  if ("Category " %in% names(who_table) && !"aware_category" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(aware_category = `Category `)
  } else if ("Category" %in% names(who_table) && !"aware_category" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(aware_category = Category)
  }

  if (!all(c("antibiotic_name", "antibiotic_class") %in% names(who_table))) {
    stop("WHO table must have 'antibiotic_name' and 'antibiotic_class' columns")
  }

  n_unique_before <- dplyr::n_distinct(data[[antibiotic_col]], na.rm = TRUE)

  message(sprintf(
    "Normalizing %d unique antibiotic names against WHO reference (%d antibiotics)...",
    n_unique_before, nrow(who_table)
  ))

  extract_abx_keywords <- function(name) {
    if (is.na(name) || name == "") return(character(0))
    cleaned <- tolower(gsub("[^a-z0-9\\s/-]", " ", name))
    words <- unlist(strsplit(cleaned, "[\\s/-]+"))
    filler <- c("iv", "oral", "injection", "tablet", "mg", "strip", "mic", "done", "by")
    words <- words[!words %in% filler & nchar(words) > 2]
    return(unique(words))
  }

  who_table$keywords <- lapply(who_table$antibiotic_name, extract_abx_keywords)
  who_table$ref_lower <- tolower(trimws(who_table$antibiotic_name))

  data$temp_abx_input <- trimws(as.character(data[[antibiotic_col]]))
  data$antibiotic_normalized <- NA_character_

  unique_inputs <- unique(data$temp_abx_input[!is.na(data$temp_abx_input) & data$temp_abx_input != ""])
  antibiotic_map <- setNames(rep(NA_character_, length(unique_inputs)), unique_inputs)

  for (input_abx in unique_inputs) {
    input_keywords <- extract_abx_keywords(input_abx)

    if (length(input_keywords) == 0) {
      antibiotic_map[input_abx] <- tolower(input_abx)
      next
    }

    scores <- sapply(1:nrow(who_table), function(i) {
      ref_keywords <- who_table$keywords[[i]]
      if (length(ref_keywords) == 0) return(0)
      overlap <- sum(input_keywords %in% ref_keywords)
      if (tolower(input_abx) == who_table$ref_lower[i]) return(1000)
      if (grepl(who_table$ref_lower[i], tolower(input_abx), fixed = TRUE) ||
        grepl(tolower(input_abx), who_table$ref_lower[i], fixed = TRUE)) {
        overlap <- overlap + 2
      }
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) return(0)
      return(overlap / total_unique)
    })

    max_score <- max(scores)
    if (max_score > 0.3) {
      top_matches <- which(scores == max_score)
      if (length(top_matches) == 1) {
        antibiotic_map[input_abx] <- who_table$ref_lower[top_matches[1]]
      } else {
        input_lower <- tolower(input_abx)
        tie_distances <- sapply(top_matches, function(idx) {
          adist(input_lower, who_table$ref_lower[idx], ignore.case = TRUE)[1, 1]
        })
        best_tie_idx <- top_matches[which.min(tie_distances)]
        antibiotic_map[input_abx] <- who_table$ref_lower[best_tie_idx]
      }
    } else {
      input_lower <- tolower(input_abx)
      distances <- sapply(who_table$ref_lower, function(ref) {
        adist(input_lower, ref, ignore.case = TRUE)[1, 1]
      })
      min_dist_idx <- which.min(distances)
      min_dist <- distances[min_dist_idx]
      if (min_dist <= 3) {
        antibiotic_map[input_abx] <- who_table$ref_lower[min_dist_idx]
      } else {
        antibiotic_map[input_abx] <- tolower(input_abx)
      }
    }
  }

  data$antibiotic_normalized <- antibiotic_map[data$temp_abx_input]
  data$antibiotic_normalized[is.na(data$temp_abx_input) | data$temp_abx_input == ""] <- NA_character_
  data$antibiotic_normalized[!is.na(data$antibiotic_normalized) &
    trimws(data$antibiotic_normalized) == ""] <- NA_character_

  if (add_class || add_aware) {
    join_cols <- c("antibiotic_normalized" = "ref_lower")
    select_cols <- "ref_lower"

    if (add_class) select_cols <- c(select_cols, "antibiotic_class")
    if (add_aware && "aware_category" %in% names(who_table)) {
      select_cols <- c(select_cols, "aware_category")
    }

    data <- data %>%
      dplyr::left_join(
        who_table %>% dplyr::select(dplyr::all_of(select_cols)),
        by = join_cols
      )
  }

  orig_empty <- is.na(data$temp_abx_input) | data$temp_abx_input == ""
  if (add_class && "antibiotic_class" %in% names(data)) {
    data$antibiotic_class[orig_empty] <- NA_character_
  }
  if (add_aware && "aware_category" %in% names(data)) {
    data$aware_category[orig_empty] <- NA_character_
  }

  if (add_class && "antibiotic_class" %in% names(data)) {
    orig_lower <- tolower(data$temp_abx_input)
    data$antibiotic_class[!is.na(orig_lower) & orig_lower == "colistin"] <- "Colistin"
    mupirocin_hl <- !is.na(orig_lower) &
      grepl("mupirocin", orig_lower) &
      grepl("high", orig_lower) &
      grepl("level", orig_lower)
    data$antibiotic_class[mupirocin_hl] <- "Mupirocin High level"
    data$antibiotic_class[!is.na(orig_lower) & grepl("nalidixic", orig_lower)] <- "Nalidixic acid"
    data$antibiotic_class[!is.na(orig_lower) & grepl("polymyxin.*b", orig_lower)] <- "Polymyxins"
  }

  data$temp_abx_input <- NULL

  n_unique_after <- dplyr::n_distinct(data$antibiotic_normalized, na.rm = TRUE)
  n_matched <- sum(!is.na(data$antibiotic_normalized))

  message(sprintf(
    "Normalized: %d unique names -> %d",
    n_unique_before, n_unique_after
  ))

  if (add_class) {
    n_with_class <- sum(!is.na(data$antibiotic_class))
    message(sprintf(
      "With antibiotic_class: %d (%.1f%%)",
      n_with_class, 100 * n_with_class / nrow(data)
    ))
  }

  if (add_aware && "aware_category" %in% names(data)) {
    n_with_aware <- sum(!is.na(data$aware_category))
    message(sprintf(
      "With AWaRe category: %d (%.1f%%)",
      n_with_aware, 100 * n_with_aware / nrow(data)
    ))
  }

  return(data)
}


#' Derive Age from Date of Birth
#'
#' Calculates age in years from date of birth and a reference date.
#' Implements proxy logic when Age is already present.
#'
#' @param data Data frame
#' @param dob_col Character. Name of DOB column. Default "DOB".
#' @param reference_date_col Character. Reference date column (usually culture date).
#'   Default "date_of_culture".
#' @param force Logical. If TRUE, recalculates even if Age present. Default FALSE.
#'
#' @return Data frame with Age column added/updated and Age_derived flag
#' @export
prep_derive_age <- function(data, dob_col = "DOB",
                       reference_date_col = "date_of_culture",
                       force = FALSE) {
  if ("Age" %in% names(data) && !force) {
    n_present <- sum(!is.na(data$Age))
    if (n_present > 0) {
      message(sprintf(
        "Age already present for %d records. Use force=TRUE to recalculate.",
        n_present
      ))
      data$Age_derived <- FALSE
      return(data)
    }
  }

  if (!dob_col %in% names(data)) {
    stop(sprintf("DOB column '%s' not found", dob_col))
  }

  if (!reference_date_col %in% names(data)) {
    stop(sprintf("Reference date column '%s' not found", reference_date_col))
  }

  data$Age <- as.numeric(
    difftime(data[[reference_date_col]], data[[dob_col]], units = "days")
  ) / 365.25

  data$Age <- ifelse(data$Age < 0, 0, data$Age)

  data$Age_derived <- TRUE
  data$Age_method <- "calculated_from_dob"

  n_derived <- sum(!is.na(data$Age))
  message(sprintf("Derived Age for %d records from DOB", n_derived))

  return(data)
}


#' Assign Age Bins
#'
#' Categorizes age into bins for stratification.
#'
#' @param data Data frame with Age column
#' @param age_col Character. Name of age column. Default "Age".
#' @param bins Character or numeric vector. Either "GBD_standard", "pediatric",
#'   "geriatric", or custom bin breaks. Default "GBD_standard".
#'
#' @return Data frame with Age_bin column (factor)
#' @export
prep_assign_age_bins <- function(data, age_col = "Age", bins = "GBD_standard") {
  if (!age_col %in% names(data)) {
    stop(sprintf("Age column '%s' not found", age_col))
  }

  if (is.character(bins) && length(bins) == 1) {
    bin_labels <- get_age_bins(bins)
  } else {
    bin_labels <- bins
  }

  breaks <- parse_age_bin_labels(bin_labels)

  data$Age_bin <- cut(
    data[[age_col]],
    breaks = breaks$breaks,
    labels = breaks$labels,
    right = FALSE,
    include.lowest = TRUE
  )

  n_binned <- sum(!is.na(data$Age_bin))
  n_unbinned <- sum(is.na(data$Age_bin) & !is.na(data[[age_col]]))

  message(sprintf(
    "Assigned age bins: %d binned, %d unbinned",
    n_binned, n_unbinned
  ))

  return(data)
}


#' Extract Genus from Organism Name
#'
#' Extracts bacterial genus (first word) from normalized organism name.
#'
#' @param data Data frame
#' @param organism_col Character. Organism column name.
#'   Default "organism_normalized".
#'
#' @return Data frame with org_genus column added
#' @export
prep_extract_genus <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Organism column '%s' not found", organism_col))
  }

  data$org_genus <- stringr::str_extract(
    data[[organism_col]],
    "^[a-z]+"
  )

  n_extracted <- sum(!is.na(data$org_genus))
  n_unique <- length(unique(data$org_genus[!is.na(data$org_genus)]))

  message(sprintf(
    "Extracted genus: %d records, %d unique genera",
    n_extracted, n_unique
  ))

  return(data)
}


#' Extract Species from Organism Name
#'
#' Extracts bacterial species (second word) from normalized organism name.
#'
#' @param data Data frame
#' @param organism_col Character. Organism column name.
#'   Default "organism_normalized".
#'
#' @return Data frame with org_species column added
#' @export
prep_extract_species <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Organism column '%s' not found", organism_col))
  }

  data$org_species <- stringr::str_extract(
    data[[organism_col]],
    "(?<=\\s)[a-z]+"
  )

  data$org_species <- ifelse(
    data$org_species == "spp",
    "species",
    data$org_species
  )

  n_extracted <- sum(!is.na(data$org_species))
  n_unique <- length(unique(data$org_species[!is.na(data$org_species)]))

  message(sprintf(
    "Extracted species: %d records, %d unique species",
    n_extracted, n_unique
  ))

  return(data)
}


#' Assign Organism Group
#'
#' Assigns organisms to taxonomic groups (Enterobacterales, Gram-positive, etc.)
#'
#' @param data Data frame
#' @param organism_col Character. Normalized organism column.
#'   Default "organism_normalized".
#'
#' @return Data frame with org_group column added
#' @export
prep_assign_organism_group <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Organism column '%s' not found", organism_col))
  }

  taxonomy <- get_organism_taxonomy()

  data <- data %>%
    dplyr::left_join(
      taxonomy,
      by = stats::setNames("organism_name", organism_col)
    )

  data$org_group <- ifelse(
    is.na(data$org_group),
    "Other",
    data$org_group
  )

  group_counts <- table(data$org_group)
  message("Organism group distribution:")
  print(group_counts)

  return(data)
}


#' Classify Antibiotic to WHO Class
#'
#' Maps antibiotic names to WHO antibiotic classes.
#'
#' @param data Data frame
#' @param antibiotic_col Character. Normalized antibiotic column.
#'   Default "antibiotic_normalized".
#' @param who_table Data frame. WHO classification table. If NULL, uses
#'   built-in mapping.
#'
#' @return Data frame with antibiotic_class column added
#' @export
prep_classify_antibiotic_class <- function(data,
                                      antibiotic_col = "antibiotic_normalized",
                                      who_table = NULL) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Antibiotic column '%s' not found", antibiotic_col))
  }

  if (is.null(who_table)) {
    message("Using built-in WHO class mapping")
    data$antibiotic_class <- NA_character_
    data$class_source <- "needs_who_table"
    warning("WHO table not provided. Please supply or package with data.")
    return(data)
  }

  data <- data %>%
    dplyr::left_join(
      who_table %>% dplyr::select(Antibiotic, Class),
      by = stats::setNames("Antibiotic", antibiotic_col)
    ) %>%
    dplyr::rename(antibiotic_class = Class)

  n_classified <- sum(!is.na(data$antibiotic_class))
  pct_classified <- 100 * n_classified / nrow(data)

  message(sprintf(
    "Classified antibiotics: %d/%d (%.1f%%)",
    n_classified, nrow(data), pct_classified
  ))

  return(data)
}


#' Classify AWaRe Category
#'
#' Assigns WHO AWaRe (Access, Watch, Reserve) categories to antibiotics.
#'
#' @param data Data frame
#' @param antibiotic_col Character. Antibiotic column.
#'   Default "antibiotic_normalized".
#' @param who_table Data frame. WHO AWaRe table. If NULL, uses built-in.
#'
#' @return Data frame with aware_category column added
#' @export
prep_classify_aware <- function(data,
                           antibiotic_col = "antibiotic_normalized",
                           who_table = NULL) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Antibiotic column '%s' not found", antibiotic_col))
  }

  if (is.null(who_table)) {
    message("Using built-in AWaRe mapping")
    data$aware_category <- NA_character_
    data$aware_source <- "needs_who_table"
    warning("WHO AWaRe table not provided.")
    return(data)
  }

  data <- data %>%
    dplyr::left_join(
      who_table %>% dplyr::select(Antibiotic, Category),
      by = stats::setNames("Antibiotic", antibiotic_col)
    ) %>%
    dplyr::rename(aware_category = Category)

  n_classified <- sum(!is.na(data$aware_category))

  message(sprintf("Classified AWaRe: %d records", n_classified))

  aware_dist <- table(data$aware_category, useNA = "ifany")
  print(aware_dist)

  return(data)
}


#' Calculate Length of Stay
#'
#' Calculates hospital length of stay in days.
#'
#' @param data Data frame
#' @param admission_col Character. Admission date column.
#' @param outcome_col Character. Outcome/discharge date column.
#'
#' @return Data frame with Length_of_stay column
#' @export
prep_calculate_los <- function(data,
                          admission_col = "date_of_admission",
                          outcome_col = "date_of_final_outcome") {
  if (!all(c(admission_col, outcome_col) %in% names(data))) {
    stop("Admission and outcome date columns not found")
  }

  data$Length_of_stay <- as.numeric(
    difftime(data[[outcome_col]], data[[admission_col]], units = "days")
  )

  data$LOS_suspicious <- data$Length_of_stay < 0 | data$Length_of_stay > 365

  n_suspicious <- sum(data$LOS_suspicious, na.rm = TRUE)
  if (n_suspicious > 0) {
    warning(sprintf(
      "%d records have suspicious LOS (<0 or >365 days)",
      n_suspicious
    ))
  }

  message(sprintf(
    "Calculated LOS: median=%.1f days, mean=%.1f days",
    median(data$Length_of_stay, na.rm = TRUE),
    mean(data$Length_of_stay, na.rm = TRUE)
  ))

  return(data)
}


#' Fill Missing Age Values
#'
#' Enriches Age data by deriving from DOB or retaining provided values,
#' tracking the method and confidence for each record.
#'
#' @param data Data frame
#' @param age_col Character. Age column name. Default "Age".
#' @param dob_col Character. Date of birth column. Default "DOB".
#' @param date_col Character. Reference date column. Default "date_of_culture".
#' @param overwrite Logical. Recalculate even if present. Default FALSE.
#'
#' @return Data frame with Age enriched
#' @export
#'
#' @examples
#' \dontrun{
#' data_enriched <- prep_fill_age(data)
#' }
prep_fill_age <- function(data,
                       age_col = "Age",
                       dob_col = "DOB",
                       date_col = "date_of_culture",
                       overwrite = FALSE) {
  has_age <- age_col %in% names(data)
  has_dob <- dob_col %in% names(data)
  has_date <- date_col %in% names(data)

  if (!has_age) {
    data[[age_col]] <- NA_real_
  }

  if (!"age_method" %in% names(data)) {
    data$age_method <- NA_character_
    data$age_confidence <- NA_character_
  }

  if (has_dob && has_date) {
    n_before_missing <- sum(is.na(data[[age_col]]))

    data <- data %>%
      dplyr::mutate(
        calculated_age = dplyr::case_when(
          !is.na(!!rlang::sym(dob_col)) & !is.na(!!rlang::sym(date_col)) ~
            as.numeric(difftime(!!rlang::sym(date_col), !!rlang::sym(dob_col), units = "days")) / 365.25,
          TRUE ~ NA_real_
        ),
        !!age_col := dplyr::case_when(
          overwrite & !is.na(calculated_age) ~ calculated_age,
          is.na(!!rlang::sym(age_col)) & !is.na(calculated_age) ~ calculated_age,
          TRUE ~ !!rlang::sym(age_col)
        ),
        age_method = dplyr::case_when(
          !is.na(calculated_age) & (overwrite | is.na(age_method)) ~ "calculated_from_dob",
          !is.na(age_method) ~ age_method,
          !is.na(!!rlang::sym(age_col)) ~ "provided",
          TRUE ~ NA_character_
        ),
        age_confidence = dplyr::case_when(
          age_method == "calculated_from_dob" ~ "high",
          age_method == "provided" ~ "high",
          TRUE ~ age_confidence
        )
      ) %>%
      dplyr::select(-calculated_age)

    n_after_missing <- sum(is.na(data[[age_col]]))
    n_enriched <- n_before_missing - n_after_missing

    if (n_enriched > 0) {
      message(sprintf(
        "Enriched Age: %d rows filled using DOB calculation",
        n_enriched
      ))
    }
  } else {
    data <- data %>%
      dplyr::mutate(
        age_method = dplyr::case_when(
          !is.na(!!rlang::sym(age_col)) & is.na(age_method) ~ "provided",
          TRUE ~ age_method
        ),
        age_confidence = dplyr::case_when(
          age_method == "provided" ~ "high",
          TRUE ~ age_confidence
        )
      )
  }

  age_summary <- data %>%
    dplyr::count(age_method, age_confidence) %>%
    dplyr::arrange(dplyr::desc(n))

  message("\nAge enrichment summary:")
  print(age_summary)

  n_still_missing <- sum(is.na(data[[age_col]]))
  if (n_still_missing > 0) {
    message(sprintf(
      "[!] Warning: %d rows still missing Age (%.1f%%)",
      n_still_missing,
      100 * n_still_missing / nrow(data)
    ))
  }

  return(data)
}


#' Fill Missing Length of Stay
#'
#' Derives Length of Stay (LOS) from admission and outcome dates.
#'
#' @param data Data frame
#' @param los_col Character. LOS column name. Default "Length_of_stay".
#' @param admission_col Character. Admission date. Default "date_of_admission".
#' @param outcome_col Character. Outcome date. Default "date_of_final_outcome".
#' @param overwrite Logical. Recalculate even if present. Default FALSE.
#'
#' @return Data frame with LOS enriched
#' @export
#'
#' @examples
#' \dontrun{
#' data_enriched <- prep_fill_los(data)
#' }
prep_fill_los <- function(data,
                       los_col = "Length_of_stay",
                       admission_col = "date_of_admission",
                       outcome_col = "date_of_final_outcome",
                       overwrite = FALSE) {
  has_los <- los_col %in% names(data)
  has_admission <- admission_col %in% names(data)
  has_outcome <- outcome_col %in% names(data)

  if (!has_los) {
    data[[los_col]] <- NA_real_
  }

  if (!has_admission || !has_outcome) {
    message(sprintf(
      "[!] Cannot calculate LOS: missing '%s' or '%s'",
      admission_col, outcome_col
    ))
    return(data)
  }

  n_before_missing <- sum(is.na(data[[los_col]]))

  data <- data %>%
    dplyr::mutate(
      calculated_los = dplyr::case_when(
        !is.na(!!rlang::sym(admission_col)) & !is.na(!!rlang::sym(outcome_col)) ~
          as.numeric(difftime(!!rlang::sym(outcome_col), !!rlang::sym(admission_col), units = "days")),
        TRUE ~ NA_real_
      ),
      !!los_col := dplyr::case_when(
        overwrite & !is.na(calculated_los) ~ calculated_los,
        is.na(!!rlang::sym(los_col)) & !is.na(calculated_los) ~ calculated_los,
        TRUE ~ !!rlang::sym(los_col)
      ),
      los_method = dplyr::case_when(
        !is.na(calculated_los) ~ "calculated_from_dates",
        !is.na(!!rlang::sym(los_col)) ~ "provided",
        TRUE ~ NA_character_
      ),
      los_confidence = dplyr::case_when(
        los_method %in% c("calculated_from_dates", "provided") ~ "high",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-calculated_los)

  negative_los <- data %>%
    dplyr::filter(!is.na(!!rlang::sym(los_col)) & !!rlang::sym(los_col) < 0)

  if (nrow(negative_los) > 0) {
    warning(sprintf(
      "[!] Warning: %d rows have negative LOS (date inconsistency)",
      nrow(negative_los)
    ))
  }

  n_after_missing <- sum(is.na(data[[los_col]]))
  n_enriched <- n_before_missing - n_after_missing

  if (n_enriched > 0) {
    message(sprintf(
      "Enriched LOS: %d rows filled",
      n_enriched
    ))
  }

  los_stats <- data %>%
    dplyr::filter(!is.na(!!rlang::sym(los_col))) %>%
    dplyr::summarise(
      mean_los = mean(!!rlang::sym(los_col), na.rm = TRUE),
      median_los = median(!!rlang::sym(los_col), na.rm = TRUE),
      min_los = min(!!rlang::sym(los_col), na.rm = TRUE),
      max_los = max(!!rlang::sym(los_col), na.rm = TRUE)
    )

  message("\nLOS summary:")
  print(los_stats)

  return(data)
}


#' Infer Hospital Department
#'
#' Infers hospital department from contextual data (specimen type,
#' diagnosis, patient demographics).
#'
#' @param data Data frame
#' @param department_col Character. Department column. Default "hospital_department".
#' @param specimen_col Character. Specimen type column. Default "specimen_type".
#' @param diagnosis_col Character. Diagnosis column. Default "diagnosis_1".
#' @param age_col Character. Age column. Default "Age".
#' @param overwrite Logical. Recalculate even if present. Default FALSE.
#'
#' @return Data frame with hospital_department enriched
#' @export
#'
#' @examples
#' \dontrun{
#' data_enriched <- prep_infer_department(data)
#' }
prep_infer_department <- function(data,
                                       department_col = "hospital_department",
                                       specimen_col = "specimen_type",
                                       diagnosis_col = "diagnosis_1",
                                       age_col = "Age",
                                       overwrite = FALSE) {
  has_department <- department_col %in% names(data)
  has_specimen <- specimen_col %in% names(data)
  has_diagnosis <- diagnosis_col %in% names(data)
  has_age <- age_col %in% names(data)

  if (!has_department) {
    data[[department_col]] <- NA_character_
  }

  if (!has_specimen && !has_diagnosis && !has_age) {
    message("[!] Cannot infer department: no contextual data available")
    return(data)
  }

  n_before_missing <- sum(is.na(data[[department_col]]))

  message("Inferring hospital department from contextual data...")

  if (!has_age) data[[age_col]] <- NA_real_
  if (!has_specimen) data[[specimen_col]] <- NA_character_
  if (!has_diagnosis) data[[diagnosis_col]] <- NA_character_

  data <- data %>%
    dplyr::mutate(
      inferred_dept = dplyr::case_when(
        !is.na(!!rlang::sym(age_col)) & !!rlang::sym(age_col) < 18 ~ "Pediatrics",
        !is.na(!!rlang::sym(specimen_col)) &
          grepl("blood|csf|broncho|balf", tolower(!!rlang::sym(specimen_col))) ~ "ICU",
        !is.na(!!rlang::sym(specimen_col)) &
          grepl("peritoneal|abscess|tissue|wound", tolower(!!rlang::sym(specimen_col))) ~ "Surgery",
        !is.na(!!rlang::sym(diagnosis_col)) &
          grepl("pregnancy|obstetric|maternal", tolower(!!rlang::sym(diagnosis_col))) ~ "Obstetrics",
        TRUE ~ NA_character_
      ),
      !!department_col := dplyr::case_when(
        overwrite & !is.na(inferred_dept) ~ inferred_dept,
        is.na(!!rlang::sym(department_col)) & !is.na(inferred_dept) ~ inferred_dept,
        TRUE ~ !!rlang::sym(department_col)
      ),
      department_method = dplyr::case_when(
        !is.na(inferred_dept) ~ "inferred_heuristic",
        !is.na(!!rlang::sym(department_col)) ~ "provided",
        TRUE ~ NA_character_
      ),
      department_confidence = dplyr::case_when(
        department_method == "provided" ~ "high",
        department_method == "inferred_heuristic" ~ "low",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-inferred_dept)

  n_after_missing <- sum(is.na(data[[department_col]]))
  n_enriched <- n_before_missing - n_after_missing

  if (n_enriched > 0) {
    message(sprintf(
      "Enriched department: %d rows filled (LOW confidence - heuristic)",
      n_enriched
    ))
  }

  dept_summary <- data %>%
    dplyr::count(!!rlang::sym(department_col), department_confidence) %>%
    dplyr::arrange(dplyr::desc(n))

  message("\nDepartment distribution:")
  print(dept_summary)

  return(data)
}


#' Clean Optional Columns
#'
#' Cleans and standardizes all optional columns in a single pass.
#' Handles value normalization, whitespace trimming, and consistency checks.
#'
#' @param data Data frame
#' @param optional_cols Character vector. Optional column names to groom.
#'   If NULL, grooms all known optional columns. Default NULL.
#'
#' @return Data frame with groomed optional columns
#' @export
#'
#' @examples
#' \dontrun{
#' data_groomed <- prep_clean_optional_columns(data)
#'
#' # Groom specific columns
#' data_groomed <- prep_clean_optional_columns(
#'   data,
#'   optional_cols = c("hospital_department", "unit_type", "comorbidities")
#' )
#' }
prep_clean_optional_columns <- function(data,
                                   optional_cols = NULL) {
  known_optional <- c(
    "infection_type", "device_inserted", "hospital_department",
    "aware_category", "unit_type", "comorbidities", "hospital_location",
    "previous_history", "specimen_category"
  )

  if (is.null(optional_cols)) {
    optional_cols <- intersect(known_optional, names(data))
  } else {
    optional_cols <- intersect(optional_cols, names(data))
  }

  if (length(optional_cols) == 0) {
    message("No optional columns to groom")
    return(data)
  }

  message(sprintf(
    "Grooming %d optional columns: %s",
    length(optional_cols),
    paste(optional_cols[1:min(3, length(optional_cols))], collapse = ", ")
  ))

  for (col in optional_cols) {
    if (is.character(data[[col]])) {
      data[[col]] <- trimws(data[[col]])
      data[[col]][data[[col]] == ""] <- NA_character_

      if (!col %in% c("comorbidities", "previous_history")) {
        data[[col]] <- tools::toTitleCase(tolower(data[[col]]))
      }

      if (col == "infection_type") {
        data[[col]] <- dplyr::case_when(
          grepl("^comm|^cai", tolower(data[[col]])) ~ "CAI",
          grepl("^hosp|^hai|^nosoco", tolower(data[[col]])) ~ "HAI",
          TRUE ~ data[[col]]
        )
      }

      if (col == "unit_type") {
        data[[col]] <- dplyr::case_when(
          grepl("icu|intensive|critical", tolower(data[[col]])) ~ "ICU",
          grepl("ward|general", tolower(data[[col]])) ~ "Ward",
          grepl("er|emergency", tolower(data[[col]])) ~ "Emergency",
          TRUE ~ data[[col]]
        )
      }
    }
  }

  completeness <- sapply(optional_cols, function(col) {
    sum(!is.na(data[[col]])) / nrow(data) * 100
  })

  completeness_df <- data.frame(
    column = names(completeness),
    completeness_pct = as.numeric(completeness),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(completeness_pct))

  message("\nOptional column completeness:")
  print(completeness_df)

  return(data)
}
# validate.R
# validate.R
# Data validation and quality control functions

#' Validate Required Fields
#'
#' Checks that all required columns are present and contain valid data.
#' Returns validation results and optionally stops on failure.
#'
#' @param data Data frame to validate
#' @param required_cols Character vector. Required column names.
#' @param stop_on_failure Logical. If TRUE, stops execution on validation failure.
#'   If FALSE, returns validation report. Default TRUE.
#' @param allow_na Logical. If TRUE, allows NA values in required columns.
#'   Default FALSE.
#' @param min_completeness Numeric. Minimum proportion of non-NA values required
#'   (0-1). Default 0.8 (80 percent completeness).
#'
#' @return List with validation results:
#' \itemize{
#'   \item \code{valid}: Logical. Overall validation status
#'   \item \code{missing_cols}: Character vector of missing columns
#'   \item \code{incomplete_cols}: Data frame of columns below min_completeness
#'   \item \code{messages}: Character vector of validation messages
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Stop on failure (default)
#' validate_required_fields(
#'   data,
#'   required_cols = c("patient_id", "organism_normalized", "antibiotic_normalized")
#' )
#'
#' # Get validation report without stopping
#' validation <- validate_required_fields(
#'   data,
#'   required_cols = c("patient_id", "organism_normalized"),
#'   stop_on_failure = FALSE
#' )
#' }
validate_required_fields <- function(data,
                                     required_cols,
                                     stop_on_failure = TRUE,
                                     allow_na = FALSE,
                                     min_completeness = 0.8) {
  validation_msgs <- character()
  missing_cols <- character()
  incomplete_cols <- data.frame()

  # Check 1: Column existence
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    msg <- sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    )
    validation_msgs <- c(validation_msgs, msg)
  }

  # Check 2: Completeness (for existing columns)
  existing_cols <- intersect(required_cols, names(data))

  if (length(existing_cols) > 0 && !allow_na) {
    completeness <- sapply(existing_cols, function(col) {
      sum(!is.na(data[[col]])) / nrow(data)
    })

    incomplete <- completeness < min_completeness

    if (any(incomplete)) {
      incomplete_cols <- data.frame(
        column = names(completeness)[incomplete],
        completeness = completeness[incomplete],
        n_missing = sapply(names(completeness)[incomplete], function(col) {
          sum(is.na(data[[col]]))
        }),
        stringsAsFactors = FALSE
      )

      msg <- sprintf(
        "Columns below %.0f%% completeness: %s",
        min_completeness * 100,
        paste(incomplete_cols$column, collapse = ", ")
      )
      validation_msgs <- c(validation_msgs, msg)
    }
  }

  # Determine overall validity
  is_valid <- length(missing_cols) == 0 && nrow(incomplete_cols) == 0

  # Create validation result
  result <- list(
    valid = is_valid,
    missing_cols = missing_cols,
    incomplete_cols = incomplete_cols,
    messages = validation_msgs,
    n_rows = nrow(data),
    n_cols_checked = length(required_cols)
  )

  # Print messages
  if (is_valid) {
    message(sprintf(
      "[v] Validation passed: All %d required columns present and complete",
      length(required_cols)
    ))
  } else {
    message("[x] Validation failed:")
    for (msg in validation_msgs) {
      message(sprintf("  - %s", msg))
    }

    if (nrow(incomplete_cols) > 0) {
      message("\nCompleteness details:")
      print(incomplete_cols)
    }
  }

  # Stop or return
  if (!is_valid && stop_on_failure) {
    stop("Data validation failed. See messages above.")
  }

  return(result)
}


#' Validate Data Quality
#'
#' Runs a set of quality checks on a data frame: minimum row count,
#' required column presence, and maximum missing-value percentage per column.
#' Returns a character vector of quality issues found, or stops if
#' \code{stop_on_failure} is \code{TRUE} and issues are detected.
#'
#' @param data             Data frame to validate.
#' @param min_rows         Integer. Minimum acceptable number of rows.
#'   Default \code{10}.
#' @param max_missing_pct  Numeric. Maximum acceptable percentage of missing
#'   values per column (0-100). Default \code{50}.
#' @param required_cols    Character vector. Column names that must be present.
#'   Default \code{c("patient_id", "organism_normalized")}.
#' @param stop_on_failure  Logical. If \code{TRUE}, stop with an error when
#'   any quality issue is found. Default \code{FALSE}.
#'
#' @return Invisibly returns the input \code{data} frame. Prints a quality
#'   report; issues are collected and (if \code{stop_on_failure}) raised as
#'   an error.
#' @export
validate_data_quality <- function(data,
                                  min_rows = 10,
                                  max_missing_pct = 50,
                                  required_cols = c("patient_id", "organism_normalized"),
                                  stop_on_failure = FALSE) {
  quality_issues <- character()

  # Check 1: Minimum rows
  n_rows <- nrow(data)
  if (n_rows < min_rows) {
    quality_issues <- c(
      quality_issues,
      sprintf("Dataset too small: %d rows (minimum: %d)", n_rows, min_rows)
    )
  }

  # Check 2: Required columns
  missing_req <- setdiff(required_cols, names(data))
  if (length(missing_req) > 0) {
    quality_issues <- c(
      quality_issues,
      sprintf("Missing required columns: %s", paste(missing_req, collapse = ", "))
    )
  }

  # Check 3: Completeness per column
  completeness <- sapply(names(data), function(col) {
    100 * sum(!is.na(data[[col]])) / n_rows
  })

  col_completeness <- data.frame(
    column = names(completeness),
    completeness_pct = as.numeric(completeness),
    n_missing = sapply(names(data), function(col) sum(is.na(data[[col]]))),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(completeness_pct)

  # Identify columns with excessive missingness
  poor_cols <- col_completeness %>%
    dplyr::filter(completeness_pct < (100 - max_missing_pct))

  if (nrow(poor_cols) > 0) {
    quality_issues <- c(
      quality_issues,
      sprintf(
        "%d columns exceed %.0f%% missing threshold: %s",
        nrow(poor_cols),
        max_missing_pct,
        paste(poor_cols$column[1:min(5, nrow(poor_cols))], collapse = ", ")
      )
    )
  }

  # Check 4: Overall completeness
  overall_completeness <- sum(!is.na(data)) / (nrow(data) * ncol(data))

  if (overall_completeness < 0.5) {
    quality_issues <- c(
      quality_issues,
      sprintf(
        "Overall completeness too low: %.1f%% (target: >=50%%)",
        overall_completeness * 100
      )
    )
  }

  # Determine pass/fail
  passes_quality <- length(quality_issues) == 0

  # Create result
  result <- list(
    passes_quality = passes_quality,
    n_rows = n_rows,
    n_cols = ncol(data),
    overall_completeness = overall_completeness,
    column_completeness = col_completeness,
    quality_issues = quality_issues
  )

  # Print results
  if (passes_quality) {
    message(sprintf(
      "[v] Quality check passed: %d rows x %d cols, %.1f%% complete",
      n_rows, ncol(data), overall_completeness * 100
    ))
  } else {
    message("[x] Quality issues detected:")
    for (issue in quality_issues) {
      message(sprintf("  - %s", issue))
    }
  }

  # Stop if requested
  if (!passes_quality && stop_on_failure) {
    stop("Data quality validation failed")
  }

  return(result)
}


#' Check Logical Consistency
#'
#' Validates logical relationships in the data (e.g., date sequences,
#' age consistency, valid value ranges).
#'
#' @param data Data frame
#' @param checks Character vector. Which checks to perform:
#'   - "date_sequence": admission < culture < outcome
#'   - "age_range": Age between 0-120
#'   - "age_dob_match": Age matches DOB
#'   - "outcome_consistency": Died patients have outcome date
#'   - "all": All checks (default)
#' @param stop_on_failure Logical. Stop on inconsistency. Default FALSE.
#'
#' @return List with consistency check results:
#'   - consistent: Logical
#'   - issues_found: Data frame of inconsistent rows
#'   - summary: Character vector of issue summaries
#'
#' @export
#'
#' @examples
#' \dontrun{
#' consistency <- check_logical_consistency(data, checks = "all")
#'
#' if (!consistency$consistent) {
#'   print(consistency$summary)
#'   View(consistency$issues_found)
#' }
#' }
check_logical_consistency <- function(data,
                                      checks = "all",
                                      stop_on_failure = FALSE) {
  # Expand "all" to all check types
  all_checks <- c("date_sequence", "age_range", "age_dob_match", "outcome_consistency")
  if ("all" %in% checks) {
    checks <- all_checks
  }

  issues <- list()
  summary_msgs <- character()

  # Check 1: Date sequence
  if ("date_sequence" %in% checks) {
    date_cols <- c("date_of_admission", "date_of_culture", "date_of_final_outcome")
    if (all(date_cols %in% names(data))) {
      date_issues <- data %>%
        dplyr::filter(
          !is.na(date_of_admission) & !is.na(date_of_culture) &
            (date_of_admission > date_of_culture |
              (!is.na(date_of_final_outcome) & date_of_culture > date_of_final_outcome))
        )

      if (nrow(date_issues) > 0) {
        issues$date_sequence <- date_issues
        summary_msgs <- c(
          summary_msgs,
          sprintf("Date sequence violations: %d rows", nrow(date_issues))
        )
      }
    }
  }

  # Check 2: Age range
  if ("age_range" %in% checks && "Age" %in% names(data)) {
    age_issues <- data %>%
      dplyr::filter(!is.na(Age) & (Age < 0 | Age > 120))

    if (nrow(age_issues) > 0) {
      issues$age_range <- age_issues
      summary_msgs <- c(
        summary_msgs,
        sprintf("Age out of range (0-120): %d rows", nrow(age_issues))
      )
    }
  }

  # Check 3: Age-DOB match
  if ("age_dob_match" %in% checks && all(c("Age", "DOB", "date_of_culture") %in% names(data))) {
    age_dob_issues <- data %>%
      dplyr::filter(
        !is.na(Age) & !is.na(DOB) & !is.na(date_of_culture)
      ) %>%
      dplyr::mutate(
        calculated_age = as.numeric(difftime(date_of_culture, DOB, units = "days")) / 365.25,
        age_diff = abs(Age - calculated_age)
      ) %>%
      dplyr::filter(age_diff > 2) # Allow 2-year tolerance

    if (nrow(age_dob_issues) > 0) {
      issues$age_dob_match <- age_dob_issues
      summary_msgs <- c(
        summary_msgs,
        sprintf("Age-DOB mismatch: %d rows (>2 year difference)", nrow(age_dob_issues))
      )
    }
  }

  # Check 4: Outcome consistency
  if ("outcome_consistency" %in% checks && all(c("final_outcome", "date_of_final_outcome") %in% names(data))) {
    outcome_issues <- data %>%
      dplyr::filter(
        final_outcome == "Died" & is.na(date_of_final_outcome)
      )

    if (nrow(outcome_issues) > 0) {
      issues$outcome_consistency <- outcome_issues
      summary_msgs <- c(
        summary_msgs,
        sprintf("Died patients without outcome date: %d rows", nrow(outcome_issues))
      )
    }
  }

  # Combine all issues
  is_consistent <- length(issues) == 0

  if (length(issues) > 0) {
    all_issues <- dplyr::bind_rows(lapply(names(issues), function(check_name) {
      issues[[check_name]] %>%
        dplyr::mutate(consistency_check = check_name)
    }))
  } else {
    all_issues <- data.frame()
  }

  # Create result
  result <- list(
    consistent = is_consistent,
    issues_found = all_issues,
    summary = summary_msgs,
    n_checks_performed = length(checks),
    n_issues = length(issues)
  )

  # Print results
  if (is_consistent) {
    message(sprintf(
      "[v] Logical consistency: All %d checks passed",
      length(checks)
    ))
  } else {
    message("[x] Logical inconsistencies found:")
    for (msg in summary_msgs) {
      message(sprintf("  - %s", msg))
    }
  }

  # Stop if requested
  if (!is_consistent && stop_on_failure) {
    stop("Logical consistency validation failed")
  }

  return(result)
}
# data.R
# Data objects and mappings for AMR preprocessing

#' Default Column Name Mappings
#'
#' Named list mapping standard column names to common aliases found in
#' AMR surveillance datasets from different sources.
#'
#' @format A named list with 15 elements
#' @export
default_column_mappings <- list(
  patient_id = c(
    "PatientInformation_id", "PatientInformation", "patient_ID",
    "Patient_ID", "PatientID", "Subject_ID", "MRN",
    "medical_record_number", "UHID", "patient_no", "Case_ID"
  ),
  gender = c(
    "Gender", "sex", "Sex", "patient_gender", "Patient_Gender",
    "gender_code", "sex_code"
  ),
  state = c(
    "State", "patient_state", "state_name", "region",
    "State_Name", "state_code"
  ),
  location = c(
    "Location", "city", "hospital_city", "hospital_location",
    "site", "City", "Hospital_Location"
  ),
  DOB = c(
    "date_of_birth", "DateOfBirth", "birth_date", "dob", "DOB",
    "Date_of_Birth", "BirthDate"
  ),
  Age = c("age", "patient_age", "Age_Years", "age_years", "AGE"),
  date_of_admission = c(
    "admission_date", "Date.of.admission",
    "Date_of_admission_in_hospital", "AdmissionDate",
    "hospital_admission_date", "admit_date",
    "Date.of.admission.in.hospital"
  ),
  date_of_culture = c(
    "date_of_event", "Date.of.event", "culture_date",
    "collection_date", "sample_date", "specimen_date",
    "event_date", "culture_collection_date",
    "Date_of_culture", "CultureDate"
  ),
  date_of_final_outcome = c(
    "Date.of.14.day.outcome", "outcome_date",
    "discharge_date", "death_date", "Date_of_outcome",
    "final_outcome_date", "DateOfOutcome",
    "date_of_discharge"
  ),
  final_outcome = c(
    "Final.outcome", "outcome", "patient_outcome",
    "status", "final_status", "discharge_status",
    "Outcome", "Status"
  ),
  organism_name = c(
    "Organism", "organism", "pathogen", "pathogen_name",
    "bacteria", "microorganism", "organism_identified",
    "OrganismName", "isolated_organism"
  ),
  antibiotic_name = c(
    "antibiotic", "drug", "drug_name", "antimicrobial",
    "antimicrobial_name", "antibiotic_tested",
    "drug_tested", "Antibiotic", "AntibioticName"
  ),
  antibiotic_value = c(
    "antibiotic_result", "susceptibility", "resistance",
    "result", "remarks", "antibiotic_remarks",
    "susceptibility_result", "Remarks", "Result",
    "susceptibility_status"
  ),
  specimen_type = c(
    "specimen", "sample_type", "sample", "Sample_type1_name",
    "source", "specimen_source", "sample_source",
    "culture_source", "SpecimenType", "SampleType"
  ),
  diagnosis = c(
    "Diagnosis", "diagnosis_1", "primary_diagnosis",
    "clinical_diagnosis", "Diagnosis_1", "ICD_code",
    "diagnosis_code"
  )
)


#' Beta-Lactam Class Hierarchy
#'
#' Ordered vector of beta-lactam classes for resistance class selection.
#' Order represents clinical hierarchy (most to least important).
#'
#' @return Character vector
#' @export
get_beta_lactam_hierarchy <- function() {
  c(
    "Carbapenems",
    "Fourth-generation-cephalosporins",
    "Third-generation-cephalosporins",
    "Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal",
    "Beta-lactam/beta-lactamase-inhibitor",
    "Aminopenicillins",
    "Penicillins"
  )
}


#' Standard Age Bins (GBD Compatible)
#'
#' Returns age bin boundaries for stratification, compatible with GBD
#' (Global Burden of Disease) methodology.
#'
#' @param type Character. "GBD_standard" for 5-year bins, "pediatric" for
#'   child-focused bins, or "geriatric" for elderly-focused bins.
#'
#' @return Character vector of age bin labels
#' @export
get_age_bins <- function(type = "GBD_standard") {
  switch(type,
    GBD_standard = c(
      "<1", "1-5", "5-10", "10-15", "15-20", "20-25", "25-30",
      "30-35", "35-40", "40-45", "45-50", "50-55", "55-60",
      "60-65", "65-70", "70-75", "75-80", "80-85", "85+"
    ),
    pediatric = c(
      "0-0.08", "0.08-1", "1-2", "2-5", "5-10", "10-15", "15-18", "18+"
    ),
    geriatric = c(
      "0-50", "50-60", "60-65", "65-70", "70-75", "75-80", "80-85",
      "85-90", "90+"
    ),
    stop("Unknown age bin type. Use 'GBD_standard', 'pediatric', or 'geriatric'")
  )
}


#' Get Magiorakos MDR/XDR Thresholds
#'
#' Returns pathogen-specific MDR/XDR classification criteria based on
#' Magiorakos et al. 2012 (Clin Microbiol Infect).
#'
#' @return Data frame with MDR/XDR thresholds per organism group
#' @export
#' @references
#' Magiorakos AP, Srinivasan A, Carey RB, et al. Multidrug-resistant,
#' extensively drug-resistant and pandrug-resistant bacteria: an international
#' expert proposal for interim standard definitions for acquired resistance.
#' Clin Microbiol Infect. 2012;18(3):268-281.
get_magiorakos_thresholds <- function() {
  tibble::tribble(
    ~organism_group, ~mdr_threshold, ~xdr_threshold, ~total_categories,
    "Enterobacterales", 3, "all_but_2", 9,
    "Pseudomonas aeruginosa", 3, "all_but_2", 10,
    "Acinetobacter spp", 3, "all_but_2", 9,
    "Staphylococcus aureus", 3, "all_but_2", 9,
    "Enterococcus spp", 3, "all_but_2", 7,
    "Streptococcus pneumoniae", 3, "all_but_2", 5
  )
}


#' Get Antimicrobial Categories for MDR/XDR Classification
#'
#' Returns antimicrobial categories used in Magiorakos MDR/XDR definitions.
#' Categories are pathogen-specific.
#'
#' @param organism_group Character. Organism group name.
#' @return Character vector of antimicrobial categories for that organism
#' @export
get_antimicrobial_categories <- function(organism_group = "Enterobacterales") {
  categories <- list(
    Enterobacterales = c(
      "Aminoglycosides",
      "Carbapenems",
      "Cephalosporins (3rd gen)",
      "Cephalosporins (4th gen)",
      "Fluoroquinolones",
      "Monobactams",
      "Penicillins + beta-lactamase inhibitors",
      "Polymyxins",
      "Tigecycline"
    ),
    "Pseudomonas aeruginosa" = c(
      "Aminoglycosides",
      "Antipseudomonal carbapenems",
      "Antipseudomonal cephalosporins",
      "Antipseudomonal fluoroquinolones",
      "Antipseudomonal penicillins + beta-lactamase inhibitors",
      "Monobactams",
      "Phosphonic acids",
      "Polymyxins",
      "Ceftazidime-avibactam",
      "Ceftolozane-tazobactam"
    ),
    "Acinetobacter spp" = c(
      "Aminoglycosides",
      "Carbapenems",
      "Cephalosporins (extended-spectrum)",
      "Fluoroquinolones",
      "Penicillins + beta-lactamase inhibitors",
      "Polymyxins",
      "Sulbactam",
      "Tetracyclines",
      "Tigecycline"
    ),
    "Staphylococcus aureus" = c(
      "Aminoglycosides",
      "Fluoroquinolones",
      "Folate pathway inhibitors",
      "Fusidic acid",
      "Glycopeptides",
      "Linezolid",
      "Mupirocin",
      "Rifampicin",
      "Tetracyclines"
    )
  )

  if (organism_group %in% names(categories)) {
    return(categories[[organism_group]])
  } else {
    warning(sprintf("No category list for '%s', using Enterobacterales default", organism_group))
    return(categories$Enterobacterales)
  }
}


#' Get Organism Taxonomy Mapping
#'
#' Reads the organism taxonomy from inst/extdata/organisms.csv and returns
#' a data frame mapping organism names to organism groups.
#'
#' @return Data frame with columns organism_name and org_group
#' @keywords internal
get_organism_taxonomy <- function() {
  file_path <- find_extdata_file("organisms.csv")
  if (file_path == "") {
    warning("organisms.csv not found. Returning empty taxonomy.")
    return(data.frame(
      organism_name = character(),
      org_group = character(),
      stringsAsFactors = FALSE
    ))
  }
  taxonomy <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  if ("organism_group" %in% names(taxonomy) && !"org_group" %in% names(taxonomy)) {
    names(taxonomy)[names(taxonomy) == "organism_group"] <- "org_group"
  }
  taxonomy[, c("organism_name", "org_group")]
}


#' Get RR Pathogen Mapping
#'
#' Returns a mapping from normalized organism names to RR pathogen categories
#' used in burden estimation. Reads from inst/extdata/organisms.csv.
#'
#' @return Data frame with columns organism_name and rr_pathogen
#' @keywords internal
