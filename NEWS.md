# anumaan 0.1.0.9036

* Prepared for CRAN submission: every exported example now runs, internal
  helpers no longer produce Rd pages, and console tables printed during
  preprocessing are emitted with `message()` so `suppressMessages()` silences
  them.
* Pathway 1 now works with osqp >= 1.0. osqp 1.0 removed the `$solve()`
  method, so every osqp solve errored and `estimate_profiles_convex()`,
  `compute_resistance_profiles()` and `bootstrap_profiles_convex()` silently
  returned uniform profile probabilities flagged as non-converged. They now call
  `osqp::solve_osqp()`, which works across osqp versions.
* `compute_marginals_from_data()` and `compute_pairwise_from_data()` no longer
  drop antibiotic classes whose names contain hyphens, slashes or spaces (e.g.
  "Third-generation-cephalosporins"), which `preprocess_for_profiles()` had
  renamed when pivoting to wide format.
* Organism, drug and syndrome names are now compared case-insensitively through
  one lower-case key. As a result `org_group` and `rr_pathogen` are populated for
  `run_preprocess()` output (previously "Other" and `NA` for every row), and
  `panel_map`, `syndrome_name` and `pathogen_name` arguments match regardless of
  case; outputs keep the data's own labels. The derived keys `rr_pathogen` and
  the `pathogen` column of `daly_load_rr_reference()` are now lower case.
* Fixed `inst/extdata/adjustment_for_CFR.csv`: an unquoted comma in "Infections
  of bones, joints, and related organs" corrupted that row.
* Functions that take a `seed` argument now restore the caller's random number
  generator state on exit.
* `preprocess_for_profiles()` no longer errors at the deduplication step; it now
  deduplicates isolate x drug rows before collapsing to antibiotic class level.
* `daly_calc_yll_associated()` and `daly_load_rr_reference()` now read the
  bundled CSV reference tables (the `.xlsx` files they expected are no longer
  shipped).

# anumaan 0.1.0.9021

* Expanded the specimen classification reference with stewardship sample labels
  seen in 2021-era extracts, including brain abscess, instrument, lung aspirate,
  lymph node, hair, and superficial biopsy.
* Collapsed `sterile_classification` outputs from detailed reference labels to
  `Sterile`, `Non-Sterile`, and `Others/Ambiguous` for downstream analysis.
* Made EDA and burden plotting functions usable without `ggpubr`; when `ggpubr`
  is not installed, `eda_theme()` falls back to `ggplot2::theme_bw()`.
* Refreshed vignettes and installation guidance for the 0.1.0.9021 workflow.
