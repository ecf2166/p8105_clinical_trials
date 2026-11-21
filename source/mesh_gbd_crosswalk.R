# MeSH -> GBD crosswalk.
library(tidyverse)
library(stringr)

# clean and normalize MeSH terms for matching. Lowercase, trim whitespace, and collapse internal whitespace to a single space.
normalize_mesh = function(x) {
  str_squish(str_to_lower(x))
}

# Map each GBD cause to its level-2 ancestor and analysis category.
# Level-2 causes are the main disease categories used in GBD analyses. 
# The analysis category is either the level-2 cause or Injuries (687) for any injury-related cause.

## Find each GBD cause's borader category 
gbd_targets = function(hierarchy) {
  stopifnot(!anyDuplicated(hierarchy$cause_id))
  ids = as.integer(hierarchy$cause_id)
  find_level2 = function(id) {
    visited = integer()
    repeat {
      i = match(id, ids)
      if (is.na(i) || id %in% visited) return(NA_integer_)
      if (hierarchy$level[i] == 2) return(as.integer(id))
      if (hierarchy$level[i] < 2) return(NA_integer_)
      visited = c(visited, id)
      id = as.integer(hierarchy$parent_id[i])
    }
  }
  hierarchy |>
    mutate(
      level2_id = vapply(cause_id, find_level2, integer(1)),
      level2_cause = cause_name[match(level2_id, cause_id)],
      analysis_cause_id = if_else(
        level2_id %in% c(688L, 696L, 717L) | cause_id == 687L,
        687L, level2_id),
      analysis_cause = cause_name[match(analysis_cause_id, cause_id)]
    )
}

# Filter studies to a cohort of interest. Default as of now is interventional trials starting 2010-2023.
crosswalk_cohort = function(studies, start_year = 2010L, end_year = 2023L,
                            study_types = "INTERVENTIONAL") {
  stopifnot(start_year <= end_year, !anyDuplicated(studies$nct_id))
  studies |>
    mutate(start_date = as.Date(start_date)) |>
    filter(
      study_type %in% study_types,
      !is.na(start_date),
      start_date >= as.Date(sprintf("%d-01-01", start_year)),
      start_date < as.Date(sprintf("%d-01-01", end_year + 1L)))
}

# Find all MeSH terms associated with the cohort.
direct_mesh = function(browse_conditions, cohort) {
  browse_conditions |>
    filter(mesh_type_condition == "mesh-list") |>
    semi_join(cohort, by = "nct_id") |>
    transmute(nct_id, mesh_term = normalize_mesh(mesh_term_condition)) |>
    filter(!is.na(mesh_term), mesh_term != "") |>
    distinct()
}

# Generates table of proposed mapping per term. 
# Review decisions: pending, approve, exclude.
# To resolve a term, enter reviewed_gbd_cause_id and reviewer, then approve.
build_mesh_review = function(browse_conditions, cohort, hierarchy, proposed_rules,
                             previous_review = NULL) {
  gbd_categories = gbd_targets(hierarchy)
  proposed_rules = proposed_rules |>
    mutate(mesh_term = normalize_mesh(mesh_term),
                  target_cause_id = as.integer(target_cause_id))
  stopifnot(!anyDuplicated(proposed_rules$mesh_term),
            all(proposed_rules$action %in% c("map", "review", "exclude")))
  if (any(!is.na(proposed_rules$target_cause_id) &
          !proposed_rules$target_cause_id %in% gbd_categories$cause_id)) {
    stop("A proposed rule references a cause ID absent from the hierarchy.")
  }
  terms = browse_conditions |>
    filter(mesh_type_condition == "mesh-list") |>
    transmute(mesh_term = normalize_mesh(mesh_term_condition)) |>
    filter(!is.na(mesh_term), mesh_term != "") |>
    distinct()
  frequency = direct_mesh(browse_conditions, cohort) |>
    count(mesh_term, name = "n_trials_in_cohort")
  exact = gbd_categories |>
    filter(!is.na(analysis_cause_id)) |>
    transmute(mesh_term = normalize_mesh(cause_name),
                     exact_cause_id = as.integer(cause_id))
  stopifnot(!anyDuplicated(exact$mesh_term))
  review = terms |>
    left_join(frequency, by = "mesh_term") |>
    mutate(n_trials_in_cohort = replace_na(n_trials_in_cohort, 0L)) |>
    left_join(proposed_rules, by = "mesh_term") |>
    left_join(exact, by = "mesh_term") |>
    mutate(
      proposed_gbd_cause_id = if_else(
        !is.na(action), target_cause_id, exact_cause_id),
      proposed_action = case_when(
        !is.na(action) ~ action,
        !is.na(exact_cause_id) ~ "map",
        TRUE ~ "review"),
      mapping_method = case_when(
        !is.na(action) ~ "explicit_project_rule",
        !is.na(exact_cause_id) ~ "normalized_label_match",
        TRUE ~ "unresolved"),
      rationale = coalesce(rationale,
        if_else(!is.na(exact_cause_id),
          "Normalized label matches a GBD cause; equivalence still requires review.",
          "No proposed mapping. Review the mesh term and trial condition context."))) |>
    select(mesh_term, n_trials_in_cohort, proposed_action,
                  proposed_gbd_cause_id, mapping_method, rationale) |>
    left_join(gbd_categories |>
      select(proposed_gbd_cause_id = cause_id,
                    proposed_gbd_cause = cause_name,
                    proposed_level2_id = level2_id,
                    proposed_level2_cause = level2_cause,
                    proposed_analysis_cause = analysis_cause),
      by = "proposed_gbd_cause_id")
  fields = c("mesh_term", "review_decision", "reviewed_gbd_cause_id",
             "reviewer", "review_note")
  if (!is.null(previous_review)) {
    stopifnot(!anyDuplicated(previous_review$mesh_term),
              all(fields %in% names(previous_review)))
    previous_review = previous_review |>
      select(all_of(fields)) |>
      mutate(reviewed_gbd_cause_id = as.integer(reviewed_gbd_cause_id))
    # Retain previously reviewed terms even if the current raw extract lacks them.
    review = full_join(review, previous_review, by = "mesh_term")
  } else {
    review = review |>
      mutate(review_decision = "pending", reviewed_gbd_cause_id = NA_integer_,
                    reviewer = NA_character_, review_note = NA_character_)
  }
  review |>
    mutate(review_decision = replace_na(review_decision, "pending"),
                  n_trials_in_cohort = replace_na(n_trials_in_cohort, 0L)) |>
    arrange(desc(n_trials_in_cohort), mesh_term)
}

# Resolves the review decisions into a final mapping for each term. 
# Approved or excluded terms are resolved; 
# Pending terms remain unresolved unless allow_proposed is TRUE and the proposed action is "map". 
# The resolved cause must be a GBD cause with a level-2 ancestor or Injuries (687).
resolve_mesh_review = function(review, hierarchy, allow_proposed = FALSE) {
  gbd_categories = gbd_targets(hierarchy)
  stopifnot(!anyDuplicated(review$mesh_term))
  if (any(is.na(review$review_decision) |
          !review$review_decision %in% c("pending", "approve", "exclude"))) {
    stop("review_decision must be pending, approve, or exclude.")
  }
  decisions = review$review_decision != "pending"
  if (any(decisions & (is.na(review$reviewer) | trimws(review$reviewer) == ""))) {
    stop("Approved/excluded decisions require a reviewer name or initials.")
  }
  out = review |>
    mutate(
      resolved_cause_id = case_when(
        review_decision == "approve" ~ coalesce(
          as.integer(reviewed_gbd_cause_id), as.integer(proposed_gbd_cause_id)),
        review_decision == "pending" & allow_proposed & proposed_action == "map" ~
          as.integer(proposed_gbd_cause_id),
        TRUE ~ NA_integer_),
      mapping_status = case_when(
        review_decision == "exclude" ~ "excluded_reviewed",
        review_decision == "approve" ~ "mapped_reviewed",
        !is.na(resolved_cause_id) ~ "mapped_provisional",
        proposed_action == "exclude" ~ "excluded_proposed",
        TRUE ~ "needs_review"))
  invalid = !is.na(out$resolved_cause_id) &
    !out$resolved_cause_id %in% gbd_categories$cause_id[!is.na(gbd_categories$analysis_cause_id)]
  if (any(invalid) || any(out$review_decision == "approve" & is.na(out$resolved_cause_id))) {
    stop("Approved/proposed target must be a valid GBD cause with a level-2 ancestor, or Injuries (687).")
  }
  out |>
    left_join(gbd_categories |>
      select(resolved_cause_id = cause_id, resolved_cause_name = cause_name,
                    level2_id, level2_cause, analysis_cause_id, analysis_cause),
      by = "resolved_cause_id")
}

# Applies the MeSH to GBD crosswalk to the trial data
# Produces one row per trial per cause, with supporting MeSH terms and coverage info
apply_mesh_crosswalk = function(browse_conditions, cohort, review, hierarchy,
                                allow_proposed = FALSE) {
  term_map = resolve_mesh_review(review, hierarchy, allow_proposed)
  evidence = direct_mesh(browse_conditions, cohort) |>
    left_join(term_map, by = "mesh_term") |>
    mutate(mapping_status = replace_na(mapping_status, "needs_review"))
  mapped = evidence |> filter(!is.na(analysis_cause_id))
  trial_causes = mapped |>
    group_by(nct_id, analysis_cause_id, analysis_cause) |>
    summarize(
      n_supporting_terms = n_distinct(mesh_term),
      supporting_terms = paste(sort(unique(mesh_term)), collapse = "; "),
      all_terms_reviewed = all(mapping_status == "mapped_reviewed"),
      .groups = "drop") |>
    group_by(nct_id) |>
    mutate(n_analysis_causes = n()) |>
    ungroup()
  flags = evidence |>
    group_by(nct_id) |>
    summarize(
      n_direct_terms = n_distinct(mesh_term),
      n_mapped_terms = sum(!is.na(analysis_cause_id)),
      n_unresolved_terms = sum(mapping_status == "needs_review"),
      n_excluded_terms = sum(mapping_status %in% c("excluded_reviewed", "excluded_proposed")),
      .groups = "drop")
  any_mesh_ids = browse_conditions |>
    filter(!is.na(mesh_term_condition), trimws(mesh_term_condition) != "") |>
    distinct(nct_id) |> pull(nct_id)
  
# Produce a trial-level audit table with coverage status based on the number of direct terms, mapped terms, unresolved terms, and excluded terms.
  audit = cohort |>
    select(nct_id, study_type, start_date) |>
    left_join(flags, by = "nct_id") |>
    mutate(
      across(c(n_direct_terms, n_mapped_terms, n_unresolved_terms, n_excluded_terms),
                    ~ replace_na(.x, 0L)),
      coverage_status = case_when(
        n_direct_terms == 0 & nct_id %in% any_mesh_ids ~ "ancestor_only",
        n_direct_terms == 0 ~ "no_mesh",
        n_mapped_terms == 0 ~ "no_mapped_direct_term",
        n_unresolved_terms > 0 ~ "mapped_with_unresolved_terms",
        TRUE ~ "mapped_no_unresolved_terms"))
  list(trial_causes = trial_causes, term_evidence = evidence,
       trial_audit = audit, resolved_crosswalk = term_map)
}
