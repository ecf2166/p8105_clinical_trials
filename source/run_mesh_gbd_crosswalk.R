
# Run this script to generate the draft crosswalk outputs.
# Default: interventional studies starting 2010-2023.

# load packages and functions
library(dplyr)
library(readr)
library(tools)
source("source/mesh_gbd_crosswalk.R")

# default start and end year (can be changed)
start_year = 2010
end_year = 2023

# define and create output directory
base = "data/crosswalk" # stores input files here
out = file.path(base, "generated")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# load input files
hierarchy = read_csv(
  file.path(base, "gbd_2023_cause_hierarchy.csv"), 
  show_col_types = FALSE
)
proposed_rules = read_csv(
  file.path(base, "mesh_gbd_proposed_rules.csv"), 
  show_col_types = FALSE
)
studies = readRDS("data/studies.rds")
browse = readRDS("data/browse_conditions.rds")

# selects cohort of trials based on start and end year
cohort = crosswalk_cohort(studies, start_year, end_year)

# read previous review decisions file (if it exists, otherwise set to NULL)
review_path = file.path(base, "mesh_gbd_review.csv")
old = if (file.exists(review_path)) read_csv(review_path,
  col_types = cols(review_decision = col_character(),
    reviewer = col_character(), review_note = col_character(),
    reviewed_gbd_cause_id = col_integer()), show_col_types = FALSE) else NULL

# rebuild the review dataframe with any new terms, preserving old decisions
review = build_mesh_review(browse, cohort, hierarchy, proposed_rules, old)

# backup and save review table
if (!is.null(old) && any(old$review_decision != "pending")) {
  file.copy(review_path, file.path(out, paste0("review_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")))
}
write_csv(review, review_path, na = "")

# save the category reference tables
gbd_categories = gbd_targets(hierarchy)
write_csv(gbd_categories, 
          file.path(out, "gbd_hierarchy_with_targets.csv"), na = ""
)
write_csv(
  gbd_categories |>
    filter(!is.na(analysis_cause_id)) |>
    distinct(analysis_cause_id, analysis_cause),
  file.path(out, "analysis_categories.csv")
          )

# apply the crosswalk to the cohort using the review decisions and save the trial-level results
result = apply_mesh_crosswalk(
  browse, cohort, review, hierarchy, 
  allow_proposed = TRUE)
saveRDS(result$trial_causes, 
  file.path(out, "trial_causes_PROVISIONAL.rds"), compress = "gzip")
saveRDS(result$trial_audit, 
        file.path(out, "trial_mapping_audit.rds"), compress = "gzip")

# summarize coverage and save the summary and unresolved/proposed mappings
coverage = result$trial_audit |> 
  count(coverage_status, name = "n_trials") |>
  mutate(percent = 100 * n_trials / sum(n_trials))
write_csv(coverage, file.path(out, "coverage_summary.csv"))

# save trial count per cause and unresolved/proposed mappings for review
# note: these counts are not unique trials, as a trial can have multiple terms and multiple categories
# however, the counts are useful for prioritizing review of high-frequency terms and categories
# and trial-category rows are deduplicated so a trial counts once within a category
write_csv(result$trial_causes |> 
            count(analysis_cause_id, analysis_cause, name = "n_trials"),
          file.path(out, "category_counts_PROVISIONAL.csv"))

# save unresolved terms and proposed mappings for review
write_csv(result$resolved_crosswalk |> 
            filter(mapping_status == "needs_review") |>
            arrange(desc(n_trials_in_cohort)), 
          file.path(out, "unresolved_terms.csv"), na = "")
write_csv(result$resolved_crosswalk |> 
            filter(mapping_status == "mapped_provisional") |>
            arrange(desc(n_trials_in_cohort)), 
          file.path(out, "proposed_mappings_to_check.csv"), na = "")

# save a few example trials for each unresolved term to help with review and context
examples = result$trial_terms |> 
  filter(mapping_status == "needs_review") |>
  arrange(mesh_term, nct_id) |> 
  group_by(mesh_term) |> 
  slice_head(n = 3) |> 
  ungroup() |>
  select(mesh_term, nct_id) |>
  left_join(studies |> select(nct_id, brief_title), by = "nct_id") |>
  mutate(record_url = paste0("https://clinicaltrials.gov/study/", nct_id))
write_csv(examples, file.path(out, "unresolved_trial_examples.csv"), na = "")

# Check unique trial-category rows and that coverage counts account for every cohort trial
stopifnot(
  !anyDuplicated(result$trial_causes[c("nct_id", "analysis_cause_id")]
  ),
  sum(coverage$n_trials) == nrow(cohort)
)


write_csv(tibble(
  start_year, end_year, 
  study_types = "INTERVENTIONAL", 
  n_cohort = nrow(cohort),
  n_review_terms = nrow(review), 
  n_mapped_trials = n_distinct(result$trial_causes$nct_id),
  n_approved_terms = sum(review$review_decision == "approve"),
  run_time = format(Sys.time(), tz = "America/New_York", usetz = TRUE)),
  file.path(out, "run_summary.csv"))

print(coverage)
