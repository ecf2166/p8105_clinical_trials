MeSH to GBD: crosswalk review
================
Sanyu Takirambudde
2026-11-21

This draft uses the IHME GBD 2023 Cause Hierarchy and the assigned
`mesh-list` condition terms in our AACT extract. The proposed
MeSH-to-GBD mappings are **tentative** and still require review.

The example cohort is interventional studies starting in 2010–2023. This
is adjustable and not final as we can extend it to include more years.

## What we are mapping

MeSH condition term → GBD cause → broader GBD category

We first link a **MeSH condition term** to a compatible GBD cause ID,
then follow GBD’s parent IDs to its Level 2 category.
`source/tidy_gbd.R` combines the three injury Level 2 groups into
`Injuries`, so this crosswalk uses that same 20 categories for analysis

A MeSH tree branch is not equivalent to a GBD cause group. For example,
organ system headings can contain infections, cancer, congenital
conditions and injuries.

## Files and workflow

1.  `data/crosswalk/gbd_2023_cause_hierarchy.csv`: extracted from IHME’s
    Cause Hierarchy sheet, with its original cause/parent identifiers.
2.  `data/crosswalk/mesh_gbd_proposed_rules.csv`: ***suggested***
    crosswalk for matching MeSH terms with GBD causes. Some terms
    flagged for review (unresolved/excluded).
3.  `data/crosswalk/mesh_gbd_review.csv`: working table with one row per
    MeSH term, ordered by unique trial frequency in the cohort.
4.  `source/mesh_gbd_crosswalk.R`: reusable functions and definitions.
5.  `source/run_mesh_gbd_crosswalk.R`: calls those functions and
    produces the review table and related coverage outputs.

## How to run pipeline

``` r
source("source/run_mesh_gbd_crosswalk.R")
```

## How to review a term

Start with frequent terms in
`data/crosswalk/generated/proposed_mappings_to_check.csv`. These are
terms with proposed mappings that have not been approved. Next check
`data/crosswalk/generated/unresolved_terms.csv` for terms that still
need a mapping decision. Keep in mind that these counts overlap across
trials so summing term frequencies is not a count of unique trials. If
need be, open example trial records in
`data/crosswalk/generated/unresolved_trial_examples.csv` when context is
needed.

Record decisions by editing the **main** `mesh_gbd_review.csv`, not the
generated ones:

- `review_decision = pending`: default; no final decision yet.
- `review_decision = approve`: accept the proposed cause or enter a
  replacement `reviewed_gbd_cause_id` if you decide to change it. Add
  your initials/name. You can also add a review note if you want.
- `review_decision = exclude`: means the term can’t be assigned a
  disease cause. Add initials/name. Add reason if you want. This
  excludes the term, not the entire trial.

## Difficult decisions/ambiguious terms

For an unmapped term, one option is to look for its concept in GBD’s
cause hierarchy/definitions and then use a compatible parent cause when
the subtype can’t be identified. If the term crosses Level 2 causes, we
can leave it pending or develop a documented contextual rule.
`Other non-communicable diseases` is a defined GBD category, so we
should not use it as a catch-all for unmatched MeSH terms. We can then
check the high-frequency ‘rules’ and all cross-category decisions. Terms
such as hypertension, anemia, sepsis, depression and chronic pain
require particular care about risk factors, symptoms, etiology and
sequelae.

## Coverage and unresolved work

``` r
library(tidyverse)
source("source/mesh_gbd_crosswalk.R")
```

``` r
coverage_path = "data/crosswalk/generated/coverage_summary.csv" 

if (file.exists(coverage_path)) {
  read_csv(coverage_path, show_col_types = FALSE) |>
    mutate(percent = round(percent, 1)) |>
    knitr::kable()
}
```

| coverage_status              | n_trials | percent |
|:-----------------------------|---------:|--------:|
| mapped_no_unresolved_terms   |    95249 |    32.5 |
| mapped_with_unresolved_terms |    43619 |    14.9 |
| no_mapped_direct_term        |    89570 |    30.6 |
| no_mesh                      |    64660 |    22.1 |

`mapped_with_unresolved_terms` is only partial coverage: at least one
disease category was assigned, but another term remains unresolved. Even
`mapped_no_unresolved_terms` is still provisional when the mappings are
pending. `no_mapped_direct_term`, `ancestor_only` and `no_mesh` remain
separate groups to be reviewed/reported. We should not assign these
trials to “Other” or drop them without reporting how that affects the
analysis denominator.

## Reviewed mappings to be used in the analysis

``` r
hierarchy = read_csv("data/crosswalk/gbd_2023_cause_hierarchy.csv")
review = read_csv("data/crosswalk/mesh_gbd_review.csv",
  col_types = cols(review_decision = col_character(),
    reviewed_gbd_cause_id = col_integer(), reviewer = col_character(),
    review_note = col_character()))
studies = read_rds("data/studies.rds")
browse = read_rds("data/browse_conditions.rds")
cohort = crosswalk_cohort(studies, 2010, 2023)

# FALSE is the default so only approved mappings are used.
result = apply_mesh_crosswalk(browse, cohort, review, hierarchy,
                             allow_proposed = FALSE)
trial_causes = result$trial_causes
trial_audit = result$trial_audit
```

The output has one row per `nct_id × analysis_cause_id`, with the terms
corresponding to each assignment. A trial can study more than one
category. We also need to come up with a reasonable denominator and
treatment of partially mapped trials before comparing research shares
with burden shares. Next we need to link countries, deduplicate
trial–country–category combinations and decide how to count
multinational trials.

## Sources

- IHME GBD cause hierarchy:
  `IHME_GBD_2023_HIERARCHIES_Y2025M10D23.XLSX`, sheet `Cause Hierarchy`.
- [NLM MeSH tree
  structures](https://www.nlm.nih.gov/mesh/intro_trees.html).
- [GBD 2023 disease and injury burden
  study](https://pubmed.ncbi.nlm.nih.gov/41092926/).
- [AACT](https://aact.ctti-clinicaltrials.org/landing);
  `browse_conditions.rds`.
