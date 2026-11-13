Data
================
2026-11-11

# Loading AACT data

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
source("source/read_aact.R")

studies = read_aact("studies", col_select = c(nct_id, brief_title, official_title, acronym, study_type, start_date, start_date_type, completion_date, overall_status, why_stopped, phase, enrollment, enrollment_type, source, source_class)) |>
  mutate(
    start_date = as.Date(start_date),
    completion_date = as.Date(completion_date),
    enrollment = as.numeric(enrollment)
  )

countries = read_aact("countries", col_select = c(nct_id, name, removed)) |>
  rename(country = name) |>
  # AACT may store booleans as "t"/"f"; this handles both that and "true"/"false"
  mutate(removed = str_to_lower(removed) == "t")

sponsors = read_aact("sponsors", col_select = c(nct_id, name, agency_class, lead_or_collaborator)) |>
  rename(sponsor_name = name)

conditions = read_aact("conditions", col_select = c(nct_id, downcase_name)) |>
  rename(condition_name = downcase_name)

browse_conditions = read_aact("browse_conditions", col_select = c(nct_id, downcase_mesh_term, mesh_type)) |>
  rename(mesh_term_condition = downcase_mesh_term,
         mesh_type_condition = mesh_type)

browse_interventions = read_aact("browse_interventions", col_select = c(nct_id, downcase_mesh_term, mesh_type)) |>
  rename(mesh_term_intervention = downcase_mesh_term,
         mesh_type_intervention = mesh_type)


# Baseline participants (only trials that posted results)
# Largest participant count per trial = the Total group

baseline_counts = read_aact("baseline_counts", col_select = c(nct_id, units, scope, count)) |>
  filter(units == "Participants", scope == "overall") |>
  mutate(count = as.numeric(count)) |>
  group_by(nct_id) |>
  summarize(baseline_n = max(count))



# number of countries removed (removed = TRUE)
count(countries, removed)
```

    ## # A tibble: 2 × 2
    ##   removed      n
    ##   <lgl>    <int>
    ## 1 FALSE   790664
    ## 2 TRUE     36138

Saving data

``` r
# don't want to knit this every time, only run once
# for keeping data in in RDS format to save space (and so this knits)
write_rds(studies, "data/studies.rds", compress = "gz")
write_rds(countries, "data/countries.rds", compress = "gz")
write_rds(sponsors, "data/sponsors.rds", compress = "gz")
write_rds(conditions, "data/conditions.rds", compress = "gz")
write_rds(browse_conditions, "data/browse_conditions.rds", compress = "gz")
write_rds(browse_interventions, "data/browse_interventions.rds", compress = "gz")
write_rds(baseline_counts, "data/baseline_counts.rds", compress = "gz")
```

# Conditions

``` r
# distinguishing top-level MeSH disease categories (https://www.nlm.nih.gov/mesh/intro_trees.html, https://meshb.nlm.nih.gov/treeView)
top_level_dz = c(
  "infections",                                   # C01
  "neoplasms",                                    # C04
  "musculoskeletal diseases",                     # C05
  "digestive system diseases",                    # C06
  "stomatognathic diseases",                      # C07
  "respiratory tract diseases",                   # C08
  "otorhinolaryngologic diseases",                # C09
  "nervous system diseases",                      # C10
  "eye diseases",                                 # C11
  "urogenital diseases",                          # C12
  "cardiovascular diseases",                      # C14
  "hemic and lymphatic diseases",                 # C15
  "congenital, hereditary, and neonatal diseases and abnormalities",  # C16
  "skin and connective tissue diseases",          # C17
  "nutritional and metabolic diseases",           # C18
  "endocrine system diseases",                    # C19
  "immune system diseases",                       # C20
  "disorders of environmental origin",            # C21
  "wounds and injuries",                          # C26
  "chemically-induced disorders",                 # C25
  "occupational diseases",                        # C24
  "mental disorders"                              # F03
)

# potentially-important second-level MeSH terms (based on WHO DALY categories)
second_level = c(
  "female urogenital diseases and pregnancy complications",
  "respiratory tract infections"
  # add others as the WHO crosswalk needs them
)

browse_conditions |>
  filter(mesh_term_condition %in% top_level_dz) |>
  distinct(nct_id, mesh_term_condition) |>
  count(mesh_term_condition, sort = TRUE)
```

    ## # A tibble: 22 × 2
    ##    mesh_term_condition                     n
    ##    <chr>                               <int>
    ##  1 neoplasms                          107074
    ##  2 nervous system diseases             76732
    ##  3 cardiovascular diseases             64877
    ##  4 urogenital diseases                 54851
    ##  5 respiratory tract diseases          52459
    ##  6 nutritional and metabolic diseases  50446
    ##  7 mental disorders                    49608
    ##  8 immune system diseases              47561
    ##  9 digestive system diseases           46498
    ## 10 infections                          46462
    ## # ℹ 12 more rows

``` r
classified = browse_conditions |>
  filter(mesh_term_condition %in% top_level_dz) |>
  distinct(nct_id)

coverage = studies |>
  mutate(group = case_when(
    nct_id %in% pull(classified, nct_id) ~ "Classified into at least one category",
    nct_id %in% pull(browse_conditions, nct_id) ~ "Have MeSH terms, but none in our categories",
    TRUE ~ "No MeSH condition terms at all"
  )) |>
  count(group) |>
  mutate(pct = n / sum(n)) |>
  arrange(desc(n))

coverage |>
  mutate(n = scales::comma(n),
         pct = scales::percent(pct)) |>
  knitr::kable(col.names = c("Group", "Trials", "Share"))
```

| Group                                       | Trials  | Share |
|:--------------------------------------------|:--------|:------|
| Classified into at least one category       | 439,522 | 72%   |
| No MeSH condition terms at all              | 129,732 | 21%   |
| Have MeSH terms, but none in our categories | 41,979  | 7%    |

``` r
unique_conditions = unique(browse_conditions$mesh_term)
```

    ## Warning: Unknown or uninitialised column: `mesh_term`.

``` r
browse_conditions |> 
  count(mesh_term_condition, sort=T) |> 
  slice_tail(n=50)
```

    ## # A tibble: 50 × 2
    ##    mesh_term_condition                              n
    ##    <chr>                                        <int>
    ##  1 toxocara canis infection (canine roundworms)     1
    ##  2 tracheobronchomegaly                             1
    ##  3 tremor hereditary essential, 1                   1
    ##  4 trichohepatoenteric syndrome                     1
    ##  5 trichophyton infection                           1
    ##  6 tricuspid valve stenosis                         1
    ##  7 trimethylaminuria                                1
    ##  8 triphalangeal thumb                              1
    ##  9 triploidy                                        1
    ## 10 tropical calcific pancreatitis                   1
    ## # ℹ 40 more rows

``` r
browse_conditions |> 
  count(nct_id) |> 
  summarize(
    mean = mean(n)
  ) |> 
  nrow()
```

    ## [1] 1

``` r
rows = nrow(browse_conditions)

browse_conditions |> 
  filter(mesh_type_condition == "mesh-ancestor") |>
  count(mesh_term_condition, sort = T) |> 
  head(30)
```

    ## # A tibble: 30 × 2
    ##    mesh_term_condition                              n
    ##    <chr>                                        <int>
    ##  1 pathological conditions, signs and symptoms 144522
    ##  2 neoplasms                                    96256
    ##  3 pathologic processes                         77370
    ##  4 nervous system diseases                      75975
    ##  5 signs and symptoms                           75907
    ##  6 neoplasms by site                            72772
    ##  7 cardiovascular diseases                      59509
    ##  8 urogenital diseases                          54793
    ##  9 respiratory tract diseases                   52049
    ## 10 nutritional and metabolic diseases           50375
    ## # ℹ 20 more rows

``` r
unique(sponsors$agency_class)
```

    ##  [1] "INDUSTRY"  "OTHER"     "OTHER_GOV" "NIH"       "UNKNOWN"   "FED"      
    ##  [7] "NETWORK"   "AMBIG"     "INDIV"     NA

``` r
studies <- read_aact("studies") 
```

Next steps:

- download and merge the WHO DALY data
