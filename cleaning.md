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
  ) |> 
  filter(!is.na(start_date))

countries = read_aact("countries", col_select = c(nct_id, name, removed)) |>
  rename(country = name) |>
  mutate(removed = str_to_lower(removed) == "t") |> 
  filter(removed == FALSE) |> 
  select(-removed)

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

baseline_counts = read_aact("baseline_counts", col_select = c(nct_id, units, count)) |>
  filter(units == "Participants") |>
  mutate(count = as.numeric(count)) |>
  group_by(nct_id) |>
  summarize(baseline_n = max(count))
```

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
| Classified into at least one category       | 435,908 | 72%   |
| No MeSH condition terms at all              | 128,107 | 21%   |
| Have MeSH terms, but none in our categories | 41,835  | 7%    |

``` r
unique_conditions = unique(browse_conditions$mesh_term_condition)

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
  )
```

    ## # A tibble: 1 × 1
    ##    mean
    ##   <dbl>
    ## 1  9.19

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
studies |>
  mutate(group = case_when(
    nct_id %in% pull(classified, nct_id) ~ "Classified into at least one category",
    nct_id %in% pull(browse_conditions, nct_id) ~ "Have MeSH terms, but none in our categories",
    TRUE ~ "No MeSH condition terms at all"
  )) |> 
  filter(group == "No MeSH condition terms at all") |>
  select(nct_id, brief_title,enrollment, group)
```

    ## # A tibble: 128,107 × 4
    ##    nct_id      brief_title                                      enrollment group
    ##    <chr>       <chr>                                                 <dbl> <chr>
    ##  1 NCT03124173 Cognitive Control of Language                            23 No M…
    ##  2 NCT06401486 DOuble-Lumen Intubation With VIdeolaryngoscopy          916 No M…
    ##  3 NCT02173197 Effectiveness of Occupational Therapy Intervent…         20 No M…
    ##  4 NCT06632483 An Observational Study to Learn More About How …        604 No M…
    ##  5 NCT04910815 Feasibility of a New Diagnostic Device to Asses…        110 No M…
    ##  6 NCT06763536 Oral Nicotine Pouches Versus Nicotine Replaceme…       1000 No M…
    ##  7 NCT06955403 Acuvue Oaysis Contact Lens Wearers Being Refit …         40 No M…
    ##  8 NCT06479330 Evaluating the Impact of Incremental Doses of a…         60 No M…
    ##  9 NCT02556541 Ultrasound-guided Peripheral Vascular Access in…         50 No M…
    ## 10 NCT01255306 The Influence of Silicone Oil on Nerve Fiber La…         60 No M…
    ## # ℹ 128,097 more rows

Exploring studies with missing MeSH terms for conditions

``` r
eligibilities = read_aact("eligibilities", col_select = c(nct_id, healthy_volunteers))
designs       = read_aact("designs", col_select = c(nct_id, primary_purpose))

study_groups = studies |>
  mutate(group = case_when(
    nct_id %in% pull(classified, nct_id)        ~ "Classified",
    nct_id %in% pull(browse_conditions, nct_id) ~ "Other MeSH",
    TRUE                                        ~ "No MeSH"
  )) |>
  left_join(eligibilities, by = "nct_id") |>
  left_join(designs, by = "nct_id")

# Compare distributions across groups, one variable at a time
study_groups |>
  count(group, study_type) |>
  group_by(group) |>
  mutate(pct = n / sum(n)) |>
  ungroup()
```

    ## # A tibble: 8 × 4
    ##   group      study_type           n      pct
    ##   <chr>      <chr>            <int>    <dbl>
    ## 1 Classified EXPANDED_ACCESS     53 0.000122
    ## 2 Classified INTERVENTIONAL  330453 0.758   
    ## 3 Classified OBSERVATIONAL   105402 0.242   
    ## 4 No MeSH    EXPANDED_ACCESS     13 0.000101
    ## 5 No MeSH    INTERVENTIONAL   99261 0.775   
    ## 6 No MeSH    OBSERVATIONAL    28833 0.225   
    ## 7 Other MeSH INTERVENTIONAL   33824 0.809   
    ## 8 Other MeSH OBSERVATIONAL     8011 0.191

``` r
no_mesh_ids = study_groups |> filter(group == "No MeSH") |> pull(nct_id)

conditions |>
  filter(nct_id %in% no_mesh_ids) |>
  count(condition_name, sort = TRUE) |>
  head(40)
```

    ## # A tibble: 40 × 2
    ##    condition_name            n
    ##    <chr>                 <int>
    ##  1 healthy                8587
    ##  2 healthy volunteers     2274
    ##  3 anesthesia              878
    ##  4 healthy subjects        777
    ##  5 hiv                     750
    ##  6 pregnancy               647
    ##  7 aging                   612
    ##  8 healthy participants    598
    ##  9 healthy volunteer       558
    ## 10 advanced solid tumors   521
    ## # ℹ 30 more rows

We wanted to investigate whether there were certain countries with more
missing MeSH terms for conditions. It appears the ocuntries with the
highest percent of missing MeSH terms are in the Middle East.

``` r
study_groups |>
  filter(study_type == "INTERVENTIONAL") |>
  inner_join(filter(countries), by = "nct_id") |>
  group_by(country) |>
  summarize(
    n_trials     = n(),
    pct_no_mesh  = mean(group == "No MeSH")
  ) |>
  filter(n_trials >= 100) |>
  arrange(desc(pct_no_mesh))
```

    ## # A tibble: 106 × 3
    ##    country    n_trials pct_no_mesh
    ##    <chr>         <int>       <dbl>
    ##  1 Egypt         13774       0.379
    ##  2 Syria           281       0.320
    ##  3 Iraq            429       0.296
    ##  4 Jordan          553       0.280
    ##  5 Bangladesh      733       0.270
    ##  6 Rwanda          130       0.269
    ##  7 Ethiopia        267       0.262
    ##  8 China         42429       0.261
    ##  9 Pakistan       5793       0.247
    ## 10 Lebanon         658       0.243
    ## # ℹ 96 more rows

``` r
study_groups |>
  filter(study_type == "INTERVENTIONAL",
         start_date_type == "ACTUAL",
         !is.na(start_date))
```

    ## # A tibble: 271,846 × 18
    ##    nct_id      brief_title          official_title acronym study_type start_date
    ##    <chr>       <chr>                <chr>          <chr>   <chr>      <date>    
    ##  1 NCT03124173 Cognitive Control o… Cognitive Con… LANGUA… INTERVENT… 2015-01-14
    ##  2 NCT06053320 Spatial-Motor Strok… Novel Spatial… <NA>    INTERVENT… 2023-04-11
    ##  3 NCT06464237 Metacognitive Inter… Metacognitive… Mio     INTERVENT… 2025-08-01
    ##  4 NCT07250165 Combined Therapy of… Prospective, … <NA>    INTERVENT… 2024-08-10
    ##  5 NCT06029426 Study to Evaluate t… A Phase 3, Ra… <NA>    INTERVENT… 2023-09-20
    ##  6 NCT05505630 The Effects of Seda… The Effects o… SED-TU… INTERVENT… 2023-01-31
    ##  7 NCT04835142 Comparison of A140 … Phase III Cli… <NA>    INTERVENT… 2020-12-31
    ##  8 NCT07643714 Effects of ACT-Base… The Effect of… <NA>    INTERVENT… 2025-10-06
    ##  9 NCT06745479 Decision-Making in … Promoting Ada… <NA>    INTERVENT… 2025-01-07
    ## 10 NCT02527031 A Comparative Study… A Comparative… APACAR2 INTERVENT… 2016-03-29
    ## # ℹ 271,836 more rows
    ## # ℹ 12 more variables: start_date_type <chr>, completion_date <date>,
    ## #   overall_status <chr>, why_stopped <chr>, phase <chr>, enrollment <dbl>,
    ## #   enrollment_type <chr>, source <chr>, source_class <chr>, group <chr>,
    ## #   healthy_volunteers <chr>, primary_purpose <chr>

Next steps:

- download and merge the WHO DALY data
