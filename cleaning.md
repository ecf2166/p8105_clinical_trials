Data
================
2026-11-11

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

studies = read_aact("studies")
countries = read_aact("countries")
sponsors = read_aact("sponsors")
conditions = read_aact("conditions")
browse_conditions = read_aact("browse_conditions")
interventions = read_aact("browse_interventions")
baseline = read_aact("baseline_counts")
```

``` r
# don't want to knit this every time, only run once
# for keeping data in in RDS format to save space (and so this knits)
write_rds(studies, "data/studies.rds", compress = "gz")
write_rds(countries, "data/countries.rds", compress = "gz")
write_rds(sponsors, "data/sponsors.rds", compress = "gz")
write_rds(conditions, "data/conditions.rds", compress = "gz")
write_rds(browse_conditions, "data/browse_conditions.rds", compress = "gz")
```

``` r
unique_conditions = unique(browse_conditions$mesh_term)

browse_conditions |> 
  count(mesh_term, sort=T) |> 
  slice_tail(n=50)
```

    ## # A tibble: 50 × 2
    ##    mesh_term                                       n
    ##    <chr>                                       <int>
    ##  1 Uniparental Disomy                              1
    ##  2 Uniparental disomy, paternal, chromosome 14     1
    ##  3 Urachal Cyst                                    1
    ##  4 Ureterocele                                     1
    ##  5 Usher syndrome, type 1C                         1
    ##  6 Usher syndrome, type 1F                         1
    ##  7 Uterine Didelphys                               1
    ##  8 Uterine Retroversion                            1
    ##  9 VATER association                               1
    ## 10 Vagina, absence of                              1
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
  filter(mesh_type == "mesh-ancestor") |>
  count(mesh_term, sort = T) |> 
  head(30)
```

    ## # A tibble: 30 × 2
    ##    mesh_term                                        n
    ##    <chr>                                        <int>
    ##  1 Pathological Conditions, Signs and Symptoms 144522
    ##  2 Neoplasms                                    96256
    ##  3 Pathologic Processes                         77370
    ##  4 Nervous System Diseases                      75975
    ##  5 Signs and Symptoms                           75907
    ##  6 Neoplasms by Site                            72772
    ##  7 Cardiovascular Diseases                      59509
    ##  8 Urogenital Diseases                          54793
    ##  9 Respiratory Tract Diseases                   52049
    ## 10 Nutritional and Metabolic Diseases           50375
    ## # ℹ 20 more rows

``` r
unique(sponsors$agency_class)
```

    ##  [1] "INDUSTRY"  "OTHER"     "OTHER_GOV" "NIH"       "UNKNOWN"   "FED"      
    ##  [7] "NETWORK"   "AMBIG"     "INDIV"     NA

``` r
studies <- read_aact("studies") 
```
