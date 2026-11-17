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
  filter(!is.na(start_date), start_date_type == "ACTUAL")

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
| Classified into at least one category       | 254,846 | 70%   |
| No MeSH condition terms at all              | 79,987  | 22%   |
| Have MeSH terms, but none in our categories | 26,853  | 7%    |

We will use the classified mesh terms and list the “have mesh terms but
none in our categories” as “other”. For analyses based on condition, we
will exclude those with no Mesh terms, otherwise we will keep them in.

Initial analyses

``` r
# unique conditions
unique_conditions = unique(browse_conditions$mesh_term_condition)

# 50 least common mesh terms
# 6,052 total condition mesh terms (ancestor + list)
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
# mean rows (number of mesh terms) per trial
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
# looking at mesh ancestors (broader terms --> what we're more interested in)
# 1,854 ancestors 
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

### We’re conflicted on whether or not to include trials that are observational in our analyses.

The argument against it is that all clinical trials are registered on
clinicaltrials.gov, whereas not all observational studies are. By
including observational studies, we may have a lot of bias

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

    ## # A tibble: 79,987 × 4
    ##    nct_id      brief_title                                      enrollment group
    ##    <chr>       <chr>                                                 <dbl> <chr>
    ##  1 NCT03124173 Cognitive Control of Language                            23 No M…
    ##  2 NCT06401486 DOuble-Lumen Intubation With VIdeolaryngoscopy          916 No M…
    ##  3 NCT06632483 An Observational Study to Learn More About How …        604 No M…
    ##  4 NCT04910815 Feasibility of a New Diagnostic Device to Asses…        110 No M…
    ##  5 NCT06763536 Oral Nicotine Pouches Versus Nicotine Replaceme…       1000 No M…
    ##  6 NCT06955403 Acuvue Oaysis Contact Lens Wearers Being Refit …         40 No M…
    ##  7 NCT06479330 Evaluating the Impact of Incremental Doses of a…         60 No M…
    ##  8 NCT05409859 Adverse Events Associated With Manual Therapies…         74 No M…
    ##  9 NCT05770024 Defining, Evaluating, and Sharing Methodologies…       1500 No M…
    ## 10 NCT07631949 Osteopathy in Postpartum Lumbopelvic Pain                45 No M…
    ## # ℹ 79,977 more rows

Exploring studies with missing MeSH terms for conditions

``` r
eligibilities = read_aact("eligibilities", col_select = c(nct_id, healthy_volunteers))
designs       = read_aact("designs", col_select = c(nct_id, primary_purpose))

write_rds(eligibilities, "data/eligibilities.rds", compress = "gz")
write_rds(designs, "data/designs.rds", compress = "gz")

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

    ## # A tibble: 6 × 4
    ##   group      study_type          n   pct
    ##   <chr>      <chr>           <int> <dbl>
    ## 1 Classified INTERVENTIONAL 188554 0.740
    ## 2 Classified OBSERVATIONAL   66292 0.260
    ## 3 No MeSH    INTERVENTIONAL  61743 0.772
    ## 4 No MeSH    OBSERVATIONAL   18244 0.228
    ## 5 Other MeSH INTERVENTIONAL  21549 0.802
    ## 6 Other MeSH OBSERVATIONAL    5304 0.198

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
    ##  1 healthy                3996
    ##  2 healthy volunteers     1453
    ##  3 anesthesia              555
    ##  4 healthy participants    525
    ##  5 aging                   427
    ##  6 healthy subjects        409
    ##  7 advanced solid tumor    371
    ##  8 advanced solid tumors   364
    ##  9 solid tumor             356
    ## 10 hiv                     352
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

    ## # A tibble: 91 × 3
    ##    country              n_trials pct_no_mesh
    ##    <chr>                   <int>       <dbl>
    ##  1 Egypt                   11073       0.383
    ##  2 Jordan                    460       0.298
    ##  3 Iraq                      356       0.292
    ##  4 Syria                     241       0.290
    ##  5 Bangladesh                483       0.286
    ##  6 Ethiopia                  193       0.280
    ##  7 China                   29648       0.279
    ##  8 Lebanon                   468       0.252
    ##  9 Pakistan                 4938       0.248
    ## 10 United Arab Emirates      270       0.241
    ## # ℹ 81 more rows

# Importing WHO DALY Data (2010, 2015, 2019, 2023)

``` r
library(readxl)

source("source/read_daly.R")

daly_df = bind_rows(
  read_daly(2010),
  read_daly(2015),
  read_daly(2019),
  read_daly(2023)
)

# how should we merge? --> thinking by country, year, and cause (need to match causes to mesh terms)
```

# Trials_df

``` r
trials_df = studies |> 
  mutate(start_year = lubridate::year(start_date)) |> 
  distinct(nct_id, .keep_all=TRUE)

trials_df |> 
  count(study_type, sort=T)
```

    ## # A tibble: 2 × 2
    ##   study_type          n
    ##   <chr>           <int>
    ## 1 INTERVENTIONAL 271846
    ## 2 OBSERVATIONAL   89840

``` r
unique(countries$country)
```

    ##   [1] "China"                               
    ##   [2] "South Korea"                         
    ##   [3] "Taiwan"                              
    ##   [4] "Thailand"                            
    ##   [5] "Russia"                              
    ##   [6] "United States"                       
    ##   [7] "Bangladesh"                          
    ##   [8] "Norway"                              
    ##   [9] "Canada"                              
    ##  [10] "Turkey (Türkiye)"                    
    ##  [11] "Egypt"                               
    ##  [12] "India"                               
    ##  [13] "Italy"                               
    ##  [14] "United Kingdom"                      
    ##  [15] "Finland"                             
    ##  [16] "Brazil"                              
    ##  [17] "Greece"                              
    ##  [18] "Tanzania"                            
    ##  [19] "Germany"                             
    ##  [20] "Poland"                              
    ##  [21] "Portugal"                            
    ##  [22] "Sweden"                              
    ##  [23] "Argentina"                           
    ##  [24] "Belgium"                             
    ##  [25] "France"                              
    ##  [26] "Hong Kong"                           
    ##  [27] "Spain"                               
    ##  [28] "Bulgaria"                            
    ##  [29] "Czechia"                             
    ##  [30] "Hungary"                             
    ##  [31] "Lithuania"                           
    ##  [32] "Malaysia"                            
    ##  [33] "Mexico"                              
    ##  [34] "Netherlands"                         
    ##  [35] "Romania"                             
    ##  [36] "Slovakia"                            
    ##  [37] "Switzerland"                         
    ##  [38] "Saudi Arabia"                        
    ##  [39] "Austria"                             
    ##  [40] "Australia"                           
    ##  [41] "South Africa"                        
    ##  [42] "Ukraine"                             
    ##  [43] "Chile"                               
    ##  [44] "Indonesia"                           
    ##  [45] "Israel"                              
    ##  [46] "Japan"                               
    ##  [47] "Pakistan"                            
    ##  [48] "Peru"                                
    ##  [49] "Philippines"                         
    ##  [50] "Singapore"                           
    ##  [51] "Iraq"                                
    ##  [52] "Belarus"                             
    ##  [53] "Serbia"                              
    ##  [54] "Denmark"                             
    ##  [55] "Puerto Rico"                         
    ##  [56] "Croatia"                             
    ##  [57] "Iran"                                
    ##  [58] "Kenya"                               
    ##  [59] "Nigeria"                             
    ##  [60] "Uganda"                              
    ##  [61] "Guinea"                              
    ##  [62] "Bosnia and Herzegovina"              
    ##  [63] "Costa Rica"                          
    ##  [64] "Latvia"                              
    ##  [65] "New Zealand"                         
    ##  [66] "Vietnam"                             
    ##  [67] "Kazakhstan"                          
    ##  [68] "North Macedonia"                     
    ##  [69] "Georgia"                             
    ##  [70] "Moldova"                             
    ##  [71] "Cambodia"                            
    ##  [72] "Cameroon"                            
    ##  [73] "Ethiopia"                            
    ##  [74] "Cyprus"                              
    ##  [75] "Colombia"                            
    ##  [76] "Dominican Republic"                  
    ##  [77] "Guatemala"                           
    ##  [78] "Panama"                              
    ##  [79] "Ireland"                             
    ##  [80] "Luxembourg"                          
    ##  [81] "Armenia"                             
    ##  [82] "Paraguay"                            
    ##  [83] "Zambia"                              
    ##  [84] "Sri Lanka"                           
    ##  [85] "Sudan"                               
    ##  [86] "Mali"                                
    ##  [87] "Lebanon"                             
    ##  [88] "Jamaica"                             
    ##  [89] "Nepal"                               
    ##  [90] "Monaco"                              
    ##  [91] "Jordan"                              
    ##  [92] "Morocco"                             
    ##  [93] "Democratic Republic of the Congo"    
    ##  [94] "Papua New Guinea"                    
    ##  [95] "Iceland"                             
    ##  [96] "Tunisia"                             
    ##  [97] "Fiji"                                
    ##  [98] "Laos"                                
    ##  [99] "Mozambique"                          
    ## [100] "Senegal"                             
    ## [101] "United Arab Emirates"                
    ## [102] "Estonia"                             
    ## [103] "Venezuela"                           
    ## [104] "Benin"                               
    ## [105] "Serbia and Montenegro"               
    ## [106] "Slovenia"                            
    ## [107] "Mongolia"                            
    ## [108] "Martinique"                          
    ## [109] "Mauritius"                           
    ## [110] "Malawi"                              
    ## [111] "Burkina Faso"                        
    ## [112] "Zimbabwe"                            
    ## [113] "Honduras"                            
    ## [114] "Trinidad and Tobago"                 
    ## [115] "Ecuador"                             
    ## [116] "Netherlands Antilles"                
    ## [117] "Kyrgyzstan"                          
    ## [118] "Kosovo"                              
    ## [119] "Syria"                               
    ## [120] "Uruguay"                             
    ## [121] "The Gambia"                          
    ## [122] "Yemen"                               
    ## [123] "Guadeloupe"                          
    ## [124] "Reunion"                             
    ## [125] "Botswana"                            
    ## [126] "Albania"                             
    ## [127] "Uzbekistan"                          
    ## [128] "Republic of the Congo"               
    ## [129] "Ghana"                               
    ## [130] "Qatar"                               
    ## [131] "Côte d’Ivoire"                       
    ## [132] "Bolivia"                             
    ## [133] "Rwanda"                              
    ## [134] "French Guiana"                       
    ## [135] "Palestinian Territories"             
    ## [136] "Sierra Leone"                        
    ## [137] NA                                    
    ## [138] "Bahrain"                             
    ## [139] "Kuwait"                              
    ## [140] "Oman"                                
    ## [141] "Solomon Islands"                     
    ## [142] "Burma"                               
    ## [143] "Haiti"                               
    ## [144] "American Samoa"                      
    ## [145] "Algeria"                             
    ## [146] "Niger"                               
    ## [147] "Afghanistan"                         
    ## [148] "Djibouti"                            
    ## [149] "Cuba"                                
    ## [150] "Gabon"                               
    ## [151] "Nicaragua"                           
    ## [152] "Malta"                               
    ## [153] "South Sudan"                         
    ## [154] "Eswatini"                            
    ## [155] "El Salvador"                         
    ## [156] "Libya"                               
    ## [157] "Madagascar"                          
    ## [158] "Greenland"                           
    ## [159] "Central African Republic"            
    ## [160] "Guinea-Bissau"                       
    ## [161] "Guam"                                
    ## [162] "North Korea"                         
    ## [163] "Holy See"                            
    ## [164] "Seychelles"                          
    ## [165] "Belize"                              
    ## [166] "The Bahamas"                         
    ## [167] "Maldives"                            
    ## [168] "Angola"                              
    ## [169] "Azerbaijan"                          
    ## [170] "Togo"                                
    ## [171] "Montenegro"                          
    ## [172] "Bhutan"                              
    ## [173] "Namibia"                             
    ## [174] "Liechtenstein"                       
    ## [175] "Antigua and Barbuda"                 
    ## [176] "Lesotho"                             
    ## [177] "Macau"                               
    ## [178] "Andorra"                             
    ## [179] "Anguilla"                            
    ## [180] "Chad"                                
    ## [181] "Brunei"                              
    ## [182] "Somalia"                             
    ## [183] "Saint Martin"                        
    ## [184] "Cabo Verde"                          
    ## [185] "Samoa"                               
    ## [186] "French Polynesia"                    
    ## [187] "Liberia"                             
    ## [188] "Grenada"                             
    ## [189] "Burundi"                             
    ## [190] "Barbados"                            
    ## [191] "Virgin Islands"                      
    ## [192] "Suriname"                            
    ## [193] "Mauritania"                          
    ## [194] "Cayman Islands"                      
    ## [195] "New Caledonia"                       
    ## [196] "San Marino"                          
    ## [197] "Comoros"                             
    ## [198] "Federal Republic of Yugoslavia"      
    ## [199] "Jersey"                              
    ## [200] "Mayotte"                             
    ## [201] "Saint Kitts and Nevis"               
    ## [202] "Northern Mariana Islands"            
    ## [203] "Guyana"                              
    ## [204] "Faroe Islands"                       
    ## [205] "Bermuda"                             
    ## [206] "Aland Islands"                       
    ## [207] "Timor-Leste"                         
    ## [208] "Antarctica"                          
    ## [209] "Bonaire, Saint Eustatius and Saba "  
    ## [210] "Tajikistan"                          
    ## [211] "Niue"                                
    ## [212] "Saint Lucia"                         
    ## [213] "Vanuatu"                             
    ## [214] "Montserrat"                          
    ## [215] "Dominica"                            
    ## [216] "Saint Vincent and the Grenadines"    
    ## [217] "United States Minor Outlying Islands"
    ## [218] "Christmas Island"                    
    ## [219] "Eritrea"                             
    ## [220] "French Southern and Antarctic Lands" 
    ## [221] "Equatorial Guinea"                   
    ## [222] "Micronesia"                          
    ## [223] "Kiribati"                            
    ## [224] "Gibraltar"                           
    ## [225] "Aruba"                               
    ## [226] "Turkmenistan"                        
    ## [227] "Curacao"

``` r
unique(daly_df$country)
```

    ##   [1] "Afghanistan"                                             
    ##   [2] "Albania"                                                 
    ##   [3] "Algeria"                                                 
    ##   [4] "Angola"                                                  
    ##   [5] "Antigua and Barbuda"                                     
    ##   [6] "Argentina"                                               
    ##   [7] "Armenia"                                                 
    ##   [8] "Australia"                                               
    ##   [9] "Austria"                                                 
    ##  [10] "Azerbaijan"                                              
    ##  [11] "Bahamas"                                                 
    ##  [12] "Bahrain"                                                 
    ##  [13] "Bangladesh"                                              
    ##  [14] "Barbados"                                                
    ##  [15] "Belarus"                                                 
    ##  [16] "Belgium"                                                 
    ##  [17] "Belize"                                                  
    ##  [18] "Benin"                                                   
    ##  [19] "Bhutan"                                                  
    ##  [20] "Bolivia (Plurinational State of)"                        
    ##  [21] "Bosnia and Herzegovina"                                  
    ##  [22] "Botswana"                                                
    ##  [23] "Brazil"                                                  
    ##  [24] "Brunei Darussalam"                                       
    ##  [25] "Bulgaria"                                                
    ##  [26] "Burkina Faso"                                            
    ##  [27] "Burundi"                                                 
    ##  [28] "Cabo Verde"                                              
    ##  [29] "Cambodia"                                                
    ##  [30] "Cameroon"                                                
    ##  [31] "Canada"                                                  
    ##  [32] "Central African Republic"                                
    ##  [33] "Chad"                                                    
    ##  [34] "Chile"                                                   
    ##  [35] "China"                                                   
    ##  [36] "Colombia"                                                
    ##  [37] "Comoros"                                                 
    ##  [38] "Congo"                                                   
    ##  [39] "Costa Rica"                                              
    ##  [40] "Côte d'Ivoire"                                           
    ##  [41] "Croatia"                                                 
    ##  [42] "Cuba"                                                    
    ##  [43] "Cyprus"                                                  
    ##  [44] "Czechia"                                                 
    ##  [45] "Democratic People's Republic of Korea"                   
    ##  [46] "Democratic Republic of the Congo"                        
    ##  [47] "Denmark"                                                 
    ##  [48] "Djibouti"                                                
    ##  [49] "Dominican Republic"                                      
    ##  [50] "Ecuador"                                                 
    ##  [51] "Egypt"                                                   
    ##  [52] "El Salvador"                                             
    ##  [53] "Equatorial Guinea"                                       
    ##  [54] "Eritrea"                                                 
    ##  [55] "Estonia"                                                 
    ##  [56] "Eswatini"                                                
    ##  [57] "Ethiopia"                                                
    ##  [58] "Fiji"                                                    
    ##  [59] "Finland"                                                 
    ##  [60] "France"                                                  
    ##  [61] "Gabon"                                                   
    ##  [62] "Gambia"                                                  
    ##  [63] "Georgia"                                                 
    ##  [64] "Germany"                                                 
    ##  [65] "Ghana"                                                   
    ##  [66] "Greece"                                                  
    ##  [67] "Grenada"                                                 
    ##  [68] "Guatemala"                                               
    ##  [69] "Guinea"                                                  
    ##  [70] "Guinea-Bissau"                                           
    ##  [71] "Guyana"                                                  
    ##  [72] "Haiti"                                                   
    ##  [73] "Honduras"                                                
    ##  [74] "Hungary"                                                 
    ##  [75] "Iceland"                                                 
    ##  [76] "India"                                                   
    ##  [77] "Indonesia"                                               
    ##  [78] "Iran (Islamic Republic of)"                              
    ##  [79] "Iraq"                                                    
    ##  [80] "Ireland"                                                 
    ##  [81] "Israel"                                                  
    ##  [82] "Italy"                                                   
    ##  [83] "Jamaica"                                                 
    ##  [84] "Japan"                                                   
    ##  [85] "Jordan"                                                  
    ##  [86] "Kazakhstan"                                              
    ##  [87] "Kenya"                                                   
    ##  [88] "Kiribati"                                                
    ##  [89] "Kuwait"                                                  
    ##  [90] "Kyrgyzstan"                                              
    ##  [91] "Lao People's Democratic Republic"                        
    ##  [92] "Latvia"                                                  
    ##  [93] "Lebanon"                                                 
    ##  [94] "Lesotho"                                                 
    ##  [95] "Liberia"                                                 
    ##  [96] "Libya"                                                   
    ##  [97] "Lithuania"                                               
    ##  [98] "Luxembourg"                                              
    ##  [99] "Madagascar"                                              
    ## [100] "Malawi"                                                  
    ## [101] "Malaysia"                                                
    ## [102] "Maldives"                                                
    ## [103] "Mali"                                                    
    ## [104] "Malta"                                                   
    ## [105] "Mauritania"                                              
    ## [106] "Mauritius"                                               
    ## [107] "Mexico"                                                  
    ## [108] "Micronesia (Federated States of)"                        
    ## [109] "Mongolia"                                                
    ## [110] "Montenegro"                                              
    ## [111] "Morocco"                                                 
    ## [112] "Mozambique"                                              
    ## [113] "Myanmar"                                                 
    ## [114] "Namibia"                                                 
    ## [115] "Nepal"                                                   
    ## [116] "Netherlands"                                             
    ## [117] "New Zealand"                                             
    ## [118] "Nicaragua"                                               
    ## [119] "Niger"                                                   
    ## [120] "Nigeria"                                                 
    ## [121] "North Macedonia"                                         
    ## [122] "Norway"                                                  
    ## [123] "occupied Palestinian territory, including east Jerusalem"
    ## [124] "Oman"                                                    
    ## [125] "Pakistan"                                                
    ## [126] "Panama"                                                  
    ## [127] "Papua New Guinea"                                        
    ## [128] "Paraguay"                                                
    ## [129] "Peru"                                                    
    ## [130] "Philippines"                                             
    ## [131] "Poland"                                                  
    ## [132] "Portugal"                                                
    ## [133] "Puerto Rico"                                             
    ## [134] "Qatar"                                                   
    ## [135] "Republic of Korea"                                       
    ## [136] "Republic of Moldova"                                     
    ## [137] "Romania"                                                 
    ## [138] "Russian Federation"                                      
    ## [139] "Rwanda"                                                  
    ## [140] "Saint Lucia"                                             
    ## [141] "Saint Vincent and the Grenadines"                        
    ## [142] "Samoa"                                                   
    ## [143] "Sao Tome and Principe"                                   
    ## [144] "Saudi Arabia"                                            
    ## [145] "Senegal"                                                 
    ## [146] "Serbia"                                                  
    ## [147] "Seychelles"                                              
    ## [148] "Sierra Leone"                                            
    ## [149] "Singapore"                                               
    ## [150] "Slovakia"                                                
    ## [151] "Slovenia"                                                
    ## [152] "Solomon Islands"                                         
    ## [153] "Somalia"                                                 
    ## [154] "South Africa"                                            
    ## [155] "South Sudan"                                             
    ## [156] "Spain"                                                   
    ## [157] "Sri Lanka"                                               
    ## [158] "Sudan"                                                   
    ## [159] "Suriname"                                                
    ## [160] "Sweden"                                                  
    ## [161] "Switzerland"                                             
    ## [162] "Syrian Arab Republic"                                    
    ## [163] "Tajikistan"                                              
    ## [164] "Thailand"                                                
    ## [165] "Timor-Leste"                                             
    ## [166] "Togo"                                                    
    ## [167] "Tonga"                                                   
    ## [168] "Trinidad and Tobago"                                     
    ## [169] "Tunisia"                                                 
    ## [170] "Türkiye"                                                 
    ## [171] "Turkmenistan"                                            
    ## [172] "Uganda"                                                  
    ## [173] "Ukraine"                                                 
    ## [174] "United Arab Emirates"                                    
    ## [175] "United Kingdom"                                          
    ## [176] "United Republic of Tanzania"                             
    ## [177] "United States of America"                                
    ## [178] "Uruguay"                                                 
    ## [179] "Uzbekistan"                                              
    ## [180] "Vanuatu"                                                 
    ## [181] "Venezuela (Bolivarian Republic of)"                      
    ## [182] "Viet Nam"                                                
    ## [183] "Yemen"                                                   
    ## [184] "Zambia"                                                  
    ## [185] "Zimbabwe"

``` r
setdiff(countries$country, daly_df$country)
```

    ##  [1] "South Korea"                         
    ##  [2] "Taiwan"                              
    ##  [3] "Russia"                              
    ##  [4] "United States"                       
    ##  [5] "Turkey (Türkiye)"                    
    ##  [6] "Tanzania"                            
    ##  [7] "Hong Kong"                           
    ##  [8] "Iran"                                
    ##  [9] "Vietnam"                             
    ## [10] "Moldova"                             
    ## [11] "Monaco"                              
    ## [12] "Laos"                                
    ## [13] "Venezuela"                           
    ## [14] "Serbia and Montenegro"               
    ## [15] "Martinique"                          
    ## [16] "Netherlands Antilles"                
    ## [17] "Kosovo"                              
    ## [18] "Syria"                               
    ## [19] "The Gambia"                          
    ## [20] "Guadeloupe"                          
    ## [21] "Reunion"                             
    ## [22] "Republic of the Congo"               
    ## [23] "Côte d’Ivoire"                       
    ## [24] "Bolivia"                             
    ## [25] "French Guiana"                       
    ## [26] "Palestinian Territories"             
    ## [27] NA                                    
    ## [28] "Burma"                               
    ## [29] "American Samoa"                      
    ## [30] "Greenland"                           
    ## [31] "Guam"                                
    ## [32] "North Korea"                         
    ## [33] "Holy See"                            
    ## [34] "The Bahamas"                         
    ## [35] "Liechtenstein"                       
    ## [36] "Macau"                               
    ## [37] "Andorra"                             
    ## [38] "Anguilla"                            
    ## [39] "Brunei"                              
    ## [40] "Saint Martin"                        
    ## [41] "French Polynesia"                    
    ## [42] "Virgin Islands"                      
    ## [43] "Cayman Islands"                      
    ## [44] "New Caledonia"                       
    ## [45] "San Marino"                          
    ## [46] "Federal Republic of Yugoslavia"      
    ## [47] "Jersey"                              
    ## [48] "Mayotte"                             
    ## [49] "Saint Kitts and Nevis"               
    ## [50] "Northern Mariana Islands"            
    ## [51] "Faroe Islands"                       
    ## [52] "Bermuda"                             
    ## [53] "Aland Islands"                       
    ## [54] "Antarctica"                          
    ## [55] "Bonaire, Saint Eustatius and Saba "  
    ## [56] "Niue"                                
    ## [57] "Montserrat"                          
    ## [58] "Dominica"                            
    ## [59] "United States Minor Outlying Islands"
    ## [60] "Christmas Island"                    
    ## [61] "French Southern and Antarctic Lands" 
    ## [62] "Micronesia"                          
    ## [63] "Gibraltar"                           
    ## [64] "Aruba"                               
    ## [65] "Curacao"

``` r
setdiff(daly_df$country, countries$country)
```

    ##  [1] "Bahamas"                                                 
    ##  [2] "Bolivia (Plurinational State of)"                        
    ##  [3] "Brunei Darussalam"                                       
    ##  [4] "Congo"                                                   
    ##  [5] "Côte d'Ivoire"                                           
    ##  [6] "Democratic People's Republic of Korea"                   
    ##  [7] "Gambia"                                                  
    ##  [8] "Iran (Islamic Republic of)"                              
    ##  [9] "Lao People's Democratic Republic"                        
    ## [10] "Micronesia (Federated States of)"                        
    ## [11] "Myanmar"                                                 
    ## [12] "occupied Palestinian territory, including east Jerusalem"
    ## [13] "Republic of Korea"                                       
    ## [14] "Republic of Moldova"                                     
    ## [15] "Russian Federation"                                      
    ## [16] "Sao Tome and Principe"                                   
    ## [17] "Syrian Arab Republic"                                    
    ## [18] "Tonga"                                                   
    ## [19] "Türkiye"                                                 
    ## [20] "United Republic of Tanzania"                             
    ## [21] "United States of America"                                
    ## [22] "Venezuela (Bolivarian Republic of)"                      
    ## [23] "Viet Nam"

``` r
# make countries the same across data sets
```

``` r
library(countrycode)

countries <- countries |>
  mutate(iso3 = countrycode(country, origin = "country.name", destination = "iso3c"))
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `iso3 = countrycode(country, origin = "country.name",
    ##   destination = "iso3c")`.
    ## Caused by warning:
    ## ! Some values were not matched unambiguously: Federal Republic of Yugoslavia, Kosovo, Micronesia, Netherlands Antilles, Saint Martin, Serbia and Montenegro, Virgin Islands
    ## To fix unmatched values, please use the `custom_match` argument. If you think the default matching rules should be improved, please file an issue at https://github.com/vincentarelbundock/countrycode/issues

``` r
no_iso = countries |> 
  filter(is.na(iso3)) |> 
  distinct(country)

# we will remove these countries from the country-level analyses, since we can't match them to DALY data
```
