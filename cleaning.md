Data
================

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
studies = read_rds("data/studies.rds")
countries = read_rds("data/countries.rds")
sponsors = read_rds("data/sponsors.rds")
conditions = read_rds("data/conditions.rds")
browse_conditions = read_rds("data/browse_conditions.rds")
browse_interventions = read_rds("data/browse_interventions.rds")
baseline_counts = read_rds("data/baseline_counts.rds")
eligibilities = read_rds("data/eligibilities.rds")
designs = read_rds("data/designs.rds")
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
    ##    mesh_term_condition                    n
    ##    <chr>                              <int>
    ##  1 neoplasms                          64268
    ##  2 nervous system diseases            48442
    ##  3 cardiovascular diseases            35394
    ##  4 mental disorders                   33623
    ##  5 urogenital diseases                31089
    ##  6 nutritional and metabolic diseases 30841
    ##  7 respiratory tract diseases         30167
    ##  8 digestive system diseases          27250
    ##  9 immune system diseases             26448
    ## 10 infections                         25123
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
| Classified into at least one category       | 267,759 | 69%   |
| No MeSH condition terms at all              | 89,438  | 23%   |
| Have MeSH terms, but none in our categories | 30,332  | 8%    |

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
    ##    mesh_term_condition             n
    ##    <chr>                       <int>
    ##  1 tay-sachs disease, juvenile     1
    ##  2 teratogenesis                   1
    ##  3 teratoid tumor, atypical        1
    ##  4 thrombophilia, hereditary       1
    ##  5 thumb deformity                 1
    ##  6 thyroiditis, chronic            1
    ##  7 time pressure                   1
    ##  8 timothy syndrome                1
    ##  9 tinea barbae                    1
    ## 10 tongue habits                   1
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
    ## 1  8.88

``` r
# looking at mesh ancestors (broader terms --> what we're more interested in)
# 1,854 ancestors 
browse_conditions |> 
  filter(mesh_type_condition == "mesh-ancestor") |>
  count(mesh_term_condition, sort = T) |> 
  head(30)
```

    ## # A tibble: 30 × 2
    ##    mesh_term_condition                             n
    ##    <chr>                                       <int>
    ##  1 pathological conditions, signs and symptoms 90254
    ##  2 neoplasms                                   58106
    ##  3 signs and symptoms                          51862
    ##  4 nervous system diseases                     47954
    ##  5 pathologic processes                        44668
    ##  6 neoplasms by site                           43151
    ##  7 cardiovascular diseases                     32894
    ##  8 mental disorders                            32482
    ##  9 neurologic manifestations                   31864
    ## 10 urogenital diseases                         31063
    ## # ℹ 20 more rows

``` r
unique(sponsors$agency_class)
```

    ## [1] "INDUSTRY"  "OTHER"     "OTHER_GOV" "NIH"       "UNKNOWN"   "FED"      
    ## [7] "INDIV"     "NETWORK"   "AMBIG"

### QUESTION: We’re conflicted on whether or not to include trials that are observational in our analyses.

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

    ## # A tibble: 89,438 × 4
    ##    nct_id      brief_title                                      enrollment group
    ##    <chr>       <chr>                                                 <dbl> <chr>
    ##  1 NCT03124173 Cognitive Control of Language                            23 No M…
    ##  2 NCT06401486 DOuble-Lumen Intubation With VIdeolaryngoscopy          916 No M…
    ##  3 NCT04910815 Feasibility of a New Diagnostic Device to Asses…        110 No M…
    ##  4 NCT06763536 Oral Nicotine Pouches Versus Nicotine Replaceme…       1000 No M…
    ##  5 NCT06955403 Acuvue Oaysis Contact Lens Wearers Being Refit …         40 No M…
    ##  6 NCT06479330 Evaluating the Impact of Incremental Doses of a…         60 No M…
    ##  7 NCT02556541 Ultrasound-guided Peripheral Vascular Access in…         50 No M…
    ##  8 NCT07631949 Osteopathy in Postpartum Lumbopelvic Pain                45 No M…
    ##  9 NCT06542549 Efficacy and Safety of AK104 Combined With Chem…        100 No M…
    ## 10 NCT02522832 Discovering New Treatments for Asthma and COPD.…         17 No M…
    ## # ℹ 89,428 more rows

Exploring studies with missing MeSH terms for conditions

``` r
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

    ## # A tibble: 3 × 4
    ##   group      study_type          n   pct
    ##   <chr>      <chr>           <int> <dbl>
    ## 1 Classified INTERVENTIONAL 267759     1
    ## 2 No MeSH    INTERVENTIONAL  89438     1
    ## 3 Other MeSH INTERVENTIONAL  30332     1

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
    ##  1 healthy                5877
    ##  2 healthy volunteers     1977
    ##  3 healthy subjects        606
    ##  4 anesthesia              578
    ##  5 healthy participants    576
    ##  6 hiv                     526
    ##  7 healthy volunteer       485
    ##  8 aging                   471
    ##  9 advanced solid tumors   464
    ## 10 advanced solid tumor    425
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

    ## # A tibble: 98 × 3
    ##    country    n_trials pct_no_mesh
    ##    <chr>         <int>       <dbl>
    ##  1 Egypt         13556       0.384
    ##  2 Syria           275       0.327
    ##  3 Iraq            425       0.299
    ##  4 Jordan          526       0.289
    ##  5 Bangladesh      657       0.288
    ##  6 Rwanda          117       0.282
    ##  7 China         40792       0.267
    ##  8 Ethiopia        248       0.266
    ##  9 Lebanon         577       0.260
    ## 10 Pakistan       5658       0.250
    ## # ℹ 88 more rows

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

# QUESTION: how should we merge? --> thinking by country, year, and cause (need to match causes to mesh terms)
```

# Trials_df

``` r
trials_df = studies |> 
  mutate(start_year = lubridate::year(start_date)) |> 
  distinct(nct_id, .keep_all=TRUE)

trials_df |> 
  count(study_type, sort=T)
```

    ## # A tibble: 1 × 2
    ##   study_type          n
    ##   <chr>           <int>
    ## 1 INTERVENTIONAL 387529

``` r
unique(countries$country)
```

    ##   [1] "China"                               
    ##   [2] "South Korea"                         
    ##   [3] "Taiwan"                              
    ##   [4] "Thailand"                            
    ##   [5] "Russia"                              
    ##   [6] "Bangladesh"                          
    ##   [7] "Canada"                              
    ##   [8] "United States"                       
    ##   [9] "Egypt"                               
    ##  [10] "India"                               
    ##  [11] "Italy"                               
    ##  [12] "Finland"                             
    ##  [13] "Greece"                              
    ##  [14] "Germany"                             
    ##  [15] "Poland"                              
    ##  [16] "Portugal"                            
    ##  [17] "United Kingdom"                      
    ##  [18] "Sweden"                              
    ##  [19] "Argentina"                           
    ##  [20] "Belgium"                             
    ##  [21] "France"                              
    ##  [22] "Brazil"                              
    ##  [23] "Hong Kong"                           
    ##  [24] "Turkey (Türkiye)"                    
    ##  [25] "Spain"                               
    ##  [26] "Bulgaria"                            
    ##  [27] "Czechia"                             
    ##  [28] "Hungary"                             
    ##  [29] "Lithuania"                           
    ##  [30] "Malaysia"                            
    ##  [31] "Mexico"                              
    ##  [32] "Netherlands"                         
    ##  [33] "Romania"                             
    ##  [34] "Slovakia"                            
    ##  [35] "Saudi Arabia"                        
    ##  [36] "Australia"                           
    ##  [37] "South Africa"                        
    ##  [38] "Ukraine"                             
    ##  [39] "Singapore"                           
    ##  [40] "Peru"                                
    ##  [41] "Switzerland"                         
    ##  [42] "Iraq"                                
    ##  [43] "Belarus"                             
    ##  [44] "Serbia"                              
    ##  [45] "Israel"                              
    ##  [46] "Norway"                              
    ##  [47] "Denmark"                             
    ##  [48] "Croatia"                             
    ##  [49] "Iran"                                
    ##  [50] "Pakistan"                            
    ##  [51] "Guinea"                              
    ##  [52] "Japan"                               
    ##  [53] "Vietnam"                             
    ##  [54] "New Zealand"                         
    ##  [55] "Puerto Rico"                         
    ##  [56] "Uganda"                              
    ##  [57] "Chile"                               
    ##  [58] "Indonesia"                           
    ##  [59] "Georgia"                             
    ##  [60] "Moldova"                             
    ##  [61] "Cyprus"                              
    ##  [62] "Ethiopia"                            
    ##  [63] "Austria"                             
    ##  [64] "Nigeria"                             
    ##  [65] "Armenia"                             
    ##  [66] "Paraguay"                            
    ##  [67] "Zambia"                              
    ##  [68] "Ireland"                             
    ##  [69] "Jamaica"                             
    ##  [70] "Lebanon"                             
    ##  [71] "North Macedonia"                     
    ##  [72] "Jordan"                              
    ##  [73] "Morocco"                             
    ##  [74] "Latvia"                              
    ##  [75] "Iceland"                             
    ##  [76] "Tunisia"                             
    ##  [77] "Colombia"                            
    ##  [78] "Slovenia"                            
    ##  [79] "Kazakhstan"                          
    ##  [80] "Sri Lanka"                           
    ##  [81] "Estonia"                             
    ##  [82] "Mongolia"                            
    ##  [83] "Costa Rica"                          
    ##  [84] "Dominican Republic"                  
    ##  [85] "Guatemala"                           
    ##  [86] "Mauritius"                           
    ##  [87] "Philippines"                         
    ##  [88] "United Arab Emirates"                
    ##  [89] "Malawi"                              
    ##  [90] "Burkina Faso"                        
    ##  [91] "Kenya"                               
    ##  [92] "Zimbabwe"                            
    ##  [93] "Honduras"                            
    ##  [94] "Luxembourg"                          
    ##  [95] "Monaco"                              
    ##  [96] "Bosnia and Herzegovina"              
    ##  [97] "Tanzania"                            
    ##  [98] "Benin"                               
    ##  [99] "Netherlands Antilles"                
    ## [100] "Venezuela"                           
    ## [101] "Trinidad and Tobago"                 
    ## [102] "Cambodia"                            
    ## [103] "Kyrgyzstan"                          
    ## [104] "Sudan"                               
    ## [105] "Ecuador"                             
    ## [106] "Kosovo"                              
    ## [107] "Cameroon"                            
    ## [108] "Syria"                               
    ## [109] "The Gambia"                          
    ## [110] "Yemen"                               
    ## [111] "Panama"                              
    ## [112] "Democratic Republic of the Congo"    
    ## [113] "Nepal"                               
    ## [114] "Qatar"                               
    ## [115] "Republic of the Congo"               
    ## [116] "Ghana"                               
    ## [117] "Côte d’Ivoire"                       
    ## [118] "Bolivia"                             
    ## [119] "Rwanda"                              
    ## [120] "French Guiana"                       
    ## [121] "Uruguay"                             
    ## [122] "Sierra Leone"                        
    ## [123] "Mali"                                
    ## [124] "Senegal"                             
    ## [125] "Mozambique"                          
    ## [126] "Solomon Islands"                     
    ## [127] "Burma"                               
    ## [128] "American Samoa"                      
    ## [129] "Algeria"                             
    ## [130] "Niger"                               
    ## [131] "Oman"                                
    ## [132] "Afghanistan"                         
    ## [133] "Djibouti"                            
    ## [134] "Gabon"                               
    ## [135] "Malta"                               
    ## [136] "Martinique"                          
    ## [137] "Albania"                             
    ## [138] "Nicaragua"                           
    ## [139] "Eswatini"                            
    ## [140] "El Salvador"                         
    ## [141] "Botswana"                            
    ## [142] "Greenland"                           
    ## [143] "Central African Republic"            
    ## [144] "Guam"                                
    ## [145] "Kuwait"                              
    ## [146] "Laos"                                
    ## [147] "North Korea"                         
    ## [148] "Holy See"                            
    ## [149] "Guadeloupe"                          
    ## [150] "Papua New Guinea"                    
    ## [151] "Guinea-Bissau"                       
    ## [152] "Cuba"                                
    ## [153] "Haiti"                               
    ## [154] "Angola"                              
    ## [155] "Azerbaijan"                          
    ## [156] "Belize"                              
    ## [157] "Togo"                                
    ## [158] "Bahrain"                             
    ## [159] "Montenegro"                          
    ## [160] "South Sudan"                         
    ## [161] "Reunion"                             
    ## [162] "Bhutan"                              
    ## [163] "Liechtenstein"                       
    ## [164] "Antigua and Barbuda"                 
    ## [165] "Uzbekistan"                          
    ## [166] "Lesotho"                             
    ## [167] "Madagascar"                          
    ## [168] "Macau"                               
    ## [169] "Andorra"                             
    ## [170] "Anguilla"                            
    ## [171] "Chad"                                
    ## [172] "Fiji"                                
    ## [173] "Palestinian Territories"             
    ## [174] "Brunei"                              
    ## [175] "Somalia"                             
    ## [176] "Liberia"                             
    ## [177] "French Polynesia"                    
    ## [178] "Grenada"                             
    ## [179] "Burundi"                             
    ## [180] "Serbia and Montenegro"               
    ## [181] "Barbados"                            
    ## [182] "The Bahamas"                         
    ## [183] "Suriname"                            
    ## [184] "Cayman Islands"                      
    ## [185] "Namibia"                             
    ## [186] "San Marino"                          
    ## [187] "Jersey"                              
    ## [188] "Saint Kitts and Nevis"               
    ## [189] "New Caledonia"                       
    ## [190] "Guyana"                              
    ## [191] "Faroe Islands"                       
    ## [192] "Virgin Islands"                      
    ## [193] "Bermuda"                             
    ## [194] "Aland Islands"                       
    ## [195] "Mayotte"                             
    ## [196] "Timor-Leste"                         
    ## [197] "Antarctica"                          
    ## [198] "Bonaire, Saint Eustatius and Saba "  
    ## [199] "Tajikistan"                          
    ## [200] "Libya"                               
    ## [201] "Niue"                                
    ## [202] "Saint Lucia"                         
    ## [203] "Vanuatu"                             
    ## [204] "Montserrat"                          
    ## [205] "Northern Mariana Islands"            
    ## [206] NA                                    
    ## [207] "United States Minor Outlying Islands"
    ## [208] "Christmas Island"                    
    ## [209] "Comoros"                             
    ## [210] "Samoa"                               
    ## [211] "Eritrea"                             
    ## [212] "Equatorial Guinea"                   
    ## [213] "Micronesia"                          
    ## [214] "Kiribati"                            
    ## [215] "Turkmenistan"

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
    ##  [5] "Hong Kong"                           
    ##  [6] "Turkey (Türkiye)"                    
    ##  [7] "Iran"                                
    ##  [8] "Vietnam"                             
    ##  [9] "Moldova"                             
    ## [10] "Monaco"                              
    ## [11] "Tanzania"                            
    ## [12] "Netherlands Antilles"                
    ## [13] "Venezuela"                           
    ## [14] "Kosovo"                              
    ## [15] "Syria"                               
    ## [16] "The Gambia"                          
    ## [17] "Republic of the Congo"               
    ## [18] "Côte d’Ivoire"                       
    ## [19] "Bolivia"                             
    ## [20] "French Guiana"                       
    ## [21] "Burma"                               
    ## [22] "American Samoa"                      
    ## [23] "Martinique"                          
    ## [24] "Greenland"                           
    ## [25] "Guam"                                
    ## [26] "Laos"                                
    ## [27] "North Korea"                         
    ## [28] "Holy See"                            
    ## [29] "Guadeloupe"                          
    ## [30] "Reunion"                             
    ## [31] "Liechtenstein"                       
    ## [32] "Macau"                               
    ## [33] "Andorra"                             
    ## [34] "Anguilla"                            
    ## [35] "Palestinian Territories"             
    ## [36] "Brunei"                              
    ## [37] "French Polynesia"                    
    ## [38] "Serbia and Montenegro"               
    ## [39] "The Bahamas"                         
    ## [40] "Cayman Islands"                      
    ## [41] "San Marino"                          
    ## [42] "Jersey"                              
    ## [43] "Saint Kitts and Nevis"               
    ## [44] "New Caledonia"                       
    ## [45] "Faroe Islands"                       
    ## [46] "Virgin Islands"                      
    ## [47] "Bermuda"                             
    ## [48] "Aland Islands"                       
    ## [49] "Mayotte"                             
    ## [50] "Antarctica"                          
    ## [51] "Bonaire, Saint Eustatius and Saba "  
    ## [52] "Niue"                                
    ## [53] "Montserrat"                          
    ## [54] "Northern Mariana Islands"            
    ## [55] NA                                    
    ## [56] "United States Minor Outlying Islands"
    ## [57] "Christmas Island"                    
    ## [58] "Micronesia"

``` r
setdiff(daly_df$country, countries$country)
```

    ##  [1] "Bahamas"                                                 
    ##  [2] "Bolivia (Plurinational State of)"                        
    ##  [3] "Brunei Darussalam"                                       
    ##  [4] "Cabo Verde"                                              
    ##  [5] "Congo"                                                   
    ##  [6] "Côte d'Ivoire"                                           
    ##  [7] "Democratic People's Republic of Korea"                   
    ##  [8] "Gambia"                                                  
    ##  [9] "Iran (Islamic Republic of)"                              
    ## [10] "Lao People's Democratic Republic"                        
    ## [11] "Maldives"                                                
    ## [12] "Mauritania"                                              
    ## [13] "Micronesia (Federated States of)"                        
    ## [14] "Myanmar"                                                 
    ## [15] "occupied Palestinian territory, including east Jerusalem"
    ## [16] "Republic of Korea"                                       
    ## [17] "Republic of Moldova"                                     
    ## [18] "Russian Federation"                                      
    ## [19] "Saint Vincent and the Grenadines"                        
    ## [20] "Sao Tome and Principe"                                   
    ## [21] "Seychelles"                                              
    ## [22] "Syrian Arab Republic"                                    
    ## [23] "Tonga"                                                   
    ## [24] "Türkiye"                                                 
    ## [25] "United Republic of Tanzania"                             
    ## [26] "United States of America"                                
    ## [27] "Venezuela (Bolivarian Republic of)"                      
    ## [28] "Viet Nam"

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
    ## ! Some values were not matched unambiguously: Kosovo, Micronesia, Netherlands Antilles, Serbia and Montenegro, Virgin Islands
    ## To fix unmatched values, please use the `custom_match` argument. If you think the default matching rules should be improved, please file an issue at https://github.com/vincentarelbundock/countrycode/issues

``` r
no_iso = countries |> 
  filter(is.na(iso3)) |> 
  distinct(country)

# we will remove these countries from the country-level analyses, since we can't match them to DALY data
```

# Import IHME Data: All countries/territories, World Bank Income Levels, Year 2023 only

``` r
ihme_2023 = read_csv("data/IHME-GBD_2023_DATA-30ca978e-1.csv")
```

    ## Rows: 8815 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (7): population_group, measure, location, sex, age, cause, metric
    ## dbl (4): year, val, upper, lower
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

woop. will add more years tmrw maybeeee
