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

    ## # A tibble: 129,732 × 4
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
    ## # ℹ 129,722 more rows

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

    ## # A tibble: 10 × 4
    ##    group      study_type           n      pct
    ##    <chr>      <chr>            <int>    <dbl>
    ##  1 Classified EXPANDED_ACCESS    851 0.00194 
    ##  2 Classified INTERVENTIONAL  332902 0.757   
    ##  3 Classified OBSERVATIONAL   105769 0.241   
    ##  4 No MeSH    EXPANDED_ACCESS    193 0.00149 
    ##  5 No MeSH    INTERVENTIONAL   99605 0.768   
    ##  6 No MeSH    OBSERVATIONAL    28941 0.223   
    ##  7 No MeSH    <NA>               993 0.00765 
    ##  8 Other MeSH EXPANDED_ACCESS     32 0.000762
    ##  9 Other MeSH INTERVENTIONAL   33920 0.808   
    ## 10 Other MeSH OBSERVATIONAL     8027 0.191

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
    ##  1 healthy                8678
    ##  2 healthy volunteers     2277
    ##  3 anesthesia              878
    ##  4 healthy subjects        777
    ##  5 hiv                     756
    ##  6 pregnancy               654
    ##  7 aging                   612
    ##  8 healthy participants    599
    ##  9 healthy volunteer       558
    ## 10 advanced solid tumors   522
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

    ## # A tibble: 107 × 3
    ##    country    n_trials pct_no_mesh
    ##    <chr>         <int>       <dbl>
    ##  1 Egypt         13864       0.378
    ##  2 Syria           281       0.320
    ##  3 Iraq            429       0.296
    ##  4 Jordan          582       0.277
    ##  5 Bangladesh      738       0.268
    ##  6 China         42902       0.259
    ##  7 Rwanda          139       0.259
    ##  8 Ethiopia        272       0.257
    ##  9 Pakistan       5813       0.246
    ## 10 Lebanon         703       0.233
    ## # ℹ 97 more rows

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

    ## # A tibble: 4 × 2
    ##   study_type           n
    ##   <chr>            <int>
    ## 1 INTERVENTIONAL  466427
    ## 2 OBSERVATIONAL   142737
    ## 3 EXPANDED_ACCESS   1076
    ## 4 <NA>               993

``` r
unique(countries$country)
```

    ##   [1] "Austria"                             
    ##   [2] "Germany"                             
    ##   [3] "Australia"                           
    ##   [4] "Belgium"                             
    ##   [5] "Brazil"                              
    ##   [6] "France"                              
    ##   [7] "Israel"                              
    ##   [8] "Italy"                               
    ##   [9] "Japan"                               
    ##  [10] "South Africa"                        
    ##  [11] "South Korea"                         
    ##  [12] "Spain"                               
    ##  [13] "Taiwan"                              
    ##  [14] "United Kingdom"                      
    ##  [15] "Canada"                              
    ##  [16] "Georgia"                             
    ##  [17] "Russia"                              
    ##  [18] "Turkey (Türkiye)"                    
    ##  [19] "Romania"                             
    ##  [20] "Thailand"                            
    ##  [21] "United States"                       
    ##  [22] "Sweden"                              
    ##  [23] "China"                               
    ##  [24] "Bangladesh"                          
    ##  [25] "Norway"                              
    ##  [26] "Egypt"                               
    ##  [27] "India"                               
    ##  [28] "Finland"                             
    ##  [29] "Greece"                              
    ##  [30] "Tanzania"                            
    ##  [31] "Poland"                              
    ##  [32] "Portugal"                            
    ##  [33] "Czechia"                             
    ##  [34] "Argentina"                           
    ##  [35] "Serbia and Montenegro"               
    ##  [36] "Hong Kong"                           
    ##  [37] "Bulgaria"                            
    ##  [38] "Hungary"                             
    ##  [39] "Lithuania"                           
    ##  [40] "Malaysia"                            
    ##  [41] "Mexico"                              
    ##  [42] "Netherlands"                         
    ##  [43] "Slovakia"                            
    ##  [44] "Switzerland"                         
    ##  [45] "Saudi Arabia"                        
    ##  [46] "Ukraine"                             
    ##  [47] "Chile"                               
    ##  [48] "Indonesia"                           
    ##  [49] "Pakistan"                            
    ##  [50] "Peru"                                
    ##  [51] "Philippines"                         
    ##  [52] "Singapore"                           
    ##  [53] "Iraq"                                
    ##  [54] "Belarus"                             
    ##  [55] "Serbia"                              
    ##  [56] "Denmark"                             
    ##  [57] "Puerto Rico"                         
    ##  [58] "Croatia"                             
    ##  [59] "Iran"                                
    ##  [60] "Kenya"                               
    ##  [61] "Nigeria"                             
    ##  [62] "Uganda"                              
    ##  [63] "Guinea"                              
    ##  [64] "Moldova"                             
    ##  [65] "New Zealand"                         
    ##  [66] "Bosnia and Herzegovina"              
    ##  [67] "Costa Rica"                          
    ##  [68] "Latvia"                              
    ##  [69] "Vietnam"                             
    ##  [70] "Kazakhstan"                          
    ##  [71] "North Macedonia"                     
    ##  [72] "Cambodia"                            
    ##  [73] "Cameroon"                            
    ##  [74] "Ethiopia"                            
    ##  [75] "Cyprus"                              
    ##  [76] "Colombia"                            
    ##  [77] "Dominican Republic"                  
    ##  [78] "Guatemala"                           
    ##  [79] "Panama"                              
    ##  [80] "Ireland"                             
    ##  [81] "Luxembourg"                          
    ##  [82] "Armenia"                             
    ##  [83] "Paraguay"                            
    ##  [84] "Zambia"                              
    ##  [85] "Sri Lanka"                           
    ##  [86] "Sudan"                               
    ##  [87] "Mali"                                
    ##  [88] "Lebanon"                             
    ##  [89] "Jamaica"                             
    ##  [90] "Nepal"                               
    ##  [91] "Monaco"                              
    ##  [92] "Jordan"                              
    ##  [93] "Morocco"                             
    ##  [94] "Ecuador"                             
    ##  [95] "San Marino"                          
    ##  [96] "Democratic Republic of the Congo"    
    ##  [97] "Papua New Guinea"                    
    ##  [98] "Iceland"                             
    ##  [99] "Tunisia"                             
    ## [100] "Fiji"                                
    ## [101] "Laos"                                
    ## [102] "Mozambique"                          
    ## [103] "Senegal"                             
    ## [104] "United Arab Emirates"                
    ## [105] "Estonia"                             
    ## [106] "Venezuela"                           
    ## [107] "Benin"                               
    ## [108] "Slovenia"                            
    ## [109] "Mongolia"                            
    ## [110] "Martinique"                          
    ## [111] "Mauritius"                           
    ## [112] "Malawi"                              
    ## [113] "Burkina Faso"                        
    ## [114] "El Salvador"                         
    ## [115] "Uruguay"                             
    ## [116] "Zimbabwe"                            
    ## [117] "Honduras"                            
    ## [118] "Trinidad and Tobago"                 
    ## [119] "Netherlands Antilles"                
    ## [120] "Kyrgyzstan"                          
    ## [121] "South Georgia and the Islands"       
    ## [122] "Kosovo"                              
    ## [123] "Syria"                               
    ## [124] "The Gambia"                          
    ## [125] "Ghana"                               
    ## [126] "Yemen"                               
    ## [127] "Guadeloupe"                          
    ## [128] "Reunion"                             
    ## [129] "Botswana"                            
    ## [130] "Albania"                             
    ## [131] "Uzbekistan"                          
    ## [132] "Republic of the Congo"               
    ## [133] "Qatar"                               
    ## [134] "Côte d’Ivoire"                       
    ## [135] "Bolivia"                             
    ## [136] "Rwanda"                              
    ## [137] "French Guiana"                       
    ## [138] "Palestinian Territories"             
    ## [139] "The Bahamas"                         
    ## [140] "Sierra Leone"                        
    ## [141] NA                                    
    ## [142] "Bahrain"                             
    ## [143] "Kuwait"                              
    ## [144] "Oman"                                
    ## [145] "Solomon Islands"                     
    ## [146] "Macau"                               
    ## [147] "Burma"                               
    ## [148] "Haiti"                               
    ## [149] "American Samoa"                      
    ## [150] "Algeria"                             
    ## [151] "Niger"                               
    ## [152] "Afghanistan"                         
    ## [153] "Djibouti"                            
    ## [154] "Cuba"                                
    ## [155] "Gabon"                               
    ## [156] "Nicaragua"                           
    ## [157] "Malta"                               
    ## [158] "South Sudan"                         
    ## [159] "Eswatini"                            
    ## [160] "Libya"                               
    ## [161] "Madagascar"                          
    ## [162] "Greenland"                           
    ## [163] "Central African Republic"            
    ## [164] "Guinea-Bissau"                       
    ## [165] "Guam"                                
    ## [166] "North Korea"                         
    ## [167] "Holy See"                            
    ## [168] "Seychelles"                          
    ## [169] "Belize"                              
    ## [170] "Maldives"                            
    ## [171] "Angola"                              
    ## [172] "Azerbaijan"                          
    ## [173] "Togo"                                
    ## [174] "Montenegro"                          
    ## [175] "Bhutan"                              
    ## [176] "Namibia"                             
    ## [177] "Liechtenstein"                       
    ## [178] "Antigua and Barbuda"                 
    ## [179] "Lesotho"                             
    ## [180] "Andorra"                             
    ## [181] "Nauru"                               
    ## [182] "French Southern and Antarctic Lands" 
    ## [183] "Anguilla"                            
    ## [184] "Chad"                                
    ## [185] "Virgin Islands"                      
    ## [186] "Brunei"                              
    ## [187] "Somalia"                             
    ## [188] "Saint Martin"                        
    ## [189] "Cabo Verde"                          
    ## [190] "Samoa"                               
    ## [191] "French Polynesia"                    
    ## [192] "Liberia"                             
    ## [193] "Grenada"                             
    ## [194] "Burundi"                             
    ## [195] "Barbados"                            
    ## [196] "Suriname"                            
    ## [197] "Mauritania"                          
    ## [198] "Cayman Islands"                      
    ## [199] "New Caledonia"                       
    ## [200] "Comoros"                             
    ## [201] "Federal Republic of Yugoslavia"      
    ## [202] "Jersey"                              
    ## [203] "Mayotte"                             
    ## [204] "Saint Kitts and Nevis"               
    ## [205] "Gibraltar"                           
    ## [206] "Northern Mariana Islands"            
    ## [207] "United States Minor Outlying Islands"
    ## [208] "Guyana"                              
    ## [209] "Faroe Islands"                       
    ## [210] "Bermuda"                             
    ## [211] "Aland Islands"                       
    ## [212] "Timor-Leste"                         
    ## [213] "Antarctica"                          
    ## [214] "Equatorial Guinea"                   
    ## [215] "Bonaire, Saint Eustatius and Saba "  
    ## [216] "Tajikistan"                          
    ## [217] "Niue"                                
    ## [218] "Saint Lucia"                         
    ## [219] "Vanuatu"                             
    ## [220] "Montserrat"                          
    ## [221] "Dominica"                            
    ## [222] "Saint Vincent and the Grenadines"    
    ## [223] "Christmas Island"                    
    ## [224] "Eritrea"                             
    ## [225] "Palau"                               
    ## [226] "Micronesia"                          
    ## [227] "Kiribati"                            
    ## [228] "Aruba"                               
    ## [229] "Turkmenistan"                        
    ## [230] "Curacao"

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
    ##  [4] "Turkey (Türkiye)"                    
    ##  [5] "United States"                       
    ##  [6] "Tanzania"                            
    ##  [7] "Serbia and Montenegro"               
    ##  [8] "Hong Kong"                           
    ##  [9] "Iran"                                
    ## [10] "Moldova"                             
    ## [11] "Vietnam"                             
    ## [12] "Monaco"                              
    ## [13] "San Marino"                          
    ## [14] "Laos"                                
    ## [15] "Venezuela"                           
    ## [16] "Martinique"                          
    ## [17] "Netherlands Antilles"                
    ## [18] "South Georgia and the Islands"       
    ## [19] "Kosovo"                              
    ## [20] "Syria"                               
    ## [21] "The Gambia"                          
    ## [22] "Guadeloupe"                          
    ## [23] "Reunion"                             
    ## [24] "Republic of the Congo"               
    ## [25] "Côte d’Ivoire"                       
    ## [26] "Bolivia"                             
    ## [27] "French Guiana"                       
    ## [28] "Palestinian Territories"             
    ## [29] "The Bahamas"                         
    ## [30] NA                                    
    ## [31] "Macau"                               
    ## [32] "Burma"                               
    ## [33] "American Samoa"                      
    ## [34] "Greenland"                           
    ## [35] "Guam"                                
    ## [36] "North Korea"                         
    ## [37] "Holy See"                            
    ## [38] "Liechtenstein"                       
    ## [39] "Andorra"                             
    ## [40] "Nauru"                               
    ## [41] "French Southern and Antarctic Lands" 
    ## [42] "Anguilla"                            
    ## [43] "Virgin Islands"                      
    ## [44] "Brunei"                              
    ## [45] "Saint Martin"                        
    ## [46] "French Polynesia"                    
    ## [47] "Cayman Islands"                      
    ## [48] "New Caledonia"                       
    ## [49] "Federal Republic of Yugoslavia"      
    ## [50] "Jersey"                              
    ## [51] "Mayotte"                             
    ## [52] "Saint Kitts and Nevis"               
    ## [53] "Gibraltar"                           
    ## [54] "Northern Mariana Islands"            
    ## [55] "United States Minor Outlying Islands"
    ## [56] "Faroe Islands"                       
    ## [57] "Bermuda"                             
    ## [58] "Aland Islands"                       
    ## [59] "Antarctica"                          
    ## [60] "Bonaire, Saint Eustatius and Saba "  
    ## [61] "Niue"                                
    ## [62] "Montserrat"                          
    ## [63] "Dominica"                            
    ## [64] "Christmas Island"                    
    ## [65] "Palau"                               
    ## [66] "Micronesia"                          
    ## [67] "Aruba"                               
    ## [68] "Curacao"

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
