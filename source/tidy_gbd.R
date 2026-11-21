# Read and tidy IHME Global Burden of Disease DALY exports

library(tidyverse)
library(countrycode)

# GBD level-2 injury causes are combined into one analysis category, because
# MeSH condition terms describe the injury broadly, not its specific cause
gbd_injury_causes = c("Transport injuries", "Unintentional injuries",
                      "Self-harm and interpersonal violence")

gbd_aggregate_locations = c("World Bank Income Levels", "Global")

read_gbd = function(path) {
  read_csv(path, show_col_types = FALSE) |>
    rename(any_of(c(measure  = "measure_name",  location = "location_name",
                    sex      = "sex_name",      age      = "age_name",
                    cause    = "cause_name",    metric   = "metric_name"))) |>
    filter(str_detect(measure, "^DALYs"), sex == "Both", age == "All ages",
           metric == "Number") |>
    transmute(location, year = as.integer(year), gbd_cause = cause,
              daly = val, daly_lower = lower, daly_upper = upper)
}

# One row per country x year x GBD level-2 cause (22 causes)
tidy_gbd_countries = function(paths) {
  map(paths, read_gbd) |>
    list_rbind() |>
    filter(!location %in% gbd_aggregate_locations,
           !str_detect(location, "^World Bank")) |>
    distinct() |>
    mutate(
      iso3 = countrycode(location, origin = "country.name", destination = "iso3c"),
      analysis_cause = if_else(gbd_cause %in% gbd_injury_causes, "Injuries", gbd_cause)
    ) |>
    relocate(iso3, location, year)
}

# One row per country x year x analysis cause (20 categories), with each cause's
# share of the country's DALYs that year.
summarize_gbd_analysis = function(gbd_countries) {
  gbd_countries |>
    group_by(iso3, location, year, analysis_cause) |>
    summarize(daly = sum(daly), .groups = "drop") |>
    group_by(iso3, year) |>
    mutate(daly_share = daly / sum(daly)) |>
    ungroup()
}

# One row per World Bank income group x year x analysis cause.
tidy_gbd_income = function(path) {
  read_gbd(path) |>
    mutate(
      income_group = str_remove(location, "^World Bank ") |>
        factor(levels = c("Low Income", "Lower Middle Income",
                          "Upper Middle Income", "High Income")),
      analysis_cause = if_else(gbd_cause %in% gbd_injury_causes, "Injuries", gbd_cause)
    ) |>
    group_by(income_group, year, analysis_cause) |>
    summarize(daly = sum(daly), .groups = "drop") |>
    group_by(income_group, year) |>
    mutate(daly_share = daly / sum(daly)) |>
    ungroup()
}

# Country -> current World Bank income group
# Note: this is the *current* classification and some countries have moved groups since 2008.
income_group_lookup = function() {
  WDI::WDI_data$country |>
    as_tibble() |>
    filter(income != "Aggregates", income != "Not classified") |>
    transmute(iso3 = iso3c,
              income_group = factor(
                str_to_title(income),
                levels = c("Low Income", "Lower Middle Income",
                           "Upper Middle Income", "High Income")))
}
