read_daly = function(year) {
  path = list.files("data", pattern = paste0("^daly.*_", year, "[_.].*xlsx$"), full.names = TRUE)
  
  raw = read_xlsx(path, sheet = "All ages", skip = 6, .name_repair = "unique_quiet")
  
  iso_lookup = raw |>
    slice(1) |>
    select(-(1:7)) |>
    mutate(across(everything(), as.character)) |>
    pivot_longer(everything(), names_to = "country", values_to = "iso3")
  
  raw |>
    rename(sex = 1, ghe_code = 2, col3 = 3, col4 = 4, col5 = 5, col6 = 6, col7 = 7) |>
    filter(sex == "Persons", !is.na(ghe_code)) |>
    mutate(
      level = case_when(
        ghe_code == 0 ~ 0,
        !is.na(col7)  ~ 4,
        !is.na(col6)  ~ 3,
        !is.na(col5)  ~ 2,
        TRUE          ~ 1
      ),
      cause = case_when(
        level == 4 ~ col7,
        level == 3 ~ col6,
        level == 2 ~ col5,
        TRUE       ~ col4
      ),
      cause = str_trim(cause)
    ) |>
    select(-c(sex, col3:col7)) |>
    pivot_longer(
      -c(ghe_code, level, cause),
      names_to = "country",
      values_to = "daly_thousands"
    ) |>
    mutate(
      daly_thousands = as.numeric(daly_thousands),
      year = year
    ) |>
    left_join(iso_lookup, by = "country") |>
    relocate(iso3, .after = country)
}