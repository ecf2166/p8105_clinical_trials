read_aact = function(table, ...) {
  read_delim(
    file.path("raw_data", paste0(table, ".txt")),
    delim = "|",
    col_types = cols(.default = "c"),
    ...
  )
}