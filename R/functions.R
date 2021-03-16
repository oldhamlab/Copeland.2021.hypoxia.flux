# functions.R

read_multi_excel <- function(excel_file) {
  sheets <- readxl::excel_sheets(excel_file)
  purrr::map(sheets, ~readxl::read_excel(excel_file, sheet = .x)) %>%
    rlang::set_names(sheets)
}

clean_technical_replicates <- function(tbl) {
  tidyr::pivot_longer(
    data = tbl,
    cols = .data$a:.data$c,
    names_to = "replicate",
    values_to = "value"
  ) %>%
    dplyr::group_by(
      dplyr::across(-c(.data$replicate, .data$value))
    ) %>%
    dplyr::mutate(
      value = replace(
        .data$value,
        which(abs(.data$value - median(.data$value)) / mad(.data$value) > 2),
        NA
      )
    ) %>%
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE)) %>%
    dplyr::ungroup()
}

data_file <-
    system.file(
      "extdata/dna-per-cell-number.xlsx",
      package = "Copeland.2021.hypoxia.flux"
    )

read_multi_excel(data_file) %>%
  purrr::map(clean_technical_replicates) %>%
  bind_rows(.id = "id") %>%
  separate(id, c("cell_type", "volume"), sep = "_", convert = TRUE) %>%
  mutate(cells = 1000 * cells)
