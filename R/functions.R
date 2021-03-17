# functions.R

path_to_data <- function(nm) {
  system.file(
    stringr::str_c("extdata/", nm),
    package = "Copeland.2021.hypoxia.flux"
  )
}


path_to_reports <- function(nm) {
  system.file(
    stringr::str_c("analysis/", nm),
    package = "Copeland.2021.hypoxia.flux"
  )
}


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


make_std_curves <- function(df) {
  df %>%
    dplyr::filter(!is.na(.data$conc)) %>%
    dplyr::select(where(~all(!is.na(.)))) %>%
    dplyr::group_by(dplyr::across(-c(.data$conc, .data$value))) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      title = stringr::str_c(!!!rlang::syms(dplyr::groups(.)), sep = "_")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      model = furrr::future_map(
        .data$data,
        ~lm(value ~ conc, data = .x, na.action = modelr::na.warn)
      ),
      summary = furrr::future_map(.data$model, broom::glance),
      plots = furrr::future_map2(.data$data, .data$title, make_std_plots)
    ) %>%
    dplyr::group_by(
      dplyr::across(
        -c(.data$data, .data$title, .data$model, .data$summary, .data$plots)
      )
    )
}


make_std_plots <- function(df, title = NULL) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = .data$conc,
      y = .data$value
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      color = "gray20",
      se = FALSE
    ) +
    ggplot2::geom_point(
      size = 3,
      alpha = 0.3,
      color = "blue"
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8,
      color = "blue"
    ) +
    ggplot2::labs(
      x = "Concentration",
      y = "Value",
      title = title
    )
}


interp_data <- function(tbl, std) {
  tbl %>%
    dplyr::filter(is.na(.data$conc)) %>%
    dplyr::select(where(~all(!is.na(.)))) %>%
    dplyr::group_by(dplyr::across(dplyr::group_vars(std))) %>%
    tidyr::nest() %>%
    dplyr::left_join(dplyr::select(std, .data$model)) %>%
    dplyr::mutate(conc = purrr::map2(.data$data, .data$model, wmo::interpolate)) %>%
    tidyr::unnest(c(.data$data, .data$conc)) %>%
    dplyr::select(-c(.data$model, .data$value))
}


interpolate <- function(new_df, model) {
  x <- stats::model.frame(model)[[deparse(model$terms[[3]])]]
  p <- polynom::polynomial(stats::coefficients(model))
  new_y <- as.list(new_df[[deparse(model$terms[[2]])]])
  new_x <- unlist(lapply(new_y, function(y) {
    roots <- solve(p, y)
    roots <- round(roots, digits = 8)
    root <- roots[which(Im(roots) == 0 & Re(roots) >= 0 & Re(roots) <= 1.25 * max(x))]
    ifelse(identical(root, numeric(0)), NA, Re(root))
  }))
  new_x
}


clean_dna_per_cell <- function(filename) {
  filename %>%
    read_multi_excel() %>%
    purrr::map(clean_technical_replicates) %>%
    dplyr::bind_rows(.id = "id") %>%
    tidyr::separate(id, c("cell_type", "volume"), sep = "_", convert = TRUE) %>%
    dplyr::mutate(cells = 1000 * cells)
}


calculate_cells_per_dna <- function(tbl) {
  tbl %>%
    dplyr::filter(!(cell_type == "lf" & volume == "100" & cells == 300000)) %>%
    dplyr::filter(!(cell_type == "pasmc" & volume == "200" & cells == 400000)) %>%
    dplyr::group_by(cell_type, volume) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = map(data, ~lm(conc ~ 0 + cells, data = .x, na.action = modelr::na.warn)),
      glance = map(model, broom::tidy)
    ) %>%
    tidyr::unnest(c(glance)) %>%
    dplyr::select(cell_type, volume, slope = estimate) %>%
    dplyr::mutate(slope = 1/slope)
}


