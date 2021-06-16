# mid.R

# calculate_ratios --------------------------------------------------------

calculate_ratios <- function(file_name) {
  readr::read_csv(file_name) %>%
    dplyr::filter(!is.na(area)) %>%
    tidyr::separate(
      .data$filename,
      into = c(NA, "window", "replicate"),
      sep = c(1, 2),
      convert = TRUE
    ) %>%
    tidyr::separate(
      .data$ion,
      into = c("metabolite", "isotope"),
      sep = " ",
      convert = TRUE
    ) %>%
    dplyr::mutate(carbons = case_when(
      .data$metabolite %in% c("citrate") ~ 6,
      .data$metabolite %in% c("2HG", "2OG", "glutamate", "glutamine") ~ 5,
      .data$metabolite %in% c("aspartate", "malate", "succinate") ~ 4,
      .data$metabolite %in% c("lactate", "pyruvate", "alanine", "serine") ~ 3
    )) %>%
    dplyr::filter(.data$window <= .data$carbons) %>%
    tidyr::pivot_wider(names_from = .data$isotope, values_from = .data$area) %>%
    dplyr::mutate(ratio = .data$M1 / .data$M0) %>%
    dplyr::group_by(.data$metabolite)
}


# import_qbias ------------------------------------------------------------

import_qbias <- function(file_list) {
  file_list[stringr::str_detect(file_list, "\\.csv$")] %>%
    rlang::set_names(
      stringr::str_extract(basename(.), pattern = "(?<=_)\\w(?=\\.csv)")
    ) %>%
    purrr::map_dfr(calculate_ratios, .id = "batch") %>%
    dplyr::group_by(batch, metabolite) %>%
    dplyr::arrange(metabolite)
}


# calculate_predicted_ratios ----------------------------------------------

calculate_predicted_ratios <- function() {
  predicted_ratios <-
    isotope_library %>%
    dplyr::mutate(table = purrr::map2(
      formula,
      polarity,
      ~mzrtools::mz_iso_quant(molecule = .x, polarity = .y)[["prob_matrix"]])
    ) %>%
    dplyr::mutate(pred_ratio = purrr::map_dbl(table, ~ .x[[2, 1]]/.x[[1, 1]])) %>%
    dplyr::select(.data$metabolite, .data$pred_ratio)
}

# calculate_correction_factors --------------------------------------------

calculate_correction_factors <- function(qbias_ratios, pred_ratios) {
  qbias_ratios %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(.data$data, ~MASS::rlm(ratio ~ poly(window, 3), data = .x)),
      predict = purrr::map2(.data$model, .data$data, predict)
    ) %>%
    tidyr::unnest(c(data, predict)) %>%
    dplyr::select(.data$batch, .data$metabolite, .data$window, .data$carbons, .data$predict) %>%
    dplyr::distinct() %>%
    dplyr::left_join(pred_ratios, by = "metabolite") %>%
    dplyr::mutate(cf = .data$predict / .data$pred_ratio) %>%
    dplyr::select(.data$metabolite, "M" = .data$window, .data$carbons, .data$cf) %>%
    dplyr::filter(.data$M < .data$carbons) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(M = M + 1) %>%
    tidyr::pivot_wider(names_from = .data$M, values_from = .data$cf) %>%
    dplyr:: mutate(
      M0 = 1,
      M1 = 1 / .data$`1` * .data$M0,
      M2 = 1 / .data$`2` * .data$M1,
      M3 = 1 / .data$`3` * .data$M2,
      M4 = 1 / .data$`4` * .data$M3,
      M5 = 1 / .data$`5` * .data$M4,
      M6 = 1 / .data$`6` * .data$M5
    ) %>%
    dplyr::select(.data$batch, .data$metabolite, tidyselect::matches("M[0-9]+")) %>%
    tidyr::pivot_longer(
      cols = matches("M[0-9]+"),
      names_to = "M",
      values_to = "cf",
      values_drop_na = TRUE
    ) %>%
    dplyr::arrange(.data$batch, .data$metabolite)
}

# clean_mids --------------------------------------------------------------

clean_mids <- function(mid_files) {
  read_function <- function(filename) {
    readr::read_csv(filename) %>%
      replace(is.na(.), 0)
  }

  mid_files %>%
    rlang::set_names(sub("\\..*$", "", basename(.))) %>%
    purrr::map_dfr(read_function, .id = "experiment") %>%
    tidyr::pivot_longer(
      cols = .data$`2HG M0`:tidyselect::last_col(),
      names_to = "ion",
      values_to = "peak_area",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(
      batch = stringr::str_extract(.data$experiment, "a|b|c"),
      method = stringr::str_extract(.data$experiment, "fs|sim"),
      cell_type = stringr::str_extract(.data$experiment, "lf|pasmc"),
      date = stringr::str_extract(.data$experiment, "\\d{4}-\\d{2}-\\d{2}"),
      .before = .data$experiment
    ) %>%
    dplyr::select(-.data$experiment) %>%
    tidyr::separate(
      .data$id,
      c("tracer", "treatment", "time", "well"),
      sep = " "
    ) %>%
    tidyr::separate(.data$ion, into = c("metabolite", "isotope")) %>%
    dplyr::mutate(
      peak_area = tidyr::replace_na(peak_area, 0),
      oxygen = dplyr::case_when(
        .data$treatment %in% c("21%", "dmso", "bay") ~ "21%",
        .data$treatment == "0.5%" ~ "0.5%",
        .data$treatment == "0.2%" ~ "0.2%"
      ),
      treatment = dplyr::case_when(
        .data$treatment %in% c("21%", "0.5%", "0.2%") ~ "None",
        TRUE ~ toupper(.data$treatment)
      ),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%", "0.2%")),
      treatment = factor(.data$treatment, levels = c("None", "DMSO", "BAY")),
      time = as.double(.data$time),
      time = dplyr::case_when(
        .data$cell_type == "lf" ~ time * 24,
        .data$cell_type == "pasmc" ~ time * 12
      )
    ) %>%
    dplyr::relocate(.data$oxygen, .before = .data$treatment) %>%
    dplyr::select(-filename)
}

# correct_mid -------------------------------------------------------------

correct_mid <- function(mid_clean){
  correction_factors <-
    dplyr::mutate(correction_factors, method = "sim")

  mid_q <-
    mid_clean %>%
    dplyr::left_join(
      correction_factors,
      by = c("batch", "method", "metabolite", "isotope" = "M")
    ) %>%
    dplyr::mutate(
      area_corr = dplyr::case_when(
        .data$method == "sim" ~ peak_area * cf,
        .data$method == "fs" ~ peak_area * 1
      )
    ) %>%
    dplyr::select(-c(.data$peak_area, .data$cf))

  correction_matrices <-
    isotope_library %>%
    dplyr::mutate(
      matrix = purrr::map2(formula, polarity, mzrtools::mz_iso_quant),
      matrix = purrr::map(matrix, purrr::pluck, "prob_matrix")
    )

  mmult <- function(m, df) {
    mid <- df$mid
    if (nrow(m) > length(mid)) {
      m <- m[1:length(mid), 1:length(mid)]
    }
    mid_corr <- mzrtools::mz_iso_correct(m, mid)
    dplyr::bind_cols(df, mid_corr = mid_corr)
  }

  mid_na <-
    mid_q %>%
    dplyr::group_by(
      .data$method,
      .data$cell_type,
      .data$date,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$time,
      .data$well,
      .data$metabolite
    ) %>%
    dplyr::filter(.data$metabolite %nin% c("palmitate", "sedoheptulose")) %>%
    dplyr::mutate(mid = .data$area_corr / sum(.data$area_corr)) %>%
    dplyr::filter(.data$mid != "NaN") %>%
    dplyr::select(-.data$area_corr) %>%
    tidyr::nest() %>%
    dplyr::left_join(
      dplyr::select(correction_matrices, -c(.data$formula, .data$polarity)),
      by = "metabolite"
    ) %>%
    dplyr::mutate(data = purrr::map2(.data$matrix, .data$data, mmult)) %>%
    tidyr::unnest(c(.data$data)) %>%
    dplyr::filter(.data$isotope %in% stringr::str_c("M", 0:6)) %>%
    dplyr::select(
      .data$method:.data$metabolite,
      .data$isotope,
      .data$tracer:.data$well,
      .data$mid,
      .data$mid_corr
    )

}

# remove_mid_outliers -----------------------------------------------------

remove_mid_outliers <- function(mid_correct) {
  mid_correct %>%
    dplyr::group_by(
      dplyr::across(c(.data$method:.data$time, .data$metabolite:.data$isotope))
    ) %>%
    dplyr::mutate(
      outlier = replace_outliers(.data$mid_corr)
    ) %>%
    dplyr::summarise(mid = mean(mid_corr, na.rm = TRUE)) %>%
    dplyr::arrange(
      .data$cell_type,
      .data$metabolite,
      .data$oxygen,
      .data$treatment,
      .data$tracer
    )
}

# replace_outliers --------------------------------------------------------

replace_outliers <- function(vec) {
  if (mad(vec, na.rm = TRUE) == 0) return (vec)
  replace(
    vec,
    abs(vec - median(vec, na.rm = TRUE)) / mad(vec, na.rm = TRUE) > 2,
    NA
  )
}

# plot_mid_curves ---------------------------------------------------------

plot_mid_curves <- function(mids) {
  mids %>%
    dplyr::group_by(
      .data$method,
      .data$cell_type,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$metabolite
    ) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      title = stringr::str_c(
        .data$metabolite,
        .data$method,
        .data$tracer,
        .data$cell_type,
        stringr::str_replace(.data$oxygen, "%", ""),
        .data$treatment,
        sep = "_"
      ),
      plots = purrr::map2(
        .data$data,
        .data$title,
        ~ggplot2::ggplot(.x) +
          ggplot2::aes(
            x = .data$time,
            y = .data$mid,
            color = .data$isotope
          ) +
          ggplot2::geom_point(
            alpha = 0.3
          ) +
          ggplot2::stat_summary(
            fun = "mean",
            geom = "point",
            alpha = 0.8
          ) +
          ggplot2::stat_summary(
            fun = "mean",
            geom = "line",
            alpha = 0.8,
            show.legend = FALSE
          ) +
          ggplot2::scale_x_continuous(
            name = "Time (h)",
            breaks = seq(24, 96, by = 24)
          ) +
          ggplot2::scale_color_manual(values = viridis::viridis(7)) +
          ggplot2::labs(
            y = "Mole fraction",
            title = .y
          )
      )
    )
}
