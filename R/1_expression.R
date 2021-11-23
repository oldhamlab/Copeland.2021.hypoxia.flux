# expression.R

read_data <- function(data_files) {
  data_files[stringr::str_detect(data_files, "\\.csv$")] %>%
    rlang::set_names(stringr::str_extract(., "(lf|pasmc)_(02|05-bay|05|bay|)")) %>%
    purrr::map_dfr(read_csv, .id = "experiment") %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%", "0.2%"), ordered = TRUE),
      treatment = factor(treatment, levels = c("None", "DMSO", "BAY"), ordered = TRUE)
    )
}

plot_blot <- function(blot_image, scale = 1.3, hjust = 0.15, vjust = 0.1) {
  # blot_image <- magick::image_read_pdf(blot_image)
  blot_image <- magick::image_read(blot_image)

  cowplot::ggdraw() +
    cowplot::draw_image(
      blot_image,
      scale = scale,
      hjust = hjust,
      vjust = vjust
    ) +
    wmo::theme_wmo(
      base_family = "Calibri",
      base_size = 8
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      panel.border = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold")
    )
}

normalize_densities <- function(blot_raw)
  {
  blot_raw %>%
    dplyr::filter(.data$time < 96) %>%
    dplyr::filter(!(experiment == "lf_05-bay" & gel %in% c("b", "e", "f"))) %>%
    tidyr::pivot_longer(
      blot:tidyselect::last_col(),
      names_to = "protein",
      values_to = "value",
      values_drop_na = TRUE
    ) %>%
    dplyr::group_by(.data$experiment, .data$gel, .data$protein) %>%
    dplyr::mutate(norm = value / mean(value, na.rm = TRUE)) %>%
    dplyr::group_by(dplyr::across(.data$experiment:.data$time)) %>%
    dplyr::mutate(density = norm / norm[protein == "blot"]) %>%
    dplyr::filter(protein != "blot") %>%
    dplyr::group_by(.data$experiment, .data$protein) %>%
    dplyr::mutate(
      fold_change = density /
        mean(density[oxygen == min(oxygen) & treatment %in% c("None", "DMSO") & time == min(time)]),
      fold_change = replace(fold_change, experiment == "lf_bay" & protein == "hif1a", sqrt(fold_change)),
      group = dplyr::case_when(
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "None" & oxygen == "21%" ~ "21%",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "DMSO" ~ "DMSO",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "BAY" ~ "BAY",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & oxygen == "0.5%" ~ "0.5%",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) %>%
    dplyr::group_by(experiment, oxygen, treatment, group, time, protein) %>%
    wmo::remove_nested_outliers(fold_change, remove = TRUE) %>%
    dplyr::relocate(group, .after = treatment)
}

normalize_qpcr <- function(raw_mrna) {
  raw_mrna %>%
    dplyr::mutate(gene = tolower(gene)) %>%
    dplyr::group_by(dplyr::across(c(experiment:gene))) %>%
    dplyr::summarize(ct = mean(ct, na.rm = TRUE)) %>%
    tidyr::pivot_wider(names_from = gene, values_from = ct) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dplyr::across(-c(experiment:time, actin), ~ . - actin)
    ) %>%
    dplyr::select(-actin) %>%
    tidyr::pivot_longer(
      -c(experiment:time),
      names_to = "protein",
      values_to = "dct"
    ) %>%
    dplyr::filter(!is.na(dct)) %>%
    dplyr::group_by(.data$experiment, .data$protein) %>%
    dplyr::mutate(
      ddct = dct - mean(
        dct[oxygen == min(oxygen) & treatment %in% c("None", "DMSO") & time == 0]
      ),
      fold_change = 2 ^ -ddct,
      group = dplyr::case_when(
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "None" & oxygen == "21%" ~ "21%",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "DMSO" ~ "DMSO",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "BAY" ~ "BAY",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & oxygen == "0.5%" ~ "0.5%",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) %>%
    dplyr::group_by(experiment, oxygen, treatment, group, time, protein) %>%
    wmo::remove_nested_outliers(fold_change, remove = TRUE) %>%
    dplyr::relocate(group, .after = treatment)
}

plot_expression <- function(
  df,
  exp,
  prot = c("ldha", "hif1a"),
  ylab = prot
) {

  df %>%
    dplyr::filter(experiment %in% exp) %>%
    dplyr::filter(protein == prot) %>%
    plot_time_lines(y = fold_change, ylab = ylab, clr = "group")
}

plot_time_lines <- function(
  df,
  y,
  ylab = prot,
  clr = c("oxygen", "treatment", "group")
) {

  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = time,
      y = {{y}},
      color = .data[[clr]],
      fill = .data[[clr]]
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.data = "mean_se",
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "line",
      fun.data = "mean_se",
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      size = 2,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Time (h)",
      y = ylab
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, 24)) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots()
}
