# flux-figures.R

annot_pairwise <- function(fluxes) {
  fluxes %>%
    dplyr::filter(experiment %in% c("02", "05", "bay")) %>%
    dplyr::group_by(experiment, cell_type, metabolite, abbreviation) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      ttest = purrr::map_dbl(data, ~t.test(
        flux ~ interaction(oxygen, treatment),
        data = .x,
        paired = TRUE
      )$p.value),
      pval = annot_p(ttest),
      y_pos = Inf,
      vjust = 1.5
    )
}

plot_growth_curve <- function(
  df,
  cell = c("lf", "df", "pasmc"),
  exper = c("02", "05", "bay")
) {
  df %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exper &
        metabolite == "cells" &
        time < 96
    ) %>%
    dplyr::group_by(
      date,
      group,
      time
    ) %>%
    dplyr::summarize(count = mean(conc, na.rm = TRUE)) %>%
    plot_time_lines(
      y = count / 1000,
      ylab = expression(paste("Cell count (x", 10^3, ")")),
      clr = "group"
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, NA),
      clip = "off"
    )
}

annot_fluxes_main <- function(x, formula) {
  x <-
    lmerTest::lmer(
      formula,
      data = x
    ) %>%
    emmeans::emmeans(~ group) %>%
    pairs() %>%
    broom::tidy() %>%
    dplyr::select(contrast, tidyselect::contains("p.value")) %>%
    dplyr::rename_with(~ "pval", .cols = tidyselect::contains("p.value")) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        contrast == "21% - 0.5%" ~ "0.5%",
        contrast == "21% - 0.2%" ~ "0.2%",
        contrast == "DMSO - BAY" ~ "BAY",
        contrast == "0.5% - BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY")),
      label = dplyr::case_when(
        pval < 0.05 & contrast == "21% - 0.5%" ~ "*",
        pval < 0.05 & contrast == "21% - 0.2%" ~ "*",
        pval < 0.05 & contrast == "DMSO - BAY" ~ "†",
        pval < 0.05 & contrast == "0.5% - BAY" ~ "‡"
      ),
      vjust = dplyr::case_when(
        label == "*" ~ 1.5,
        label == "†" ~ 1.5,
        label == "‡" ~ 3
      ),
      y = Inf
    ) %>%
    dplyr::filter(!is.na(label))

  y <-
    setdiff(c("21%", "0.5%", "0.2%", "DMSO", "BAY"), x$group) %>%
    factor(levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY")) %>%
    list() %>%
    rlang::set_names("group")

  dplyr::bind_rows(x, y)
}

plot_growth_rates <- function(
  df,
  cell = c("lf", "pasmc"),
  exper = c("02", "05", "bay"),
  annot = NULL
) {
  x <-
    df %>%
    dplyr::filter(cell_type %in% cell & experiment %in% exper)

  annot <-
    annot_fluxes_main(x, mu ~ group + (1 | date)) %>%
    dplyr::filter(group %in% x$group)

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = group,
      y = mu
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      width = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = "mean_se",
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        label = label,
        y = y,
        vjust = vjust
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Treatment",
      y = "Growth rate (/h)"
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    theme_plots() +
    ggplot2::coord_cartesian(
      ylim = c(0, 0.03),
      clip = "off"
    ) +
    NULL
}

plot_high_fluxes <- function(
  df,
  cell = c("lf", "pasmc"),
  exper = c("02", "05", "bay")
) {

  x <-
    df %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exper &
        metabolite %in% c("lactate", "glucose")
    )

  annot <-
    x %>%
    dplyr::group_by(abbreviation) %>%
    tidyr::nest() %>%
    dplyr::mutate(annot = purrr::map(data, annot_fluxes_main, formula = flux ~ group + (1 | date))) %>%
    tidyr::unnest(c(annot)) %>%
    dplyr::filter(group %in% x$group)

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = abbreviation,
        y = y,
        vjust = vjust,
        label = label,
        group = group
      ),
      position = ggplot2::position_dodge(width = 0.6),
      hjust = 0.5,
      family = "Calibri",
      size = 6/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = 0.2),
      breaks = scales::extended_breaks(n = 7)
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_low_fluxes <- function(
  df,
  cell = c("lf", "pasmc"),
  exper = c("02", "05", "bay")
) {

  x <-
    df %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exper &
        metabolite %nin% c("lactate", "glucose")
    )

  annot <-
    x %>%
    dplyr::group_by(abbreviation) %>%
    tidyr::nest() %>%
    dplyr::mutate(annot = purrr::map(data, annot_fluxes_main, formula = flux ~ group + (1 | date))) %>%
    tidyr::unnest(c(annot)) %>%
    dplyr::filter(group %in% x$group)

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = abbreviation,
        y = y,
        vjust = vjust,
        label = label,
        group = group
      ),
      position = ggplot2::position_dodge(width = 0.6),
      family = "Calibri",
      size = 6/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      trans = ggallin::pseudolog10_trans,
      breaks = c(-100, -10, 0, 10, 100),
      limits = c(-250, 250)
    ) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

clean_viability <- function(viability_file) {
  readr::read_csv(viability_file) %>%
    dplyr::mutate(
      viability = 100 * live / (dead + live),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) %>%
    dplyr::group_by(time, oxygen) %>%
    wmo::remove_nested_outliers(viability, remove = TRUE)
}

plot_cells_per_dna <- function(dna_per_cell_clean) {
  dna_per_cell_clean %>%
    dplyr::filter(.data$volume == 200 & cells < 400000) %>%
    ggplot2::ggplot() +
    ggplot2::facet_wrap(~cell_type, labeller = ggplot2::as_labeller(toupper)) +
    ggplot2::aes(
      x = cells/1000,
      y = conc
    ) +
    ggplot2::geom_smooth(
      formula = y ~ 0 + x,
      method = "lm",
      color = clrs[[2]],
      size = 0.5,
      se = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.data = "mean_se",
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      fill = "black",
      size = 2,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = expression(paste("Cell count (x", 10^3, ")")),
      y = "DNA (ng)"
    ) +
    theme_plots() +
    ggplot2::coord_cartesian(xlim = c(0, NA), clip = "off") +
    NULL
}

clean_dna_count_hypoxia <- function(dna_count_hypoxia_file) {
  df <-
    readr::read_csv(dna_count_hypoxia_file) %>%
    dplyr::mutate(oxygen = factor(oxygen, levels = c("21%", "0.5%"))) %>%
    clean_technical_replicates()

  std <-
    df %>%
    dplyr::filter(!is.na(conc)) %>%
    lm(value ~ conc, data = .)

  df %>%
    dplyr::filter(is.na(conc)) %>%
    dplyr::mutate(conc = interpolate(., std)) %>%
    dplyr::select(-value)
}

plot_dna_count_hypoxia <- function(dna_count_hypoxia) {
  ggplot2::ggplot(dna_count_hypoxia) +
    ggplot2::aes(
      x = count / 1000,
      y = conc,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      size = 0.5,
      formula = y ~ x,
      se = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      pch = 21,
      color = "white",
      # alpha = 0.3,
      size = 2,
      show.legend = FALSE) +
    ggplot2::labs(
      x = expression(paste("Cell count (×", 10^3, ")")),
      y = "DNA (ng)"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 300, 100)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    theme_plots() +
    ggplot2::coord_cartesian(
      ylim = c(0, NA),
      xlim = c(0, 300),
      clip = "off"
    )
}

plot_k <- function(degradation_rates, k) {
  annot <-
    k %>%
    dplyr::select(-k) %>%
    dplyr::mutate(
      label = "*",
      ypos = Inf,
      vjust = 1.5
    )

  degradation_rates %>%
    dplyr::mutate(
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" ~ "0.5%",
        treatment == "DMSO" ~ "DMSO"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO"))
    ) %>%
    dplyr::left_join(annot, by = c("metabolite", "oxygen", "treatment")) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), k),
      y = k
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        color = group,
        label = label,
        y = ypos,
        vjust = vjust
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Rate constant (/h)",
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      limits = force
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_evap_data <- function(evap_clean) {
  evap_clean %>%
    dplyr::filter(experiment == "05" & cell_type == "lf") %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = volume,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
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
      y = "Volume (mL)"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, 24)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    theme_plots()
}
