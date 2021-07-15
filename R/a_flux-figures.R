# flux-figures.R

plot_growth_curve <- function(
  df,
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay")
) {
  df %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
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
  lmerTest::lmer(
    formula,
    data = x
  ) %>%
    emmeans::emmeans(~ group) %>%
    pairs() %>%
    broom::tidy() %>%
    dplyr::select(contrast, pval = adj.p.value) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        contrast == "21% - 0.5%" ~ "0.5%",
        contrast == "DMSO - BAY" ~ "BAY",
        contrast == "0.5% - BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY")),
      label = dplyr::case_when(
        pval < 0.05 & contrast == "21% - 0.5%" ~ "*",
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
    dplyr::filter(!is.na(label)) %>%
    dplyr::bind_rows(list(group = factor(
      c("21%", "0.5%", "DMSO", "BAY"),
      levels = c("21%", "0.5%", "DMSO", "BAY")))
    )
}

plot_growth_rates <- function(
  df,
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  annot = NULL
) {
  x <-
    df %>%
    dplyr::filter(cell_type %in% cell & experiment %in% exp)

  annot <-annot_fluxes_main(x, mu ~ group + (1 | date))

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
  exp = c("02", "05", "bay")
) {

  x <-
    df %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
        metabolite %in% c("lactate", "glucose")
    )

  annot <-
    x %>%
    dplyr::group_by(abbreviation) %>%
    tidyr::nest() %>%
    dplyr::mutate(annot = purrr::map(data, annot_fluxes_main, formula = flux ~ group + (1 | date))) %>%
    tidyr::unnest(c(annot))

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
      show.legend = FALSE
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
      y = "Flux (fmol/cell/h)"
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    theme_plots()
}

plot_low_fluxes <- function(
  df,
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay")
) {

  x <-
    df %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
        metabolite %nin% c("lactate", "glucose")
    )

  annot <-
    x %>%
    dplyr::group_by(abbreviation) %>%
    tidyr::nest() %>%
    dplyr::mutate(annot = purrr::map(data, annot_fluxes_main, formula = flux ~ group + (1 | date))) %>%
    tidyr::unnest(c(annot))

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
