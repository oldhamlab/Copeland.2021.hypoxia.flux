# mid-figures.R

plot_labeling_rate <- function(mids) {
  df <-
    mids %>%
    dplyr::filter(
      .data$method == "sim" &
        .data$cell_type == "lf" &
        .data$tracer == "glc6" &
        .data$metabolite %in% c("pyruvate") &
        .data$isotope == "M0"
    ) %>%
    dplyr::mutate(
      labeled = 1 - mid,
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" ~ "0.5%",
        treatment == "DMSO" ~ "DMSO",
        treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY")),
      idx = dplyr::case_when(
        group %in% c("21%", "DMSO") ~ 1,
        group %in% c("0.5%", "BAY") ~ 2
      )
    )

  a <-
    df %>%
    dplyr::group_by(.data$time, .data$group) %>%
    wmo::remove_nested_outliers(labeled, remove = TRUE) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = labeled,
      color = group,
      fill = group
    ) +
    ggplot2::geom_smooth(
      method = "nls",
      formula = y ~ SSasympOrig(x, asym, lrc),
      fullrange = TRUE,
      se = FALSE,
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
    ggplot2::scale_x_continuous(
      breaks = seq(0, 72, 24),
      limits = c(0, 72)
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::labs(
      x = "Time (h)",
      y = "Pyruvate labeled",
      color = NULL,
      fill = NULL
    ) +
    theme_plots()

  k <-
    df %>%
    dplyr::group_by(.data$group) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      m = purrr::map(data, ~nls(labeled ~ SSasympOrig(time, asym, lrc), data = .x)),
      s = purrr::map(m, broom::tidy)
    ) %>%
    tidyr::unnest(c(s)) %>%
    dplyr::filter(term == "lrc") %>%
    dplyr::mutate(
      k = exp(estimate),
      experiment = dplyr::case_when(
        group %in% c("21%", "0.5%") ~ "hypoxia",
        group %in% c("DMSO", "BAY") ~ "bay"
      ),
      sem = std.error / abs(estimate) * k
    )

  b <-
    k %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = group,
      y = k,
      fill = group
    ) +
    ggplot2::geom_col(
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = k - sem,
        ymax = k + sem
      ),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::labs(
      x = NULL,
      y = "Rate",
      color = NULL,
      fill = NULL
    ) +
    theme_plots()

  list(curve = a, rate = b)
}

plot_manuscript_mids <- function(df) {
  df %>%
    dplyr::filter(cell_type == "lf") %>%
    dplyr::mutate(
      metabolite = factor(
        metabolite,
        levels = c("FBP", "3PG", "PYR", "CIT", "AKG", "MAL"),
        labels = c("FBP", "`3PG`", "PYR", "CIT", "AKG", "MAL")
      )
    ) %>%
    dplyr::filter(tracer != "lac3" & time == 72 & !is.na(metabolite)) %>%
    plot_mids()
}

annot_mids_main <- function(a, formula) {
  lmerTest::lmer(
    mid ~ group * isotope + (1 | date),
    data = a
  ) %>%
    emmeans::emmeans(~ group * isotope) %>%
    emmeans::mvcontrast("pairwise", mult.name = "isotope") %>%
    tibble::as_tibble() %>%
    dplyr::select(contrast, tidyselect::contains("p.value")) %>%
    dplyr::rename_with(~ "pval", .cols = tidyselect::contains("p.value")) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        pval < 0.05 & contrast == "21% - 0.5%" ~ "*",
        pval < 0.05 & contrast == "21% - 0.2%" ~ "*",
        pval < 0.05 & contrast == "DMSO - BAY" ~ "†",
        pval < 0.05 & contrast == "0.5% - BAY" ~ "‡"
      ),
      order = dplyr::case_when(
        contrast == "21% - 0.5%" ~ 1,
        contrast == "21% - 0.2%" ~ 1,
        contrast == "DMSO - BAY" ~ 2,
        contrast == "0.5% - BAY" ~ 3
      )
    ) %>%
    dplyr::filter(!is.na(label)) %>%
    dplyr::arrange(order) %>%
    dplyr::pull(label) %>%
    stringr::str_c(collapse = " ")
}

plot_mids <- function(df) {

  tracer_labels <-
    c(expression(paste("[1,2-"^13, "C"[2], "] glucose")),
      expression(paste("[U-"^13, "C"[6], "] glucose")),
      expression(paste("[U-"^13, "C"[5], "] glutamine")),
      expression(paste("[U-"^13, "C"[3], "] lactate")))

  tracer_levels <-
    c("glc2", "glc6", "q5", "lac3")

  x <-
    df %>%
    dplyr::mutate(
      tracer = factor(tracer, levels = tracer_levels, labels = tracer_labels),
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" & treatment == "None" ~ "0.5%",
        oxygen == "21%" & treatment == "DMSO" ~ "DMSO",
        oxygen == "21%" & treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY"))
    )

  annot <-
    x %>%
    dplyr::group_by(tracer, metabolite) %>%
    tidyr::nest() %>%
    dplyr::mutate(annot = purrr::map(data, annot_mids_main)) %>%
    tidyr::unnest(c(annot))

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = isotope,
      y = mid
    ) +
    ggplot2::facet_grid(
      tracer ~ metabolite,
      labeller = ggplot2::label_parsed,
      # switch = "y",
      scales = "free_x",
      space = "free_x"
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
        x = -Inf,
        y = Inf,
        vjust = 1.5,
        hjust = -0.3,
        label = annot
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1.1)
    ) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_lf_mids <- function(df) {
  p <-
    df %>%
    dplyr::filter(cell_type == "lf") %>%
    dplyr::mutate(
      metabolite = factor(
        metabolite,
        levels = c("ALA", "ASP", "GLU", "GLN", "LAC", "SER")
      )
    ) %>%
    dplyr::filter(time == 72 & !is.na(metabolite)) %>%
    plot_mids() +
    theme_patchwork(widths = unit(7, "in"), heights = unit(7.5, "in"), tags = NULL)
}

plot_pasmc_mids <- function(df) {
  df %>%
    dplyr::filter(cell_type == "pasmc") %>%
    dplyr::mutate(
      metabolite = factor(
        metabolite,
        levels = c(
          "FBP",
          "3PG",
          "PYR",
          "ALA",
          "SER",
          "LAC",
          "CIT",
          "AKG",
          "GLN",
          "GLU",
          "MAL",
          "ASP"
        )
      ),
      metabolite = forcats::fct_recode(metabolite, "`3PG`" = "3PG")
    ) %>%
    dplyr::filter(time == 36 & !is.na(metabolite)) %>%
    plot_mids() +
    theme_patchwork(widths = unit(9.5, "in"), heights = unit(7.5, "in"), tags = NULL)
}

get_m5_citrate <- function(df) {
  x <-
    df %>%
    dplyr::filter(tracer == "q5" & treatment == "None" & isotope == "M5" & metabolite == "CIT") %>%
    dplyr::filter(cell_type == "pasmc" & time == 36)

  out <-
    x %>%
    dplyr::group_by(oxygen) %>%
    dplyr::summarise(mean = mean(mid), se = sd(mid)/sqrt(dplyr::n())) %>%
    tidyr::pivot_wider(names_from = oxygen, values_from = c(mean, se)) %>%
    unlist() * 100

  pval <- t.test(mid ~ oxygen, data = x, paired = TRUE)$p.value

  c(out, pval = pval) %>%
    round(digits = 2)
}

format_time_course_mids <- function(model_mids) {
  tracer_labels <-
    c(expression(paste("[1,2-"^13, "C"[2], "] glucose")),
      expression(paste("[U-"^13, "C"[6], "] glucose")),
      expression(paste("[U-"^13, "C"[5], "] glutamine")),
      expression(paste("[U-"^13, "C"[3], "] lactate")))

  tracer_levels <-
    c("glc2", "glc6", "q5", "lac3")

  model_mids %>%
    tidyr::unnest(c(data)) %>%
    dplyr::filter(metabolite %in% c("PYR", "CIT", "MAL")) %>%
    dplyr::filter(tracer != "lac3") %>%
    dplyr::filter(isotope != "M6") %>%
    dplyr::mutate(
      tracer = factor(
        tracer,
        levels = tracer_levels,
        labels = tracer_labels
      ),
      metabolite = factor(metabolite, levels = c("PYR", "CIT", "MAL"))
    )
}

plot_mid_time_course <- function(time_course_mids, cells, o2, treat, color) {
  time_course_mids %>%
    dplyr::filter(cell_type == cells & oxygen == o2 & treatment == treat) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = mean,
      color = isotope,
      fill = isotope
    ) +
    ggplot2::facet_grid(metabolite ~ tracer, labeller = label_parsed) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = mean - se,
        ymax = mean + se
      ),
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      pch = 21,
      color = "white",
      size = 2,
      show.legend = TRUE
    ) +
    ggplot2::labs(
      x = "Time (h)",
      y = "Mole fraction",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, 24)) +
    ggplot2::scale_color_viridis_d(option = color, end = 0.9) +
    ggplot2::scale_fill_viridis_d(option = color, end = 0.9) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots() +
    ggplot2::theme(
      legend.key.size = ggplot2::unit(0.5, units = "lines")
    )
}

plot_lactate_mids <- function(pruned_mids, cell) {
  pruned_mids %>%
    dplyr::filter(cell_type == cell) %>%
    dplyr::mutate(
      metabolite = factor(
        metabolite,
        levels = c("FBP", "3PG", "PYR", "CIT", "AKG", "MAL"),
        labels = c("FBP", "3PG", "PYR", "CIT", "AKG", "MAL")
      )
    ) %>%
    dplyr::filter(tracer == "lac3" & time == 72 & !is.na(metabolite)) %>%
    plot_mids() +
    ggplot2::facet_wrap(~ metabolite, scales = "free_x", nrow = 2) +
    ggplot2::theme(legend.position = "bottom") +
    theme_patchwork(
      widths = ggplot2::unit(4, "in"),
      heights = ggplot2::unit(2.5, "in"),
      tags = NULL
    )
}
