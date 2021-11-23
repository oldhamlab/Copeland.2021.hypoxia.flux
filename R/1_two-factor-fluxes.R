# two-factor-fluxes.R

plot_twoby_fluxes <- function(df, annot, metab, ylab) {
  df %>%
    dplyr::filter(metabolite == metab) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = oxygen,
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        fill = treatment
      ),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        group = treatment
      ),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, metabolite == metab),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      color = "black",
      size = 6/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      color = NULL,
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    # ggplot2::coord_cartesian(ylim = c(-750, 1250)) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

analyze_twoby_fluxes <- function(growth_rates, fluxes)
{
  df <-
    growth_rates %>%
    dplyr::filter(experiment == "05-bay") %>%
    dplyr::rename(flux = mu) %>%
    dplyr::select(-X0) %>%
    dplyr::mutate(metabolite = "growth") %>%
    dplyr::bind_rows(dplyr::filter(fluxes, experiment == "05-bay"))

  annot <-
    df %>%
    dplyr::group_by(metabolite) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(flux ~ oxygen * treatment + (1 | date), data = .x)),
      res = purrr::map(m, ~emmeans::emmeans(
        .x,
        "pairwise" ~ oxygen * treatment,
        simple = "each",
        adjust = "Tukey",
        combine = TRUE
      )[["contrasts"]]
      ),
      out = purrr::map(res, broom::tidy)
    ) %>%
    tidyr::unnest(c(out)) %>%
    dplyr::filter(oxygen != ".") %>%
    dplyr::select(metabolite, oxygen, adj.p.value) %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      y_pos = Inf,
      vjust = 1.5,
      lab = annot_p(adj.p.value)
    )

  list(data = df, annot = annot)
}

annot_twoby_densities <- function(blot_norm)
{
  blot_norm %>%
    dplyr::filter(experiment == "lf_05-bay") %>%
    dplyr::group_by(protein, oxygen, treatment) %>%
    wmo::remove_nested_outliers(fold_change, remove = TRUE) %>%
    dplyr::group_by(protein) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(fold_change ~ oxygen * treatment + (1 | gel), data = .x)),
      res = purrr::map(m, ~emmeans::emmeans(
        .x,
        "pairwise" ~ oxygen * treatment,
        simple = "each",
        adjust = "mvt",
        combine = TRUE
      )[["contrasts"]]
      ),
      out = purrr::map(res, broom::tidy)
    ) %>%
    tidyr::unnest(c(out)) %>%
    # dplyr::filter(treatment != ".") %>%
    dplyr::select(protein, oxygen, treatment, adj.p.value) %>%
    dplyr::mutate(
      oxygen = replace(oxygen, oxygen == ".", "0.5%"),
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      y_pos = Inf,
      vjust = 1.5,
      lab = dplyr::case_when(
        adj.p.value < 0.05 ~ "*",
        # adj.p.value < 0.07 ~ as.character(round(adj.p.value, 2)),
        TRUE ~ NA_character_
      )
    )
}

plot_twoby_densities <- function(df, prot, annot, ylab) {

  df %>%
    dplyr::filter(experiment == "lf_05-bay") %>%
    dplyr::filter(protein == prot) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = oxygen,
      y = fold_change
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        fill = treatment
      ),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        group = treatment,
        color = treatment
      ),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, protein == prot & is.na(treatment)),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        color = treatment,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt,
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, protein == prot & !is.na(treatment)),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        color = treatment,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      color = "Treatment"
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
    theme_plots()
}
