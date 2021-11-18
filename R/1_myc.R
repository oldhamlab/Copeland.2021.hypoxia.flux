# 1_myc.R

combine_fluxes <- function(growth_rates, fluxes)
{
  growth_rates %>%
    dplyr::filter(experiment == "05-simyc") %>%
    dplyr::rename(flux = mu) %>%
    dplyr::select(-X0) %>%
    dplyr::mutate(metabolite = "growth") %>%
    dplyr::bind_rows(dplyr::filter(fluxes, experiment == "05-simyc")) %>%
    dplyr::filter(treatment %in% c("siCTL", "siMYC"))
}

annot_fluxes <- function(df)
{
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
    dplyr::filter(treatment != ".") %>%
    dplyr::select(metabolite, treatment, adj.p.value) %>%
    dplyr::mutate(
      treatment = factor(treatment, levels = c("siCTL", "siMYC")),
      y_pos = Inf,
      vjust = 1.5,
      lab = annot_p(adj.p.value)
    )
}

plot_myc <- function(df, annot, metab, ylab)
{
  df %>%
    dplyr::filter(metabolite == metab) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = treatment,
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        fill = oxygen
      ),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        group = oxygen
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
        x = treatment,
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
      x = "Treatment",
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
