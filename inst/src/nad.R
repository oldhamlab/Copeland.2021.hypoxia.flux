# nad.R

# clean_nad ---------------------------------------------------------------

clean_nad <- function(nad_data) {
  nad_data %>%
    dplyr::bind_rows(.id = "metabolite") %>%
    dplyr::rename(rep = replicate) %>%
    clean_technical_replicates() %>%
    tidyr::separate(.data$experiment, c(NA, "date"), "_")
}

# finalize_nad ------------------------------------------------------------

finalize_nad <- function(nad_interp, cells_per_dna) {

  x <-
    dplyr::filter(cells_per_dna, cell_type == "lf" & volume == 200) %>%
    dplyr::pull(slope)

  nad_interp %>%
    dplyr::group_by(metabolite, date, oxygen, treatment, nucleotide) %>%
    dplyr::summarise(conc = mean(conc)) %>%
    dplyr::mutate(
      conc = dplyr::case_when(
        metabolite == "dna" ~ conc * x,
        metabolite == "nad" ~ conc * 720
      ),
      nucleotide = replace(nucleotide, is.na(nucleotide), "Count")
    ) %>%
    tidyr::pivot_wider(-c(metabolite), names_from = nucleotide, values_from = conc) %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("none", "DMSO", "BAY")),
      treatment = forcats::fct_recode(treatment, "None" = "none"),
      Ratio = NADH/NAD,
      dplyr::across(tidyselect::contains("NAD"), ~. / Count * 1000)
    ) %>%
    dplyr::select(-Count) %>%
    tidyr::pivot_longer(c(NAD, NADH, Ratio), names_to = "measurement", values_to = "value") %>%
    dplyr::arrange(measurement, oxygen, treatment)

}

# annot_nad ---------------------------------------------------------------

annot_nad <- function(nad_final) {
  nad_final %>%
    dplyr::filter(treatment != "None") %>%
    dplyr::group_by(measurement) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(value ~ treatment * oxygen + (1 | date), data = .x)),
      s = purrr::map(
        m,
        ~emmeans::emmeans(.x, ~ treatment * oxygen) %>%
          pairs(simple = "each", combine = TRUE) %>%
          broom::tidy()
      )
    ) %>%
    tidyr::unnest(c(s)) %>%
    dplyr::select(measurement, oxygen, treatment, pval = adj.p.value) %>%
    dplyr::mutate(
      lab = annot_p(pval),
      y = Inf,
      vjust = 1.5
    )
}

# plot_nad ----------------------------------------------------------------

plot_nad <- function(nad_final, annot, metab, ylab) {

  annot <- dplyr::filter(annot, measurement == metab)

  annot1 <-
    dplyr::filter(annot, oxygen != ".") %>%
    dplyr::mutate(
      treatment = "BAY",
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )
  annot2 <-
    dplyr::filter(annot, treatment != ".") %>%
    dplyr::mutate(
      oxygen = "21%",
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )

  nad_final %>%
    dplyr::filter(treatment != "None" & measurement == metab) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = treatment,
      y = value,
      fill = oxygen
    ) +
    ggplot2::stat_summary(
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = oxygen),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = annot1,
      ggplot2::aes(
        color = oxygen,
        y = y,
        label = lab,
        vjust = vjust
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot2,
      ggplot2::aes(
        y = y,
        label = lab,
        vjust = vjust
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt,
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::labs(
      x = "Treatment",
      y = ylab,
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
    # ggplot2::coord_cartesian(ylim = c(-750, 1250)) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )

}
