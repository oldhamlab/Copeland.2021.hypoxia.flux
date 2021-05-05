# analysis.R

# plot_growth_curve -------------------------------------------------------

plot_growth_curve <- function(
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  clr = c("oxygen", "treatment")
) {

  flux_measurements %>%
    dplyr::filter(
      cell_type == cell &
        experiment == exp &
        metabolite == "cells" &
        time < 96
    ) %>%
    dplyr::group_by(
      date,
      oxygen,
      treatment,
      time
    ) %>%
    dplyr::summarize(count = mean(conc, na.rm = TRUE)) %>%
    plot_time_lines(
      y = count/1000,
      ylab = expression(paste("Cell count (x", 10^3, ")")),
      clr = clr
    )
}

# plot_growth_rates -------------------------------------------------------

plot_growth_rates <- function(
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  clr = c("oxygen", "treatment"),
  xaxis = clr,
  annot
) {

  if(missing(annot)) {
    annot <-
      tibble::tibble(
        p = t.test(
          mu ~ !!rlang::ensym(clr),
          data = growth_rates,
          subset = c(cell_type == cell & experiment == exp),
          paired = TRUE
        )$p.value,
        y_pos = Inf,
        x_pos = 1.5,
        vjustvar = 1.5) %>%
      dplyr::mutate(lab = annot_p(p))
  }

  growth_rates %>%
    dplyr::filter(cell_type == cell & experiment == exp) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = .data[[xaxis]],
      y = mu
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = .data[[clr]]),
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
    # ggbeeswarm::geom_beeswarm(
    #   ggplot2::aes(color = .data[[clr]]),
    #   priority = "random",
    #   alpha = 0.4,
    #   cex = 4,
    #   size = 1,
    #   dodge.width = 0.7,
    #   show.legend = FALSE
    # ) +
    # ggplot2::stat_summary(
    #   color = "black",
  #   ggplot2::aes(group = .data[[clr]]),
  #   geom = "crossbar",
  #   fun = "mean",
  #   fatten = 1,
  #   width = 0.4,
  #   position = ggplot2::position_dodge(width = 0.7),
  #   show.legend = FALSE
  # ) +
  ggplot2::geom_text(
    data = annot,
    ggplot2::aes(
      label = lab,
      x = x_pos,
      y = y_pos,
      vjust = vjustvar
    ),
    family = "Calibri",
    size = 6/ggplot2::.pt
  ) +
    ggplot2::labs(
      x = "Treatment",
      y = "Growth rate (/h)"
    ) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::coord_cartesian(ylim = c(0, 0.03)) +
    theme_plots()
}

# clean_viability ---------------------------------------------------------

clean_viability <- function(viability_file) {
  readr::read_csv(viability_file) %>%
    dplyr::mutate(
      viability = 100 * live / (dead + live),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) %>%
    dplyr::group_by(time, oxygen) %>%
    wmo::remove_nested_outliers(viability, remove = TRUE)
}

# plot_blot ---------------------------------------------------------------

plot_blot <- function(blot_image) {

  # blot_image <- magick::image_read_pdf(blot_image)
  blot_image <- magick::image_read(blot_image)

  cowplot::ggdraw() +
    cowplot::draw_image(
      blot_image,
      scale = 1.3,
      hjust = 0.15,
      vjust = 0.1
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


# plot_densities ----------------------------------------------------------

plot_densities <- function(
  densities,
  exp,
  prot = c("ldha", "hif1a"),
  ylab = prot,
  clr = c("oxygen", "treatment")
) {

  densities %>%
    dplyr::filter(experiment == exp) %>%
    dplyr::filter(protein == prot) %>%
    plot_time_lines(y = fold_change, ylab = ylab, clr = clr)
}

# normalize_qpcr ----------------------------------------------------------

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
      names_to = "rna",
      values_to = "dct"
    ) %>%
    dplyr::filter(!is.na(dct)) %>%
    dplyr::group_by(.data$experiment, .data$rna) %>%
    dplyr::mutate(
      ddct = dct - mean(
        dct[oxygen == min(oxygen) & treatment %in% c("None", "DMSO") & time == 0]
      ),
      log2fc = 2 ^ -ddct
    ) %>%
    dplyr::group_by(experiment, oxygen, treatment, time, rna) %>%
    wmo::remove_nested_outliers(log2fc, remove = TRUE)
}

# plot_mrna ---------------------------------------------------------------

plot_mrna <- function(
  mrna,
  exp,
  gene = c("ldha", "hif1a"),
  ylab = prot,
  clr = c("oxygen", "treatment")
) {

  mrna %>%
    dplyr::filter(rna == gene) %>%
    dplyr::filter(experiment == exp) %>%
    plot_time_lines(y = log2fc, ylab = ylab, clr = clr)
}

# plot_high_fluxes --------------------------------------------------------

plot_high_fluxes <- function(
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  clr = c("oxygen", "treatment"),
  xaxis = clr,
  annot
) {

  annot <-
    annot %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
        metabolite %in% c("lactate", "glucose")
    )

  fluxes %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
        metabolite %in% c("lactate", "glucose")
    ) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = .data[[clr]]),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = .data[[clr]]),
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
        y = y_pos,
        vjust = vjust,
        label = pval
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Flux (fmol/cell/h)"
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    theme_plots()
}

# plot_low_fluxes ---------------------------------------------------------

plot_low_fluxes <- function(
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  clr = c("oxygen", "treatment"),
  xaxis = clr,
  annot
) {

  annot <-
    annot %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
        !(metabolite %in% c("lactate", "glucose"))
    )

  fluxes %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
        !(metabolite %in% c("lactate", "glucose"))
    ) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = .data[[clr]]),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = .data[[clr]]),
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
        y = y_pos,
        vjust = vjust,
        label = pval
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
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

# plot_viability ----------------------------------------------------------

plot_viability <- function(viability) {

  plot_time_lines(
    viability,
    y = viability,
    ylab = "Cell viability (%)",
    clr = "oxygen"
  )
}

# read_data ---------------------------------------------------------------

read_data <- function(data_files) {
  data_files[stringr::str_detect(data_files, "\\.csv$")] %>%
    rlang::set_names(stringr::str_extract(., "(lf|pasmc)_(02|05-bay|05|bay|)")) %>%
    purrr::map_dfr(read_csv, .id = "experiment") %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%", "0.2%"), ordered = TRUE),
      treatment = factor(treatment, levels = c("None", "DMSO", "BAY"), ordered = TRUE)
    )
}

# normalize_densities -----------------------------------------------------

normalize_densities <- function(blot_raw) {
  blot_raw %>%
    dplyr::filter(.data$time < 96) %>%
    dplyr::filter(!(experiment == "lf_05-bay" & gel == "b")) %>%
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
        mean(density[oxygen == min(oxygen) & treatment %in% c("None", "DMSO") & time == min(time)])
    ) %>%
    dplyr::group_by(experiment, oxygen, treatment, time, protein) %>%
    wmo::remove_nested_outliers(fold_change, remove = TRUE)
}

# plot_time_lines ---------------------------------------------------------

plot_time_lines <- function(
  df,
  y,
  ylab = prot,
  clr = c("oxygen", "treatment")
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
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots()
}

# annot_fluxes ------------------------------------------------------------

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
    ) #%>%
  # tidyr::unnest(c(data)) %>%
  # dplyr::group_by(
  #   dplyr::across(c(metabolite:experiment, oxygen, treatment, ttest:vjust))
  #   ) %>%
  # dplyr::summarize(flux = mean(flux))
}

# plot_cells_per_dna ------------------------------------------------------

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
    ggplot2::coord_cartesian(xlim = c(0, NA)) +
    theme_plots()
}

# clean_dna_count_hypoxia -------------------------------------------------

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

# plot_dna_count_hypoxia --------------------------------------------------

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
    ggplot2::coord_cartesian(ylim = c(0, NA), xlim = c(0, 300)) +
    theme_plots()
}

# plot_evap_data ----------------------------------------------------------

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

# plot_k ------------------------------------------------------------------

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
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
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
      values = clrs
    ) +
    ggplot2::scale_fill_manual(
      values = clrs
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

# plot_mids ---------------------------------------------------------------

plot_mids <- function(df) {

  tracer_labels <-
    c(expression(paste("[1,2-"^13, "C"[2], "] glucose")),
      expression(paste("[U-"^13, "C"[6], "] glucose")),
      expression(paste("[U-"^13, "C"[5], "] glutamine")),
      expression(paste("[U-"^13, "C"[3], "] lactate")))

  tracer_levels <-
    c("glc2", "glc6", "q5", "lac3")

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
    ) %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = isotope, y = mean, fill = group) +
    ggplot2::facet_grid(
      tracer ~ metabolite,
      labeller = ggplot2::label_parsed,
      # switch = "y",
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge()
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = mean - se,
        ymax = mean + se
      ),
      position = ggplot2::position_dodge(width = 0.9)
    ) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL
    ) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

# plot_manuscript_mids ----------------------------------------------------

plot_manuscript_mids <- function(model_mids) {
  model_mids %>%
    dplyr::filter(cell_type == "lf") %>%
    tidyr::unnest(c(data)) %>%
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

# plot_lf_mids ------------------------------------------------------------

plot_lf_mids <- function(model_mids) {
  p <-
    model_mids %>%
    dplyr::filter(cell_type == "lf") %>%
    tidyr::unnest(c(data)) %>%
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

# plot_pasmc_mids ---------------------------------------------------------

plot_pasmc_mids <- function(model_mids) {
  model_mids %>%
    dplyr::filter(cell_type == "pasmc") %>%
    tidyr::unnest(c(data)) %>%
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
    dplyr::filter(time == 48 & !is.na(metabolite)) %>%
    plot_mids() +
    # ggplot2::facet_grid(
    #   metabolite ~ tracer,
    #   labeller = ggplot2::label_parsed
    # ) +
    theme_patchwork(widths = unit(9.5, "in"), heights = unit(7.5, "in"), tags = NULL)
}

# clean_model_fluxes ------------------------------------------------------

clean_model_fluxes <- function(map_flux_files, model_reactions) {
  pathways <-
    model_reactions %>%
    dplyr::select(-equation) %>%
    dplyr::add_row(name = "BIOMASS", pathway = factor("biomass"), index = 0)

  map_flux_files[stringr::str_detect(map_flux_files, "(lf|pasmc)_.*_model\\.csv")] %>%
    rlang::set_names(stringr::str_extract(basename(.), pattern = ".*(?=_model\\.csv)")) %>%
    purrr::map_dfr(readr::read_csv, .id = "experiment") %>%
    tidyr::separate(
      experiment,
      c("cell_type", "treatment"),
      sep = "_"
    ) %>%
    dplyr::mutate(
      treatment = dplyr::case_when(
        treatment == "21" ~ "21%",
        treatment == "05" ~ "0.5%",
        treatment == "dmso" ~ "DMSO",
        treatment == "bay" ~ "BAY"
      ),
      treatment = factor(treatment, levels = c("21%", "0.5%", "DMSO", "BAY"))
    ) %>%

    # identify net and exchange fluxes
    tidyr::separate(id, c("id", "type"), sep = " ", fill = "right") %>%
    dplyr::mutate(type = replace(type, is.na(type), "net")) %>%
    dplyr::select(-se) %>%

    # add pathway info
    dplyr::left_join(pathways, by = c("id" = "name")) %>%
    dplyr::select(
      cell_type,
      treatment,
      pathway,
      index,
      id,
      type,
      equation,
      flux = value,
      lb,
      ub
    ) %>%
    dplyr::group_by(cell_type, treatment) %>%
    dplyr::arrange(treatment, pathway)
}

# plot_labeling_rate ------------------------------------------------------

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

# calculate_flux_differences ----------------------------------------------

calculate_flux_differences <- function(map_fluxes, cell, control, experiment) {
  map_fluxes %>%
    dplyr::filter(cell_type == cell) %>%
    dplyr::filter(treatment %in% c(control, experiment)) %>%
    dplyr::mutate(
      treatment = dplyr::case_when(
        treatment == control ~ "ctl",
        treatment == experiment ~ "exp"
      )
    ) %>%
    tidyr::pivot_wider(names_from = treatment, values_from = c(flux, lb, ub)) %>%
    dplyr::mutate(
      ratio = dplyr::if_else(
        pmax(lb_ctl, lb_exp) - pmin(ub_ctl, ub_exp) > 0,
        flux_exp/flux_ctl,
        NA_real_
      ),
      ctl = control,
      exp = experiment
    ) %>%
    dplyr::select(
      cell_type:equation,
      ctl,
      exp,
      tidyselect::contains("ctl"),
      tidyselect::contains("exp"),
      ratio
    )
}

# normalize_fluxes --------------------------------------------------------

normalize_fluxes <- function(map_fluxes, reference_flux) {

  map_fluxes %>%
    dplyr::group_by(.data$cell_type, .data$treatment) %>%
    dplyr::mutate(
      lb = lb / flux[id == reference_flux],
      ub = ub / flux[id == reference_flux],
      flux = flux / flux[id == reference_flux]
    )

}

# assemble_flux_differences -----------------------------------------------

assemble_flux_differences <- function(map_fluxes) {

  glc_norm <- normalize_fluxes(map_fluxes, "GLUT")
  growth_norm <- normalize_fluxes(map_fluxes, "BIOMASS")

  tibble::tibble(
    map_fluxes = list(map_fluxes, map_fluxes, glc_norm, glc_norm, growth_norm, growth_norm),
    cell = rep(c("lf", "lf"), 3),
    control = rep(c("21%", "DMSO"), 3),
    experiment = rep(c("0.5%", "BAY"), 3)
  ) %>%
    purrr::pmap_dfr(calculate_flux_differences, .id = "normalization") %>%
    dplyr::mutate(normalization = dplyr::case_when(
      normalization %in% 1:2 ~ "none",
      normalization %in% 3:4 ~ "glucose",
      normalization %in% 5:6 ~ "growth"
    ))

}

# make_grid ---------------------------------------------------------------

make_grid <- function(arrow, lhs, rhs) {
  forward <- tidyr::expand_grid(substrate = lhs, product = rhs)
  reverse <- tidyr::expand_grid(substrate = rhs, product = lhs)

  if (arrow == "->") {
    forward
  } else {
    reverse
  }
}

# parse_eq ----------------------------------------------------------------

parse_eq <- function(df) {

  pattern <- "\\w?[^0-9*\\s\\.]\\w+(\\.[a-z]*)?"

  df %>%

    # select relevant fluxes
    dplyr::filter(type == "net" & !stringr::str_detect(equation, ".ms")) %>%

    # parse equation
    dplyr::mutate(
      arrow = stringr::str_extract(equation, pattern = "(<->)|(<-)|(->)"),
      halved = stringr::str_split(equation, pattern = arrow),
      lhs = purrr::map(
        halved,
        ~ unlist(stringr::str_extract_all(.x[1], pattern = pattern))
      ),
      rhs = purrr::map(
        halved,
        ~ unlist(stringr::str_extract_all(.x[2], pattern = pattern))
      )
    ) %>%

    # deal with directionality
    dplyr::mutate(
      arrow = dplyr::case_when(flux >= 0 ~ "->", flux < 0 ~ "<-"),
      lb = dplyr::case_when(flux >= 0 ~ lb, flux < 0 ~ -ub),
      ub = dplyr::case_when(flux >= 0 ~ ub, flux < 0 ~ -lb),
      flux = abs(flux),
      grid = purrr::pmap(list(arrow, lhs, rhs), make_grid)
    ) %>%
    tidyr::unnest(c(grid)) %>%

    # eliminate duplicated fluxes
    dplyr::distinct()
}

# make_graph --------------------------------------------------------------

make_graph <- function(map_flux_differences, nodes, treat, normalizer) {

  edges <-
    map_flux_differences %>%
    dplyr::rename(
      treatment_ctl = ctl,
      treatment_exp = exp,
    ) %>%
    tidyr::pivot_longer(
      tidyselect::matches("_(ctl|exp)"),
      names_to = c(".value", NA),
      names_sep = "_"
    ) %>%
    dplyr::filter(normalization == normalizer & treatment == treat) %>%
    parse_eq() %>%
    dplyr::select(
      cell_type,
      treatment,
      from = substrate,
      to = product,
      flux,
      ratio
    ) %>%
    dplyr::mutate(ratio = log(ratio, base = 2)) %>%
    dplyr::distinct()

  tidygraph::tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = TRUE
  )

}

# set_graph_style ---------------------------------------------------------

set_gr_style <- function() {
  ggraph::set_graph_style(
    family = "Calibri",
    face = "plain",
    size = 8,
    text_size = 6,
    text_colour = "black",
    title_face = "plain",
    plot_margin = ggplot2::margin(5, 5, 5, 5)
  )
}

# plot_ratio_network ------------------------------------------------------

plot_ratio_network <- function(graph, caption) {
  fold <- c(Inf, 5, 3, 2, 1.5, 1.1, 0.91, 0.67, 0.5, 0.33, 0.2, 0)

  fold_cuts <- log(fold, base = 2)

  lvls <-
    tidygraph::activate(graph, edges) %>%
    dplyr::pull(ratio) %>%
    cut(
      fold_cuts,
      include.lowest = TRUE
    )

  clrs <-
    RColorBrewer::brewer.pal(11, "PRGn") %>%
    rlang::set_names(levels(lvls))

  clrs[ceiling(length(clrs)/2)] <- "grey90"

  labs <-
    c(
      "> +5x",
      "+3x",
      "+2x",
      "+1.5x",
      "+1.1x",
      "0x",
      "-1.1x",
      "-1.5x",
      "-2x",
      "-3x",
      "< -5x"
    ) %>%
    rev() %>%
    rlang::set_names(names(clrs))

  graph <-
    graph %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(Flux = dplyr::case_when(
      flux > 1000 ~ "> 1000",
      flux > 100 ~ "> 100",
      flux > 10 ~ "> 10",
      TRUE ~ "> 0"
    ))

  set_gr_style()

  ggraph::ggraph(
    graph,
    circular = FALSE,
    x = tidygraph::.N()$x,
    y = tidygraph::.N()$y
  ) +
    ggraph::geom_edge_parallel(
      ggplot2::aes(
        color = cut(ratio, fold_cuts, include.lowest = TRUE),
        start_cap = ggraph::label_rect(
          node1.name,
          fontfamily = "Calibri",
          fontsize = 6,
          padding = ggplot2::margin(1, 1, 1, 1, "mm")
        ),
        end_cap = ggraph::label_rect(
          node2.name,
          fontfamily = "Calibri",
          fontsize = 6,
          padding = ggplot2::margin(1, 1, 1, 1, "mm")
        ),
        edge_width = Flux
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type = "open"),
      linejoin = "mitre",
      lineend = "round",
      show.legend = TRUE
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      size = 6/ggplot2::.pt
    ) +
    ggraph::scale_edge_color_manual(
      name = "Ratio",
      values = clrs,
      labels = labs,
      na.translate = TRUE,
      na.value = "grey90",
      drop = FALSE
    ) +
    ggraph::scale_edge_width_manual(
      values = c(0.25, 0.5, 1, 1.5),
      guide = "legend"
      ) +
    ggplot2::labs(
      x = NULL,
      y = NULL
      # title = caption
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.position = c(0.07, 0.2),
      plot.title = ggplot2::element_text(
        size = ggplot2::rel(1),
        hjust = 0.5,
        face = "bold",
        margin = ggplot2::margin(b = 0)
      ),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      legend.text.align = 1
    ) +
    NULL

}

# plot_lactate_mids -------------------------------------------------------

plot_lactate_mids <- function(model_mids, cell) {
  model_mids %>%
    dplyr::filter(cell_type == cell) %>%
    tidyr::unnest(c(data)) %>%
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
    ggplot2::theme(legend.position = "right") +
    theme_patchwork(widths = unit(4, "in"), heights = unit(2.5, "in"), tags = NULL)
}

# format_time_course_mids -------------------------------------------------

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

# plot_mid_time_course ----------------------------------------------------

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

# plot_normoxia_network ---------------------------------------------------

plot_normoxia_network <- function(hypoxia_graph) {

  hypoxia_graph %>%
    tidygraph::activate(edges) %>%
    dplyr::filter(flux > 0.1) %>%
    ggraph::ggraph(
      circular = FALSE,
      x = tidygraph::.N()$x,
      y = tidygraph::.N()$y
    ) +
    ggraph::geom_edge_parallel(
      ggplot2::aes(
        # color = log(flux, base = 10),
        color = flux,
        start_cap = ggraph::label_rect(
          node1.name,
          fontfamily = "Calibri",
          fontsize = 6,
          padding = ggplot2::margin(1, 1, 1, 1, "mm")
        ),
        end_cap = ggraph::label_rect(
          node2.name,
          fontfamily = "Calibri",
          fontsize = 6,
          padding = ggplot2::margin(1, 1, 1, 1, "mm")
        )
      ),
      edge_width = 0.5,
      arrow = ggplot2::arrow(length = ggplot2::unit(1, "mm"), type = "open"),
      linejoin = "mitre",
      lineend = "round"
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      size = 6/ggplot2::.pt
    ) +
    ggraph::scale_edge_color_viridis(
      name = "Flux",
      trans = "log10",
      end = 0.95,
      # breaks = c(0.01, 0.1, 1, 10, 100, 1000),
      # labels = c(1, 10, 100, 1000),
      option = "plasma"
    ) +
    ggplot2::guides(
      edge_color = ggraph::guide_edge_colorbar(barwidth = 0.5, barheight = 3)
    ) +
    ggplot2::labs(
      title = "Normoxia",
      y = NULL
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.position = c(0.07, 0.14),
      plot.title = ggplot2::element_text(
        size = rel(1),
        hjust = 0.5,
        face = "bold",
        margin = ggplot2::margin(b = 0)
      ),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      legend.text.align = 1
    ) +
    NULL
}

# plot_twoby_fluxes -------------------------------------------------------

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
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
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

# analyze_twoby_fluxes ----------------------------------------------------

analyze_twoby_fluxes <- function(growth_rates, fluxes) {

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
      m = purrr::map(data, ~lmerTest::lmer(flux ~ oxygen * treatment + (1|date), data = .x)),
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

# annot_twoby_densities ---------------------------------------------------

annot_twoby_densities <- function(blot_norm) {
  blot_norm %>%
    dplyr::filter(experiment == "lf_05-bay") %>%
    dplyr::group_by(protein) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(fold_change ~ oxygen * treatment + (1|gel), data = .x)),
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
    # dplyr::filter(treatment != ".") %>%
    dplyr::select(protein, oxygen, treatment, adj.p.value) %>%
    dplyr::mutate(
      oxygen = replace(oxygen, oxygen == ".", "0.5%"),
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      y_pos = Inf,
      vjust = 1.5,
      lab = annot_p(adj.p.value)
    )
}


# plot_twoby_densities ----------------------------------------------------

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
