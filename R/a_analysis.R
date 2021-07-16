# analysis.R

read_data <- function(data_files) {
  data_files[stringr::str_detect(data_files, "\\.csv$")] %>%
    rlang::set_names(stringr::str_extract(., "(lf|pasmc)_(02|05-bay|05|bay|)")) %>%
    purrr::map_dfr(read_csv, .id = "experiment") %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%", "0.2%"), ordered = TRUE),
      treatment = factor(treatment, levels = c("None", "DMSO", "BAY"), ordered = TRUE)
    )
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

normalize_fluxes <- function(map_fluxes, reference_flux) {

  map_fluxes %>%
    dplyr::group_by(.data$cell_type, .data$treatment) %>%
    dplyr::mutate(
      lb = lb / flux[id == reference_flux],
      ub = ub / flux[id == reference_flux],
      flux = flux / flux[id == reference_flux]
    )

}

assemble_flux_differences <- function(map_fluxes) {

  glc_norm <- normalize_fluxes(map_fluxes, "GLUT")
  growth_norm <- normalize_fluxes(map_fluxes, "BIOMASS")

  tibble::tibble(
    map_fluxes = list(map_fluxes, map_fluxes, map_fluxes, glc_norm, glc_norm, glc_norm, growth_norm, growth_norm, growth_norm),
    cell = rep(c("lf", "lf", "pasmc"), 3),
    control = rep(c("21%", "DMSO", "21%"), 3),
    experiment = rep(c("0.5%", "BAY", "0.5%"), 3)
  ) %>%
    purrr::pmap_dfr(calculate_flux_differences, .id = "normalization") %>%
    dplyr::mutate(normalization = dplyr::case_when(
      normalization %in% 1:3 ~ "none",
      normalization %in% 4:6 ~ "glucose",
      normalization %in% 7:9 ~ "growth"
    ))

}

make_grid <- function(arrow, lhs, rhs) {
  forward <- tidyr::expand_grid(substrate = lhs, product = rhs)
  reverse <- tidyr::expand_grid(substrate = rhs, product = lhs)

  if (arrow == "->") {
    forward
  } else {
    reverse
  }
}

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

make_graph <- function(map_flux_differences, cell, nodes, treat, normalizer) {

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
    dplyr::filter(cell_type == cell & normalization == normalizer & treatment == treat) %>%
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
      y = NULL,
      title = caption
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

plot_normoxia_network <- function(hypoxia_graph, cell) {

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
