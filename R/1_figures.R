# figures.R

# plot setup
clrs <- c(RColorBrewer::brewer.pal(4, "Set1")[1:4], "#08306b", RColorBrewer::brewer.pal(9, "Set1")[9:8])
names(clrs) <- c("21%", "0.5%", "DMSO", "BAY", "0.2%", "siCTL", "siMYC")

# functions
theme_plots <- function() {
  list(
    wmo::theme_wmo(
      base_family = "Calibri",
      base_size = 8
    ),
    ggplot2::theme(
      panel.border = ggplot2::element_rect(size = 0.25),
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.title.y.left = ggplot2::element_text(margin = ggplot2::margin(r = 3))
    ),
    ggplot2::coord_cartesian(clip = "off")
  )
}

theme_patchwork <- function(design = NULL, widths = NULL, heights = NULL, tags = "A", ...) {
  list(
    patchwork::plot_layout(
      design = design,
      widths = widths,
      heights = heights,
      ...
    ),
    patchwork::plot_annotation(
      tag_levels = tags,
      theme = ggplot2::theme(plot.margin = ggplot2::margin(-3, -3, -3, -3))
    )
  )
}

write_figures <- function(plot, filename, path = "manuscript/figures") {
  # path <- "manuscript/figures"

  gtab <- patchwork::patchworkGrob(plot)

  overall_width <-
    grid::convertWidth(
      sum(gtab$widths),
      unitTo = "in",
      valueOnly = TRUE
    )

  overall_height <-
    grid::convertHeight(
      sum(gtab$heights),
      unitTo = "in",
      valueOnly = TRUE
    )

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = cairo_pdf,
    path = path,
    width = overall_width,
    height = overall_height,
    units = "in"
  )

  if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

  stringr::str_c(path, "/", filename)

}

arrange_fluxes <- function(a, b, c, d, e, f, g, h, i, j) {
  layout <- "
  abc
  def
  ghi
  jjj
"

  a + b + c + d + e + f + g + h + i + j +
    theme_patchwork(
      design = layout,
      widths = unit(1.25, "in"),
      heights = unit(1.25, "in"),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_s1 <- function(a, b, c, d, e) {
  layout <- "
  abb
  cd#
  eee
  "

  a + b + c + d + e +
    theme_patchwork(
      design = layout,
      widths = unit(1.25, "in"),
      heights = unit(1.25, "in")
    )
}

arrange_m2 <- function(m2ab, m2c) {

  a <- m2ab$curve
  b <- m2ab$rate
  c <- m2c

  layout <- "
  ab#
  ccc
  ccc
"

  a + b + c +
    theme_patchwork(
      design = layout,
      widths = unit(1.5, "in"),
      heights = unit(c(1.25), "in")
    )
}

arrange_s6 <- function(a, b) {
  layout <- "
    a
    b
  "

  a + b +
    theme_patchwork(
      design = layout,
      widths = unit(3.5, "in"),
      heights = unit(c(3), "in")
    )
}

arrange_m4 <- function(hypoxia_graph_ratio_plot, bay_graph_ratio_plot) {

  a <- hypoxia_graph_ratio_plot
  b <- bay_graph_ratio_plot

  a + b +
    theme_patchwork(widths = unit(3, "in"), heights = unit(3.5, "in"), guides = "collect")
}

arrange_s7 <- function(a, b, c, d) {
  layout <- "
    ab
    cd
  "

  a + b + c + d +
    theme_patchwork(
      design = layout,
      widths = unit(3, "in"),
      heights = unit(c(3.5), "in"),
      guides = "collect"
    )
}

arrange_m5 <- function(a, b, c, d, e, f, g, h, i, j, k) {
  ((a | b | c) + plot_layout(guides = "collect")) /
    (d | e) /
    (f + plot_layout(guides = "collect")) /
    ((g | h) + plot_layout(widths = c(1, 2), guides = "keep")) /
    ((i | j | k) + plot_layout(guides = "collect")) +
    theme_patchwork(
      # design = layout,
      widths = unit(5, "in"),
      heights = unit(c(1.1), "in"),
      # guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_s8 <- function(a, b, c, d, e) {
  (a | b) /
    (c | d | e) + plot_layout(widths = c(1, 2, 2)) +
    theme_patchwork(
      widths = unit(5, "in"),
      heights = unit(c(1.1), "in")
    )
}

arrange_m6 <- function(a, b, c, d, e)
{
  layout <- "
  abde
  ccde
  "
  a + b + c + d + e +
  theme_patchwork(
    design = layout,
    widths = unit(c(2, 2, 1.5, 1.5), "in"),
    heights = unit(2, "in"),
    # guides = "collect"
  ) &
  theme(legend.position = "bottom")
}

arrange_s9 <- function(a, b, c, d, e, f) {
  # layout <- "
  # aabb
  # cc##
  # d#e#
  # "

  (a | b) /
    (c | d) /
    (e | f) +

    # a + b + c + d + e +
    theme_patchwork(
      # design = layout,
      widths = unit(5, "in"),
      heights = unit(c(1.25, 2, 2.5), "in")
    )
}

create_resources <- function() {
  tibble::tribble(
    ~category, ~`REAGENT or RESOURCE`, ~SOURCE, ~IDENTIFIER,
    "Antibodies", "HIF-1α", "BD Biosciences", "610958",
    "Antibodies", "c-MYC", "Cell Signaling Technologies", "D84C12",
    "Antibodies", "LDHA", "Cell Signaling Technologies", "2012",
    "Antibodies", "HRP-α-Rabbit IgG", "Cell Signaling Technologies", "7074",
    "Antibodies", "HRP-α-Mouse IgG", "Cell Signaling Technologies", "7076",
    "Chemicals, peptides, and recombinant proteins", "[1,2-^13^C~1~] glucose", "Cambridge Isotope Labs", "CLM-504-PK",
    "Chemicals, peptides, and recombinant proteins", "[U-^13^C~6~] glucose", "Cambridge Isotope Labs", "CLM-1396-PK",
    "Chemicals, peptides, and recombinant proteins", "[U-^13^C~5~] glutamine", "Cambridge Isotope Labs", "CLM-1822-H-PK",
    "Chemicals, peptides, and recombinant proteins", "[U-^13^C~3~] lactate", "Sigma", "485926",
    "Chemicals, peptides, and recombinant proteins", "Molidustat (BAY-85-3934)", "Cayman", "15297",
    "Critical commercial assays", "Glucose colorimetric assay kit", "Cayman", "10009582",
    "Critical commercial assays", "ʟ-Lactate assay kit", "Cayman", "700510",
    "Critical commercial assays", "Pyruvate assay kit", "Cayman", "700470",
    "Depositied data", "Raw and analyzed data", "This paper", "https://github.com/oldhamlab/Copeland.2021.hypoxia.flux",
    "Depositied data", "RNA-seq reads", "This paper", "SRA: PRJNA721596",
    "Depositied data", "Summarized RNA-seq data", "This paper", "https://github.com/oldhamlab/rnaseq.lf.hypoxia.molidustat",
    "Experimental models: Cell lines", "Normal human lung fibroblasts", "Lonza", "CC-2512",
    "Experimental models: Cell lines", "Pulmonary artery smooth muscle cells", "Lonza", "CC-2581",
    "Oligonucleotides", "ACTB (Hs03023943_g1)", "Life Technologies", "4351370",
    "Oligonucleotides", "GLUT1 (Hs00892681_m1)", "Life Technologies", "4351370",
    "Oligonucleotides", "LDHA (Hs00855332_g1)", "Life Technologies", "4351370"
  ) %>%
    flextable::as_grouped_data(groups = c("category")) %>%
    flextable::as_flextable(hide_grouplabel = TRUE) %>%
    flextable::bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body") %>%
    flextable::bold(part = "header", bold = TRUE) %>%
    flextable::colformat_double(
      i = ~ is.na(category),
      j = "REAGENT or RESOURCE",
      digits = 0,
      big.mark = ""
    ) %>%
    flextable::compose(
      i = 8,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[1,2-", flextable::as_sup("13"), "C", flextable::as_sub("2"), "] glucose")
    ) %>%
    flextable::compose(
      i = 9,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("6"), "] glucose")
    ) %>%
    flextable::compose(
      i = 10,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("5"), "] glutamine")
    ) %>%
    flextable::compose(
      i = 11,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("3"), "] lactate")
    ) %>%
    flextable::font(fontname = "Calibri", part = "all") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::set_table_properties(layout = "autofit")
}

format_flux_table <- function(
  flux_differences,
  cell = c("lf", "pasmc"),
  experiment = c("0.5%", "BAY"),
  ssr_ctl = NULL,
  ssr_exp = NULL
) {

  big_border = flextable::fp_border_default(color = "black", width = 1)
  small_border = flextable::fp_border_default(color = "black", width = 0.25)

  conditions <- c(unique(flux_differences$ctl), unique(flux_differences$exp))

  df <-
    flux_differences %>%
    dplyr::ungroup() %>%
    dplyr::filter(normalization == "none" & cell_type == cell & exp == experiment)

  conditions <- c(unique(df$ctl), unique(df$exp))

  df %>%
    dplyr::select(-c(normalization, cell_type, index, ctl, exp)) %>%
    dplyr::arrange(desc(type)) %>%
    dplyr::mutate(
      pathway = stringr::str_to_sentence(pathway),
      type = toupper(type),
      dplyr::across(tidyselect::matches("ctl|exp"), ~scales::scientific(.x))
    ) %>%
    flextable::as_grouped_data(groups = c("type", "pathway")) %>%
    flextable::as_flextable(hide_grouplabel = TRUE) %>%
    flextable::border_remove() %>%
    flextable::bold(j = 1, i = ~ !is.na(pathway), bold = TRUE, part = "body") %>%
    flextable::bold(j = 1, i = ~ !is.na(type), bold = TRUE, part = "body") %>%
    flextable::add_header_row(
      values = c("", conditions, ""),
      colwidths = c(2, 3, 3, 1)
    ) %>%
    flextable::compose(
      i = 1,
      j = 3:5,
      part = "header",
      value = flextable::as_paragraph(conditions[[1]], flextable::as_sup("a"))
    ) %>%
    flextable::compose(
      i = 1,
      j = 6:8,
      part = "header",
      value = flextable::as_paragraph(conditions[[2]], flextable::as_sup("b"))
    ) %>%
    flextable::set_header_labels(
      id = "ID",
      equation = "Reaction",
      flux_ctl = "Flux",
      lb_ctl = "LB",
      ub_ctl = "UB",
      flux_exp = "Flux",
      lb_exp = "LB",
      ub_exp = "UB",
      ratio = "Ratio"
    ) %>%
    flextable::align(i = 1, part = "header", align = "center") %>%
    flextable::align(i = 2, j = 3:8, part = "header", align = "center") %>%
    flextable::merge_h(part = "header") %>%
    flextable::bold(part = "header", bold = TRUE) %>%
    flextable::colformat_double(
      digits = 2,
      big.mark = ""
    ) %>%
    flextable::hline_top(part = "header", border = big_border) %>%
    flextable::hline_bottom(part = "all", border = big_border) %>%
    flextable::hline(i = ~ !is.na(type), border = small_border) %>%
    flextable::hline(i = 1, j = c(3:5, 6:8), border = small_border, part = "header") %>%
    flextable::add_footer_lines(c("a", "b")) %>%
    flextable::compose(
      i = 1,
      part = "footer",
      value = flextable::as_paragraph(flextable::as_sup("a"), ssr_ctl)
    ) %>%
    flextable::compose(
      i = 2,
      part = "footer",
      value = flextable::as_paragraph(flextable::as_sup("b"), ssr_exp)
    ) %>%
    flextable::font(fontname = "Calibri", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::set_table_properties(layout = "autofit")
}
