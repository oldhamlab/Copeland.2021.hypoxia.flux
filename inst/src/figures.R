# figures.R

# theme_plots -------------------------------------------------------------

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

# theme_patchwork ---------------------------------------------------------

theme_patchwork <- function(design = NULL, widths = NULL, heights = NULL, tags = "A") {
  list(
    patchwork::plot_layout(
      design = design,
      widths = widths,
      heights = heights
    ),
    patchwork::plot_annotation(
      tag_levels = tags,
      theme = ggplot2::theme(plot.margin = ggplot2::margin(-3, -3, -3, -3))
    )
  )
}

# write_figures -----------------------------------------------------------

write_figures <- function(plot, filename) {
  path <-
    system.file(
      "manuscript/figures",
      package = "Copeland.2021.hypoxia.flux"
    )

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

# arrange_fluxes ----------------------------------------------------------

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
      heights = unit(1.25, "in")
    )
}

# arrange_s1 --------------------------------------------------------------

arrange_s1 <- function(a, b, c, d) {
  layout <- "
  aa#
  bc#
  ddd
  "

  a + b + c + d +
    theme_patchwork(
      design = layout,
      widths = unit(1.25, "in"),
      heights = unit(1.25, "in")
    )
}

# arrange_m3 --------------------------------------------------------------

arrange_m3 <- function(m3ab, m3c) {

  a <- m3ab$curve
  b <- m3ab$rate
  c <- m3c

  layout <- "
  aab
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

# arrange_m4 --------------------------------------------------------------

arrange_m4 <- function(hypoxia_graph_ratio_plot, bay_graph_ratio_plot, m4c) {

  a <- hypoxia_graph_ratio_plot
  b <- bay_graph_ratio_plot

  top <-
    a + b +
    plot_layout(guides = "collect") &
    theme(legend.text.align = 1)

  top / m4c +
    theme_patchwork(widths = unit(5, "in"), heights = unit(c(2.5, 2), "in"))
}

# arrange_s6 --------------------------------------------------------------

arrange_s6 <- function(a, b, c, d) {
  layout <- "
    ab
    cd
    cd"

  a + b + c + d +
    theme_patchwork(
      design = layout,
      widths = unit(5, "in"),
      heights = unit(c(2.5), "in")
    )
}

# arrange_m5 --------------------------------------------------------------

arrange_m5 <- function(a, b, c, d, e, f, g, h, i) {
  layout <- "
  abc
  deg
  ffg
  hhi
  "
  a + b + c + d + e + f + g + h + i +
    theme_patchwork(
      design = layout,
      widths = unit(1.25, "in"),
      heights = unit(c(1.25), "in")
    )
}
