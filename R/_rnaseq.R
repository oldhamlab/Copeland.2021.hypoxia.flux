# rnaseq.R

# count_rnaseq ------------------------------------------------------------

count_rnaseq <- function() {

  count_data <- lf.hypoxia.molidustat.rnaseq::count_data

  pheno_data <-
    lf.hypoxia.molidustat.rnaseq::pheno_data %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%"), labels = c("N", "H")),
      treatment = factor(treatment, labels = c("DMSO", "BAY")),
      group = stringr::str_c(oxygen, treatment, sep = "."),
      group = factor(group, levels = c("N.DMSO", "H.DMSO", "N.BAY", "H.BAY"))
    )

  design <- ~ experiment + group

  dds <-
    DESeq2::DESeqDataSetFromMatrix(
      countData = count_data,
      colData = pheno_data,
      design = design
    )

  keep <- rowSums(DESeq2::counts(dds)) > 1
  dds <- dds[keep, ]
  DESeq2::DESeq(dds)
}

# vst_rnaseq --------------------------------------------------------------

vst_rnaseq <- function(dds) {
  vsd <- DESeq2::vst(dds, blind = FALSE)
  SummarizedExperiment::assay(vsd) <-
    limma::removeBatchEffect(SummarizedExperiment::assay(vsd), vsd$experiment)

  DESeq2::plotPCA(
    vsd,
    intgroup = "group",
    returnData = TRUE
  ) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        group == "N.DMSO" ~ "21%\nDMSO",
        group == "H.DMSO" ~ "0.5%\nDMSO",
        group == "N.BAY" ~ "21%\nBAY",
        group == "H.BAY" ~ "0.5%\nBAY"
      )
    )

}

# plot_rnaseq_pca ---------------------------------------------------------

plot_rnaseq_pca <- function(pca_data) {

  percent_variance <- round(100 * attr(pca_data, "percentVar"))

  pair_clrs <- c(
    N.DMSO = "#b2df8a",
    H.DMSO = "#33a02c",
    N.BAY = "#cab2d6",
    H.BAY = "#6a3d9a"
  )

  ggplot2::ggplot(pca_data) +
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = group
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = group),
      pch = 21,
      color = "white",
      size = 2,
      show.legend = FALSE
    ) +
    ggforce::geom_mark_ellipse(
      ggplot2::aes(
        color = group,
        label = label
      ),
      expand = ggplot2::unit(2, "mm"),
      label.fontsize = 5,
      label.fontface = "plain",
      label.family = "Calibri",
      label.hjust = 0.5,
      label.buffer = ggplot2::unit(0, "mm"),
      label.margin = ggplot2::margin(-1.5, -1.5, -1.5, -1.5, "mm"),
      con.type = "none",
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = paste0("PC1: ", percent_variance[1], "% variance"),
      y = paste0("PC2: ", percent_variance[2], "% variance")
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = pair_clrs,
      labels = c("21% | DMSO", "0.5% | DMSO", "21% | BAY", "0.5% | BAY")
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = pair_clrs,
      labels = c("21% | DMSO", "0.5% | DMSO", "21% | BAY", "0.5% | BAY")
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.05)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.05)) +
    theme_plots() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(1, units = "lines")
    )

}

# identify_deg ------------------------------------------------------------

identify_deg <- function(dds, expr) {
  alpha <- 0.05
  fc <- 1

  mod_mat <- model.matrix(
    DESeq2::design(dds),
    SummarizedExperiment::colData(dds)
  )

  n.dmso <- colMeans(mod_mat[dds$group == "N.DMSO", ])
  n.bay <- colMeans(mod_mat[dds$group == "N.BAY", ])
  h.dmso <- colMeans(mod_mat[dds$group == "H.DMSO", ])
  h.bay <- colMeans(mod_mat[dds$group == "H.BAY", ])

  con <- eval(expr)
  # con <- (h.dmso - n.dmso) - (n.bay - n.dmso)

  DESeq2::results(
    dds,
    contrast = con,
    alpha = alpha,
    lfcThreshold = log(fc, base = 2),
    tidy = TRUE,
    parallel = TRUE
  ) %>%
    dplyr::rename(symbol = row) %>%
    dplyr::arrange(padj)

}

# plot_rnaseq_venn --------------------------------------------------------

plot_rnaseq_venn <- function(dds) {

  bay <-
    identify_deg(dds, expr(h.bay - n.bay)) %>%
    dplyr::filter(padj < 0.05)
  hyp <- identify_deg(dds, expr(h.dmso - n.dmso)) %>%
    dplyr::filter(padj < 0.05)

  nm <- rownames(SummarizedExperiment::assay(dds))
  bay_deg <- nm %in% bay$symbol
  hyp_deg <- nm %in% hyp$symbol

  tibble::tibble(nm, hyp_deg, bay_deg) %>%
    ggplot2::ggplot() +
    ggplot2::aes(A = hyp_deg, B = bay_deg) +
    ggvenn::geom_venn(
      set_names = c("0.5%", "BAY"),
      digits = 0,
      fill_color = clrs[c(2, 4)],
      fill_alpha = 0.5,
      stroke_size = 0.25,
      set_name_size = 6/ggplot2::.pt
    ) +
    ggplot2::coord_fixed() +
    theme_plots() +
    theme_void()

}

# plot_rnaseq_volcano -----------------------------------------------------

plot_rnaseq_volcano <- function(dds) {

  df <-
    identify_deg(dds, expr((h.dmso - n.dmso) - (n.bay - n.dmso))) %>%
    dplyr::arrange(desc(abs(stat))) %>%
    dplyr::mutate(fc = 2 ^ log2FoldChange)

  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = fc,
      y = padj
    ) +
    ggplot2::geom_hex(
      bins = 100,
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = head(df, 40),
      ggplot2::aes(label = symbol),
      size = 4/ggplot2::.pt,
      max.overlaps = 20
    ) +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::scale_y_continuous(trans = reverselog_trans(10)) +
    ggplot2::scale_x_continuous(
      trans = "log2",
      breaks = 2 ^ seq(-4, 4, 2),
      labels = scales::trans_format("log2", scales::math_format(2^.x))
    ) +
    ggplot2::labs(
      x = "ΔHypoxia/ΔBAY",
      y = "Adjusted P value"
    ) +
    theme_plots()

}

# plot_rnaseq_goi ---------------------------------------------------------

plot_rnaseq_goi <- function(dds, goi) {
  df <-
    tibble::as_tibble(
      SummarizedExperiment::assay(dds),
      rownames = "symbol"
    ) %>%
    tidyr::pivot_longer(
      -symbol,
      names_to = "sample",
      values_to = "count"
    ) %>%
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::colData(dds)),
      by = c("sample" = "id")
    )

  df %>%
    dplyr::filter(symbol %in% goi) %>%
    ggplot2::ggplot() +
    ggplot2::facet_wrap(~ symbol, scales = "free_y") +
    ggplot2::aes(
      x = treatment,
      y = count,
      color = oxygen
    ) +
    ggplot2::geom_point(
      alpha = 0.4,
      position =  ggplot2::position_dodge(width = 0.4),
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "crossbar",
      fun = "mean",
      width = 0.3,
      position =  ggplot2::position_dodge(width = 0.4),
      fatten = 2,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      x = "Treatment",
      y = "Count",
      color = "Oxygen"
    ) +
    ggplot2::ylim(c(0, NA)) +
    NULL
}
