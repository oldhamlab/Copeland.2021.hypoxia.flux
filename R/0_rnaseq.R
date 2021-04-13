# rnaseq.R

# count_rnaseq ------------------------------------------------------------

count_rnaseq <- function() {

  dds <- rnaseq.lf.hypoxia.molidustat::lf_hyp_bay_se

  SummarizedExperiment::colData(dds) <-
    SummarizedExperiment::colData(dds) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      experiment = factor(experiment),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"), labels = c("N", "H")),
      treatment = factor(treatment, labels = c("DMSO", "BAY")),
      group = stringr::str_c(oxygen, treatment, sep = "."),
      group = factor(group, levels = c("N.DMSO", "H.DMSO", "N.BAY", "H.BAY"))
    ) %>%
    as("DataFrame")

  design <- ~ experiment + group

  dds <-
    DESeq2::DESeqDataSet(
      dds,
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
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::rowData(dds), rownames = "row"),
      by = "row"
    ) %>%
    dplyr::rename(symbol = hgnc_symbol) %>%
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

plot_rnaseq_volcano <- function(results) {

  df <-
    tibble::as_tibble(results)

  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = 2 ^ log2FoldChange,
      y = padj
    ) +
    ggplot2::geom_hex(
      bins = c(250, 15),
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
    dds[
      !is.na(SummarizedExperiment::rowData(dds)$hgnc_symbol) &
        SummarizedExperiment::rowData(dds)$hgnc_symbol %in% goi,
    ] %>%
    SummarizedExperiment::assay() %>%
    tibble::as_tibble(rownames = "entrez_id") %>%
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::rowData(dds), rownames = "entrez_id"),
      by = "entrez_id"
    ) %>%
    dplyr::select(symbol = hgnc_symbol, tidyselect::matches("^V\\d+")) %>%
    tidyr::pivot_longer(
      -symbol,
      names_to = "id",
      values_to = "count",
      names_prefix = "V",
      names_transform = list(id = ~sprintf("%02s", .x))
    ) %>%
    dplyr::left_join(SummarizedExperiment::colData(dds), by = "id", copy = TRUE)

  ggplot2::ggplot(df) +
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


# run_gsea ----------------------------------------------------------------

run_gsea <- function(results) {

  pathways <-
    fgsea::gmtPathways("~/Dropbox (Partners HealthCare)/msigdb_v7.2/msigdb_v7.2_GMTs/msigdb.v7.2.entrez.gmt")

  rnks <-
    results %>%
    dplyr::select(row, stat) %>%
    dplyr::arrange(stat) %>%
    tibble::deframe()

  fgsea::fgsea(
    pathways = pathways,
    stats = rnks,
    nPermSimple = 10000,
    eps = 0,
    BPPARAM = BiocParallel::bpparam()
  ) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(desc(NES)) %>%
    tidyr::separate(pathway, c("source", "pathway"), "_", extra = "merge") %>%
    dplyr::filter(padj < 0.05)

}


# plot_gsea ---------------------------------------------------------------

plot_gsea <- function(rnaseq_gsea, sources) {
  rnaseq_gsea %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::filter(source %in% sources) %>%
    dplyr::select(source, pathway, NES) %>%
    dplyr::mutate(
      pathway = stringr::str_replace_all(pathway, "_", " "),
      pathway = tolower(pathway)
    ) %>%
    dplyr::distinct() %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = NES,
      y = reorder(pathway, NES)
    ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = NES > 0)
    ) +
    ggplot2::labs(
      x = "Normalized Enrichment Score",
      y = NULL
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::scale_fill_manual(
      name = NULL,
      labels = c("Up in BAY", "Up in Hypoxia"),
      values = unname(clrs[c(4, 2)])
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}
