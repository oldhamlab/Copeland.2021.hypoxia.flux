# rnaseq.R

count_rnaseq <- function() {

  dds <- rnaseq.lf.hypoxia.molidustat::lf_hyp_bay_rnaseq

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

  rownames(SummarizedExperiment::colData(dds)) <- SummarizedExperiment::colData(dds)$id

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
    ggplot2::geom_point(
      ggplot2::aes(fill = group),
      pch = 21,
      color = "white",
      size = 2,
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
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = c(3, 3))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(2, 2))) +
    theme_plots() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(1, units = "lines")
    )

}

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

  annots <-
    tibble::as_tibble(SummarizedExperiment::rowData(dds), rownames = "row") %>%
    dplyr::select(row, hgnc_symbol, description)

  DESeq2::results(
    dds,
    contrast = con,
    alpha = alpha,
    lfcThreshold = log(fc, base = 2),
    tidy = TRUE,
    parallel = TRUE
  ) %>%
    dplyr::left_join(annots, by = "row") %>%
    dplyr::rename(symbol = hgnc_symbol) %>%
    dplyr::relocate(symbol, description, .after = "row") %>%
    dplyr::arrange(padj)
}

find_rnaseq_overlap <- function(dds) {
  bay <-
    identify_deg(dds, expr(n.bay - n.dmso)) %>%
    dplyr::filter(padj < 0.05)
  hyp <- identify_deg(dds, expr(h.dmso - n.dmso)) %>%
    dplyr::filter(padj < 0.05)

  nm <- rownames(SummarizedExperiment::assay(dds))
  bay_deg <- nm %in% bay$row
  hyp_deg <- nm %in% hyp$row
  tibble::tibble(nm, hyp_deg, bay_deg)
}

find_gsea_overlap <- function(a, b, src = "HALLMARK") {
  dplyr::full_join(a, b, by = c("source", "pathway")) %>%
    dplyr::select(source, pathway, padj.x, padj.y) %>%
    dplyr::filter(source %in% src) %>%
    dplyr::mutate(
      hyp_deg = padj.x < 0.05,
      bay_deg = padj.y < 0.05,
      dplyr::across(tidyselect::matches("_deg"), tidyr::replace_na, FALSE)
    )
}

plot_rnaseq_venn <- function(df, title) {
  ggplot2::ggplot(df) +
    ggplot2::aes(A = hyp_deg, B = bay_deg) +
    ggvenn::geom_venn(
      set_names = c("0.5%", "BAY"),
      digits = 0,
      show_percentage = TRUE,
      fill_color = clrs[c(2, 4)],
      fill_alpha = 0.5,
      stroke_size = 0.25,
      set_name_size = 6/ggplot2::.pt,
      text_size = 5/ggplot2::.pt
    ) +
    ggplot2::labs(title = title) +
    theme_plots() +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

plot_rnaseq_volcano <- function(results, gois = NULL, xlab = NULL) {

  df <-
    tibble::as_tibble(results)

  left <-
    df %>%
    dplyr::filter(log2FoldChange < 0) %>%
    dplyr::slice_min(stat, n = 10)

  right <-
    df %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::slice_max(stat, n = 10)

  ggplot2::ggplot(df) +
    ggplot2::aes(
      # x = 2 ^ log2FoldChange,
      x = log2FoldChange,
      y = padj
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = symbol,
        color = symbol %in% gois
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = -8 - left$log2FoldChange,
      hjust = 0,
      direction = "y",
      family = "Calibri",
      segment.color = "black",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = symbol,
        color = symbol %in% gois
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = 8 - right$log2FoldChange,
      hjust = 1,
      direction = "y",
      segment.color = "black",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_hex(
      bins = c(250, 15),
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::scale_color_manual(values = c("black", "darkred")) +
    ggplot2::scale_y_continuous(trans = reverselog_trans(10)) +
    ggplot2::scale_x_continuous(
      breaks = seq(-6, 6, 2),
      labels = scales::math_format(2^.x),
      expand = ggplot2::expansion(mult = 0.3)
    ) +
    ggplot2::labs(
      x = xlab,
      y = "Adjusted P value"
    ) +
    theme_plots() +
    ggplot2::coord_cartesian(xlim = c(-6, 6), clip = "off")
}

get_unique_symbol_ids <- function(dds) {
  tibble::as_tibble(SummarizedExperiment::rowData(dds), rownames = "entrezid") %>%
    dplyr::filter(hgnc_symbol != "" & !is.na(hgnc_symbol)) %>%
    dplyr::group_by(hgnc_symbol) %>%
    dplyr::filter(baseMean == max(baseMean)) %>%
    dplyr::pull(entrezid)
}

dds_to_symbols <- function(dds, unique_symbol_ids) {
  df <- dds[rownames(dds) %in% unique_symbol_ids, ]
  rownames(df) <- SummarizedExperiment::rowData(df)$hgnc_symbol
  df
}

plot_rnaseq_goi <- function(dds, goi) {

  df <-
    dds[rownames(dds) %in% goi] %>%
    SummarizedExperiment::assay() %>%
    tibble::as_tibble(rownames = "symbol") %>%
    tidyr::pivot_longer(
      -symbol,
      names_to = "id",
      values_to = "count"
    ) %>%
    dplyr::left_join(
      dplyr::select(tibble::as_tibble(SummarizedExperiment::colData(dds)), id:group),
      by = "id",
      copy = TRUE
    ) %>%
    dplyr::mutate(oxygen = factor(oxygen, levels = c("N", "H"), labels = c("21%", "0.5%")))

  ggplot2::ggplot(df) +
    ggplot2::facet_wrap(~ symbol, scales = "free_y", nrow = 1) +
    ggplot2::aes(
      x = treatment,
      y = count/1000,
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
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::labs(
      x = "Treatment",
      y = expression(paste("Count (x", 10^3, ")")),
      fill = NULL
    ) +
    ggplot2::ylim(c(0, NA)) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}

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
    tidyr::separate(pathway, c("source", "pathway"), "_", extra = "merge")
}

plot_gsea <- function(rnaseq_gsea, sources, lbls, vals) {
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
      labels = lbls,
      values = vals
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

run_tfea <- function(df)
{
  df <- dplyr::rename(df, Genes = row)

  load("~/Dropbox (Partners HealthCare)/RM2020_GH_doubleElite.rdata")
  TFEA.ChIP::set_user_data(MetaData, Mat01)

  de_genes <-
    TFEA.ChIP::Select_genes(
      df,
      min_LFC = 0.5
    )

  ctl_genes <-
    TFEA.ChIP::Select_genes(
      df,
      min_pval = 0.5,
      max_pval = 1,
      min_LFC = -0.25,
      max_LFC = 0.25
    )

  TFEA.ChIP::contingency_matrix(de_genes, ctl_genes) %>%
    TFEA.ChIP::getCMstats() %>%
    TFEA.ChIP::rankTFs(rankMethod = "gsea")
}

plot_tfea <- function(rnaseq_tfea)
{
  df <-
    rnaseq_tfea %>%
    dplyr::mutate(padj = p.adjust(pVal, method = "fdr")) %>%
    dplyr::filter(padj < 0.05)

  ggplot2::ggplot(df) +
    ggplot2::aes(
      y = reorder(TF, ES),
      x = ES,
      fill = ES > 0
    ) +
    ggplot2::geom_col(
      show.legend = TRUE
    ) +
    ggplot2::labs(
      y = NULL,
      x = "Enrichment score"
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::scale_fill_manual(
      name = NULL,
      labels = c("With BAY", "With Hypoxia"),
      values = unname(clrs[c(4, 2)])
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = c("bottom"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}
