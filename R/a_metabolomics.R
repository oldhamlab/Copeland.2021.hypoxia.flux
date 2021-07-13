# metabolomics.R

format_metab_targeted <- function(metab_targeted_files) {
  readxl::read_excel(
    metab_targeted_files,
    sheet = 2
  ) %>%
    dplyr::select(
      id = `Raw File Name`,
      sample = `Sample ID`,
      metabolite = `Compound Name`,
      mz = `Detected Mass`,
      rt = RT,
      area = `Peak Area`
    ) %>%
    dplyr::left_join(wmo::hmdb_mappings, by = "metabolite") %>%
    dplyr::mutate(metabolite = dplyr::case_when(
      metabolite == "glyceraldehyde 3-phosphate" ~ "GAP",
      TRUE ~ metabolite
    )) %>%
    dplyr::filter(!is.na(hmdb)) %>%
    dplyr::filter(!stringr::str_detect(hmdb, "ISTD")) %>%
    dplyr::mutate(
      id = stringr::str_c("S", sprintf("%02d", id)),
      type = dplyr::case_when(
        sample == "water" ~ "blank",
        stringr::str_detect(sample, "mix") ~ "qc",
        TRUE ~ "sample"
      ),
      oxygen = stringr::str_extract(sample, "hyp|norm"),
      oxygen = factor(oxygen, levels = c("norm", "hyp"), labels = c("N", "H")),
      treatment = stringr::str_extract(sample, "bay|dmso"),
      treatment = factor(treatment, levels = c("dmso", "bay"), labels = c("DMSO", "BAY")),
      group = stringr::str_c(oxygen, treatment, sep = "."),
      group = factor(group, levels = c("N.DMSO", "N.BAY", "H.DMSO", "H.BAY")),
      replicate = stringr::str_extract(sample, "(?<=\\.)\\w{1}$"),
      mz = replace(mz, mz == "N/F", NA_real_),
      mz = as.numeric(mz),
      dplyr::across(c(rt, area), ~replace(., . == 0, NA_real_))
    ) %>%
    dplyr::select(-sample) %>%
    new_tbl_se(
      a_data = "area",
      f_names = "hmdb",
      f_data = c("metabolite", "mz", "rt"),
      s_names = "id",
      s_data = c("type", "oxygen", "treatment", "replicate", "group")
    ) %>%
    tbl_to_se()

}

new_tbl_se <- function(
  tbl,
  a_data,
  f_names,
  f_data = NULL,
  s_names,
  s_data = NULL
) {
  structure(
    tbl,
    class = c("tbl_se", class(tbl)),
    a_data = a_data,
    f_names = f_names,
    s_names = s_names,
    f_data = f_data,
    s_data = s_data
  )
}

tbl_to_se <- function(tbl_se, assay_name) {

  assay_data <-
    tbl_se %>%
    dplyr::select(
      attr(tbl_se, "f_names"),
      attr(tbl_se, "s_names"),
      attr(tbl_se, "a_data")
    ) %>%
    tidyr::pivot_wider(
      names_from = attr(tbl_se, "s_names"),
      values_from = attr(tbl_se, "a_data")
    ) %>%
    tibble::column_to_rownames(attr(tbl_se, "f_names"))

  feature_data <-
    tbl_se %>%
    dplyr::select(
      attr(tbl_se, "f_names"),
      attr(tbl_se, "f_data")
    ) %>%
    dplyr::group_by(!!rlang::sym(attr(tbl_se, "f_names"))) %>%
    dplyr::summarise(
      metabolite = unique(metabolite),
      mz = mean(mz, na.rm = TRUE),
      mz_min = min(mz, na.rm = TRUE),
      mz_max = max(mz, na.rm = TRUE),
      rt = mean(rt, na.rm = TRUE),
      rt_min = min(rt, na.rm = TRUE),
      rt_max = max(rt, na.rm = TRUE)
    ) %>%
    tibble::column_to_rownames(attr(tbl_se, "f_names")) %>%
    .[match(rownames(assay_data), rownames(.)), ]

  sample_data <-
    tbl_se %>%
    dplyr::select(
      attr(tbl_se, "s_names"),
      attr(tbl_se, "s_data")
    ) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames(attr(tbl_se, "s_names")) %>%
    .[match(colnames(assay_data), rownames(.)), ]

  SummarizedExperiment::SummarizedExperiment(
    assays = assay_data,
    rowData = feature_data,
    colData = sample_data
  )

}

remove_missing_metab <- function(raw) {

  qc <- SummarizedExperiment::assay(raw[, raw$type == "qc"])

  missing <- names(which(apply(qc, 1, function(x) sum(is.na(x))) > 0))

  raw[rownames(raw) %nin% missing, raw$type %nin% c("water", "blank")]

}

prepare_assay_data <- function(se) {
  SummarizedExperiment::assay(se) %>%
    tibble::rownames_to_column("hmdb") %>%
    tidyr::pivot_longer(-hmdb, names_to = "sample", values_to = "value") %>%
    dplyr::mutate(
      value = log(value),
      run_order = as.numeric(stringr::str_extract(sample, "\\d{2}"))
    ) %>%
    dplyr::group_by(hmdb) %>%
    tidyr::nest()
}

correct_drift <- function(missing) {

  models <-
    missing[, missing$type == "qc"] %>%
    prepare_assay_data() %>%
    dplyr::mutate(
      model = map(data, ~smooth.spline(x = .x$run_order, y = .x$value, spar = 0.2)),
      mean = map_dbl(data, ~mean(.x$value))
    ) %>%
    dplyr::select(-data)

  corrected <-
    missing %>%
    prepare_assay_data() %>%
    dplyr::left_join(models, by = "hmdb") %>%
    dplyr::mutate(pred = purrr::map2(model, data, ~predict(.x, .y$run_order)$y)) %>%
    tidyr::unnest(c(data, pred)) %>%
    dplyr::mutate(corr = value + mean - pred) %>%
    dplyr::select(hmdb, sample, corr) %>%
    tidyr::pivot_wider(names_from = sample, values_from = corr) %>%
    tibble::column_to_rownames("hmdb") %>%
    exp()

  SummarizedExperiment::assay(missing) <- corrected
  missing

}

quality_control <- function(drift) {
  rsd <-
    drift[, drift$type == "qc"] %>%
    SummarizedExperiment::assay() %>%
    apply(1, function(x) 1.4826 * mad(x, na.rm = TRUE) / median(x, na.rm = TRUE))

  mad_qc <-
    drift[, drift$type == "qc"] %>%
    SummarizedExperiment::assay() %>%
    apply(1, function(x) mad(x, na.rm = TRUE))

  ref <-
    drift[, drift$type == "qc"] %>%
    SummarizedExperiment::assay() %>%
    apply(1, function(x) median(x, na.rm = TRUE))

  mad_s <-
    drift[, drift$type == "sample"] %>%
    SummarizedExperiment::assay() %>%
    apply(1, function(x) mad(x, na.rm = TRUE))

  d_ratio <- mad_qc / mad_s

  SummarizedExperiment::rowData(drift)$rsd <- rsd
  SummarizedExperiment::rowData(drift)$rsd <- d_ratio
  SummarizedExperiment::rowData(drift)$good <- rsd < 0.2 & d_ratio < 0.4
  SummarizedExperiment::rowData(drift)$reference <- ref

  drift[SummarizedExperiment::rowData(drift)$good == TRUE, drift$type != "qc"]
}

impute_missing <- function(qc) {
  set.seed(42)
  SummarizedExperiment::assay(qc) <-
    missForest::missForest(
      t(SummarizedExperiment::assay(qc)),
      maxiter = 10
    )$ximp %>%
    t()
  qc
}

pqn <- function(imputed) {
  mat <- SummarizedExperiment::assay(imputed)
  quotients <- mat / SummarizedExperiment::rowData(imputed)$reference
  quotient_medians <- apply(quotients, 2, median)
  SummarizedExperiment::assay(imputed) <- t(t(mat) / quotient_medians)
  imputed
}

plot_metab_pca <- function(clean) {

  input <-
    SummarizedExperiment::assay(clean) %>%
    log() %>%
    t() %>%
    scale() %>%
    t()

  pheno <- SummarizedExperiment::colData(clean)

  design <- model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- stringr::str_extract(colnames(design), "(?<=group).*")

  df <-
    limma::removeBatchEffect(input, batch = pheno$replicate, design = design) %>%
    t() %>%
    pcaMethods::pca(scale = "none", center = TRUE)

  percent_variance <- round(100 * c(df@R2[[1]], df@R2[[2]]))

  pair_clrs <- c(
    N.DMSO = "#b2df8a",
    H.DMSO = "#33a02c",
    N.BAY = "#cab2d6",
    H.BAY = "#6a3d9a"
  )

  merge(pcaMethods::scores(df), SummarizedExperiment::colData(clean), by = 0) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        group == "N.DMSO" ~ "21%\nDMSO",
        group == "H.DMSO" ~ "0.5%\nDMSO",
        group == "N.BAY" ~ "21%\nBAY",
        group == "H.BAY" ~ "0.5%\nBAY"
      )
    ) %>%
    ggplot2::ggplot() +
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

fit_metab_limma <- function(clean) {

  input <-
    SummarizedExperiment::assay(clean) %>%
    log() %>%
    t() %>%
    scale() %>%
    t()

  pheno <- SummarizedExperiment::colData(clean)

  design <- model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- stringr::str_extract(colnames(design), "(?<=group).*")

  corfit <- limma::duplicateCorrelation(input, design, block = pheno$replicate)

  cm <-
    limma::makeContrasts(
      hyp_on_dmso = H.DMSO - N.DMSO,
      bay_on_norm = N.BAY - N.DMSO,
      hyp_on_bay = H.BAY - N.BAY,
      deltas = (H.DMSO - N.DMSO) - (N.BAY - N.DMSO),
      levels = design
    )

  fit <-
    limma::lmFit(
      input,
      design,
      block = pheno$replicate,
      correlation = corfit$consensus
    ) %>%
    limma::contrasts.fit(cm) %>%
    limma::eBayes()
}

metab_top_table <- function(clean, fit, contrast) {
  limma::topTable(
    fit,
    number = Inf,
    p.value = 1,
    coef = contrast
  ) %>%
    tibble::as_tibble(rownames = "hmdb") %>%
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::rowData(clean), rownames = "hmdb"),
      by = "hmdb"
    )
}

plot_metab_volcano <- function(results) {

  left <-
    results %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::slice_min(t, n = 10)

  right <-
    results %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::slice_max(t, n = 10)

  ggplot2::ggplot(results) +
    ggplot2::aes(
      # x = 2 ^ log2FoldChange,
      x = logFC,
      y = adj.P.Val
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = metabolite,
        color = metabolite %in% c("2-hydroxyglutarate","aconitate", "glyceraldehyde 3-phosphate", "hydroxyproline", "malate", "taurine")
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = -4,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = metabolite,
        color = metabolite %in% c("2-hydroxyglutarate","aconitate", "glyceraldehyde 3-phosphate", "hydroxyproline", "malate", "taurine")
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      segment.color = "black",
      nudge_x = 4.5,
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = subset(results, adj.P.Val > 0.05),
      pch = 21,
      color = "white",
      fill = "grey80"
    ) +
    ggplot2::geom_point(
      data = subset(results, logFC > 0 & adj.P.Val < 0.05),
      pch = 21,
      color = "white",
      fill = clrs[[2]]
    ) +
    ggplot2::geom_point(
      data = subset(results, logFC < 0 & adj.P.Val < 0.05),
      pch = 21,
      color = "white",
      fill = clrs[[4]]
    ) +
    ggplot2::scale_color_manual(values = c("black", "darkred")) +
    ggplot2::scale_y_continuous(trans = reverselog_trans(10)) +
    ggplot2::scale_x_continuous(
      breaks = seq(-4, 4, 2),
      labels = scales::math_format(2^.x),
      expand = ggplot2::expansion(mult = 1)
    ) +
    ggplot2::labs(
      x = "ΔHypoxia/ΔBAY",
      y = "Adjusted P value"
    ) +
    theme_plots()
}

se_to_tibble <- function(se) {
  tibble::as_tibble(SummarizedExperiment::assay(se), rownames = "feature") %>%
    tidyr::pivot_longer(-feature, names_to = "sample", values_to = "area") %>%
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::colData(se), rownames = "sample"),
      by = "sample"
    ) %>%
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::rowData(se), rownames = "feature"),
      by = "feature"
    )
}

plot_mois <- function(clean, moi) {

  se_to_tibble(clean) %>%
    dplyr::filter(metabolite %in% moi) %>%
    dplyr::group_by(metabolite) %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("N", "H"), labels = c("21%", "0.5%")),
      area = area / mean(area[oxygen == "21%" & treatment == "DMSO"]),
      metabolite = dplyr::case_when(
        metabolite == "2-hydroxyglutarate" ~ "2HG",
        metabolite == "glyceraldehyde 3-phosphate" ~ "GAP",
        TRUE ~ metabolite
      )
    ) %>%
    ggplot2::ggplot() +
    ggplot2::facet_wrap(~ metabolite, nrow = 1) +
    ggplot2::aes(
      x = treatment,
      y = area,
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
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::labs(
      x = "Treatment",
      y = "Peak area (normalized)",
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

run_msea <- function(tt, pathways) {

  stats <- tt$t
  names(stats) <- tt$hmdb

  msea_res <-
    fgsea::fgsea(pathways = pathways, stats = stats) %>%
    dplyr::filter(pval < 0.05) %>%
    tidyr::separate(pathway, c("source", "pathway"), "\\) ") %>%
    dplyr::mutate(source = stringr::str_replace(source, "\\(", "")) %>%
    dplyr::arrange(-abs(NES))
}

get_metab_pathways <- function(databases = "kegg") {
  multiGSEA::getMultiOmicsFeatures(
    dbs = databases,
    layer = "metabolome",
    organism = "hsapiens",
    returnMetabolome = "HMDB",
    useLocal = TRUE
  ) %>%
    unlist(recursive = FALSE) %>%
    rlang::set_names(purrr::map(names(.), stringr::str_extract, pattern = "(?<=metabolome\\.).*"))
}

plot_msea <- function(msea) {
  msea %>%
    dplyr::select(-padj) %>%
    dplyr::rename(padj = pval) %>%
    plot_gsea(sources = "KEGG") +
    ggplot2::theme(
      axis.text.y.right = element_text(size = 5)
    )
}

plot_leading_edge <- function(tt, pathway) {

  stats <- tt$t
  names(stats) <- tt$hmdb

  rnk <- rank(-stats)
  ord <- order(rnk)

  stats_adj <- stats[ord]
  stats_adj <- stats_adj / max(abs(stats_adj))

  pathway <- unname(as.vector(na.omit(match(pathway, names(stats_adj)))))
  pathway <- sort(pathway)

  gsea_res <-
    fgsea::calcGseaStat(
      stats_adj,
      selectedStats = pathway,
      returnAllExtremes = TRUE
    )

  bottoms <- gsea_res$bottoms
  tops <- gsea_res$tops

  n <- length(stats_adj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  nms <- as.vector(rbind(names(bottoms), names(tops)))
  df <- tibble::tibble(
    x = c(0, xs, n + 1),
    y = c(0, ys, 0),
    names = c(NA_character_, nms, NA_character_)
  ) %>%
    dplyr::left_join(wmo::hmdb_mappings, by = c("names" = "hmdb"))

  annot <-
    df %>%
    dplyr::filter(!is.na(metabolite)) %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(x = max(x)) %>%
    dplyr::mutate(metabolite = dplyr::case_when(
      metabolite == "2-oxoglutarate" ~ "AKG",
      metabolite == "phosphoenolpyruvate" ~ "PEP",
      metabolite == "malate" ~ "MAL",
      metabolite == "fumarate" ~ "FUM",
      metabolite == "aconitate" ~ "ACO",
      metabolite == "pyruvate" ~ "PYR",
      metabolite == "citrate" ~ "CIT",
      metabolite == "succinate" ~ "SUC",
      TRUE ~ metabolite
    ))

  diff <- (max(tops) - min(bottoms)) / 8

  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = x,
      y = y
    ) +
    ggplot2::geom_line(color = clrs[[3]]) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "black",
      size = 0.25
    ) +
    ggplot2::geom_segment(
      data = data.frame(x = pathway),
      ggplot2::aes(
        x = x,
        y = -diff/2,
        xend = x,
        yend = diff/2
      ),
      size = 0.1) +
    ggrepel::geom_text_repel(
      data = annot,
      ggplot2::aes(
        y = diff/2,
        label = metabolite
      ),
      angle = 90,
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_y = 0.2,
      # nudge_x = -3,
      hjust = 1,
      # vjust = 0.5,
      direction = "x",
      min.segment.length = 0.3
    ) +
    ggplot2::labs(
      x = "Rank",
      y = "Enrichment score",
      title = "KEGG: Citrate cycle"
    ) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
    theme_plots() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        margin = ggplot2::margin(b = 1),
        size = 8
      )
    ) +
    NULL
}

format_sample_info <- function(filename) {
  readr::read_csv(filename, col_types = readr::cols()) %>%
    dplyr::select(
      file = Filename,
      type = `Sample type`,
      id = `Sample ID`
    ) %>%
    dplyr::mutate(
      file = stringr::str_c("S", file, ".mzML"),
      run_order = stringr::str_extract(file, "\\d+") %>% as.numeric(),
      type = dplyr::case_when(
        id == "water" ~ "blank",
        stringr::str_detect(id, "mix") ~ "qc",
        TRUE ~ "sample"
      ),
      oxygen = dplyr::case_when(
        stringr::str_detect(id, "hyp") ~ "Hyp",
        stringr::str_detect(id, "norm") ~ "Norm",
        TRUE ~ NA_character_
      ),
      treatment = dplyr::case_when(
        stringr::str_detect(id, "dmso") ~ "DMSO",
        stringr::str_detect(id, "bay") ~ "BAY",
        TRUE ~ NA_character_
      ),
      replicate = stringr::str_extract(id, "(?<=(dmso|bay)\\.)\\w{1}$"),
      replicate = dplyr::case_when(
        replicate == "a" ~ "11",
        replicate == "b" ~ "12",
        TRUE ~ stringr::str_c(2, replicate)
      ),
      group = dplyr::case_when(
        type == "qc" ~ "QC",
        type == "sample" ~ stringr::str_c(oxygen, treatment, sep = "."),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(file, run_order, type, oxygen, treatment, replicate, group)
}

read_metab_raw <- function(data_files, polarity, samples) {
  data_files[which(stringr::str_detect(data_files, polarity))] %>%
    MSnbase::readMSData(
      files = .,
      pdata = new("NAnnotatedDataFrame", samples),
      mode = "onDisk"
    ) %>%
    MSnbase::filterRt(c(160, 1200))
}

get_tics <- function(msnexp, samples) {
  tics <- MSnbase::chromatogram(msnexp, aggregationFun = "sum")
  rt <- lapply(unlist(tics), MSnbase::rtime)
  int <- lapply(unlist(tics), MSnbase::intensity)
  purrr::map2_dfr(rt, int, ~dplyr::bind_cols(rt = .x, intensity = .y), .id = "run_order") %>%
    dplyr::mutate(run_order = as.numeric(run_order)) %>%
    dplyr::left_join(samples, by = "run_order")
}

plot_chromatograms <- function(df) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = rt,
      y = intensity,
      color = group,
      group = run_order
    ) +
    ggplot2::geom_line(size = 0.25) +
    ggplot2::labs(color = NULL)
}

plot_intensities <- function(df) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = run_order,
      y = intensity,
      group = run_order,
      fill = group
    ) +
    ggplot2::geom_boxplot(
      outlier.alpha = 0.02
    )
}

optimize_centwave_params <- function(raw, polarity) {
  parameter_list <-
    list(
      ppm = 2.5,
      min_peakwidth = c(20, 40),
      max_peakwidth = c(300, 500),
      snthresh = 1,
      prefilter_k = 3,
      prefilter_int = 4500,
      mzdiff = c(-0.01, 0.1),
      noise = 1000
    )

  out_dir <-
    stringr::str_c(
      system.file(
        "data-raw/metabolomics/lf_05-bay",
        package = "Copeland.2021.hypoxia.flux"
      ),
      "/optimize/",
      polarity
    )

  test <- MSnbase::filterFile(raw, "S10.mzML")

  IPO2::optimize_centwave(test, parameter_list, out_dir)
}

plot_metabolite <- function(xcmsnexp, samples, mz, rt) {
  mzr <- mz + c(-0.02, 0.02)
  rtr <- 60 * (rt + c(-2, 2))
  object <- filterFile(xcmsnexp, samples)
  chr_ex <- chromatogram(object, mz = mzr, rt = rtr)
  plot(chr_ex, peakType = "rectangle", peakBg = NA)
}

optimize_align_group_params <- function(merged, polarity) {
  out_dir <-
    stringr::str_c(
      system.file(
        "data-raw/metabolomics/lf_05-bay",
        package = "Copeland.2021.hypoxia.flux"
      ),
      "/optimize/",
      polarity
    )

  test <- MSnbase::filterFile(merged, seq(5, 40, 7))

  IPO2::optimize_align_group(test, out_dir = out_dir)
}

adjust_rtime <- function(object, param, center_sample = integer()) {
  centerSample(param) <- center_sample
  adjustRtime(object, param)
}

group_peaks <- function(object, param) {
  sampleGroups(param) <- object$group
  groupChromPeaks(object, param = param)
}
