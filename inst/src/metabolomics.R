# metabolomics.R

# format_metab_targeted ---------------------------------------------------

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
    dplyr::filter(!is.na(hmdb)) %>%
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
      s_data = c("type", "oxygen", "treatment", "replicate")
    ) %>%
    tbl_to_se()

}

# new_tbl_se --------------------------------------------------------------

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

# tibble_to_se ------------------------------------------------------------

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
    tibble::column_to_rownames(attr(tbl_se, "f_names"))

  sample_data <-
    tbl_se %>%
    dplyr::select(
      attr(tbl_se, "s_names"),
      attr(tbl_se, "s_data")
    ) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames(attr(tbl_se, "s_names"))

    SummarizedExperiment::SummarizedExperiment(
    assays = assay_data,
    rowData = feature_data,
    colData = sample_data
  )

}

# remove_missing_metab ----------------------------------------------------

remove_missing_metab <- function(raw) {

  qc <- SummarizedExperiment::assay(raw[, raw$type == "qc"])

  missing <- names(which(apply(qc, 1, function(x) sum(is.na(x))) > 0))

  raw[rownames(raw) %nin% missing, raw$type %nin% c("water", "blank")]

}

# correct_drift -----------------------------------------------------------

correct_drift <- function(missing) {

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

# quality_control ---------------------------------------------------------

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

# impute_missing ----------------------------------------------------------

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

# normalize ---------------------------------------------------------------

normalize <- function(imputed) {
  mat <- SummarizedExperiment::assay(imputed)
  quotients <- mat / SummarizedExperiment::rowData(imputed)$reference
  quotient_medians <- apply(quotients, 2, median)
  SummarizedExperiment::assay(imputed) <- t(t(mat) / quotient_medians)
  imputed
}


