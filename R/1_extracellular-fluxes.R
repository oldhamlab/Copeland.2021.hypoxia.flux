# extracellular-fluxes.R

clean_dna_per_cell <- function(filename) {
  filename %>%
    read_multi_excel() %>%
    purrr::map(clean_technical_replicates) %>%
    dplyr::bind_rows(.id = "id") %>%
    tidyr::separate(id, c("cell_type", "volume"), sep = "_", convert = TRUE) %>%
    dplyr::mutate(cells = 1000 * cells)
}

calculate_cells_per_dna <- function(tbl) {
  tbl %>%
    dplyr::filter(!(cell_type == "lf" & volume == "100" & cells == 300000)) %>%
    dplyr::filter(!(cell_type == "pasmc" & volume == "200" & cells == 400000)) %>%
    dplyr::group_by(cell_type, volume) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = map(data, ~lm(conc ~ 0 + cells, data = .x, na.action = modelr::na.warn)),
      glance = map(model, broom::tidy)
    ) %>%
    tidyr::unnest(c(glance)) %>%
    dplyr::select(cell_type, volume, slope = estimate) %>%
    dplyr::mutate(slope = 1/slope)
}

clean_flux_meta <- function(file_list) {
  file_list %>%
    rlang::set_names(stringr::str_extract(basename(.), "^.+(?=_)")) %>%
    purrr::map_dfr(
      readr::read_csv,
      col_types = "dcccdc",
      .id = "experiment"
    ) %>%
    tidyr::separate(experiment, c("cell_type", "experiment"), "_")
}

assemble_flux_data <- function(file_list) {
  nms <- sub("\\..*$", "", basename(file_list))
  sheets <- unique(unlist(purrr::map(file_list, readxl::excel_sheets)))
  data_list <- purrr::map(file_list, read_multi_excel)

  purrr::map(sheets, ~purrr::map(data_list, .x)) %>%
    rlang::set_names(sheets) %>%
    purrr::map(rlang::set_names, nms) %>%
    purrr::map(dplyr::bind_rows, .id = "experiment")
}

clean_fluxes <- function(data_list) {
  df <-
    data_list[c("dna", "glc", "lac", "pyr")] %>%
    dplyr::bind_rows(.id = "metabolite") %>%
    dplyr::filter(!is.na(.data$a)) %>%
    dplyr::select(where(~any(!is.na(.)))) %>%
    clean_technical_replicates() %>%
    tidyr::separate(.data$experiment, c("cell_type", "experiment", "batch", "date"), "_")

  dna <-
    df %>%
    dplyr::filter(.data$metabolite == "dna") %>%
    dplyr::mutate(volume = dplyr::if_else(.data$date <= "2018-05-25", 100, 200)) %>%
    dplyr::left_join(cells_per_dna, by = c("cell_type", "volume")) %>%
    dplyr::mutate(
      conc = .data$conc * slope,
      metabolite = "dna",
      detector = "picogreen"
    ) %>%
    dplyr::select(-c(.data$volume, .data$slope))

  others <-
    df %>%
    dplyr::filter(.data$metabolite != "dna") %>%
    dplyr::mutate(
      conc = dplyr::case_when(
        metabolite == "lac" & batch == "a" ~ .data$conc * 10.5,
        metabolite == "lac" & batch != "a" ~ .data$conc * 10,
        metabolite == "glc" ~ .data$conc * 555.074,
        metabolite == "pyr" ~ .data$conc * 20
      ),
      metabolite = dplyr::case_when(
        metabolite == "lac" ~ "lactate",
        metabolite == "glc" ~ "glucose",
        metabolite == "pyr" ~ "pyruvate"
      ),
      detector = "enzyme"
    )

  pyr <-
    data_list[["pyr"]] %>%
    dplyr::filter(!is.na(.data$pyruvate)) %>%
    tidyr::separate(.data$experiment, c("cell_type", "experiment", "batch", "date"), "_") %>%
    dplyr::mutate(istd = dplyr::coalesce(.data$KV, .data$`d8-valine`)) %>%
    dplyr::mutate(
      detector = case_when(
        !is.na(.data$KV) ~ "hplc",
        !is.na(.data$`d8-valine`) ~ "lcms",
        TRUE ~ "enzyme"
      ),
      istd = dplyr::case_when(
        experiment == "05" & batch == "a" & run == "a" & !is.na(.data$conc) ~ istd * 25,
        TRUE ~ .data$istd
      ),
      value = pyruvate / istd,
      metabolite = "pyruvate"
    ) %>%
    dplyr::select(.data$metabolite, .data$cell_type:.data$conc, .data$value, .data$detector)

  istds <- c("Norvaline", "Sarcosine")
  secondary_aa <- c("Hydroxyproline", "Proline")

  aa <-
    data_list[["aa"]] %>%
    dplyr::mutate(
      dplyr::across(
        c(tidyselect::contains("1") & !tidyselect::contains(c(istds, secondary_aa))),
        ~ . /.data$`1 Norvaline`
      ),
      dplyr::across(
        c(tidyselect::contains("2") & !tidyselect::contains(c(istds, secondary_aa))),
        ~ . /.data$`2 Norvaline`
      ),
      dplyr::across(
        c(tidyselect::contains("1") & tidyselect::contains(secondary_aa)),
        ~ . /.data$`1 Sarcosine`
      ),
      dplyr::across(
        c(tidyselect::contains("2") & tidyselect::contains(secondary_aa)),
        ~ . /.data$`2 Sarcosine`
      ),
    ) %>%
    tidyr::pivot_longer(matches("\\d .*"), names_to = "metabolite", values_to = "value") %>%
    tidyr::separate(metabolite, c("detector", "metabolite"), " ") %>%
    tidyr::separate(.data$experiment, c("cell_type", "experiment", "batch", "date"), "_") %>%
    dplyr::mutate(
      metabolite = tolower(.data$metabolite),
      detector = dplyr::if_else(.data$detector == 1, "mwd", "fld"),
      conc = replace(
        .data$conc,
        .data$experiment == "bay" &
          .data$date %in% c("2018-11-06", "2018-11-11") &
          .data$metabolite %in% c("asparagine", "glutamine", "tryptophan") &
          .data$conc == 225,
        22.5),
      conc = dplyr::case_when(
        .data$batch == "a" ~ 200 / 180 * .data$conc,
        TRUE ~ 200 / 190 * .data$conc
      )
    ) %>%
    dplyr::filter(!(.data$cell_type == "pasmc" & .data$metabolite == "glutamine"))

  gln <-
    data_list[["gln"]] %>%
    dplyr::mutate(
      value = .data$`1 Glutamine` / .data$`1 Norvaline`,
      detector = "mwd",
      conc = 20 * .data$conc,
      metabolite = "glutamine"
    ) %>%
    dplyr::select(.data$experiment:.data$conc, .data$detector, .data$value, .data$metabolite) %>%
    tidyr::separate(.data$experiment, c("cell_type", "experiment", "batch", "date"), "_")

  dplyr::bind_rows(dna, others, pyr, aa, gln) %>%
    dplyr::relocate(.data$detector, .after = .data$metabolite) %>%
    dplyr::arrange(metabolite, detector, cell_type, experiment, batch, date, run, conc, id) %>%
    dplyr::mutate(detector = tidyr::replace_na(.data$detector, "na")) %>%
    dplyr::filter(!is.na(.data$value)) %>%
    dplyr::filter(metabolite %nin% c(tolower(istds), "hydroxyproline"))
}

clean_flux_std <- function(df) {
  outliers <-
    tibble::tribble(
      ~metabolite, ~experiment, ~date, ~batch, ~run, ~detector, ~conc,
      "lactate", "bay", "2018-06-02", "a", "a", "enzyme", 10500,
      "pyruvate", "bay", "2018-11-02", "b", "a", "enzyme", 1500,
      "glycine", "05", "2017-11-06", "a", "b", "fld", 10
    )

  df %>%
    dplyr::filter(!(detector == "fld" & conc > 900)) %>%
    dplyr::anti_join(outliers, by = c(
      "metabolite",
      "experiment",
      "date",
      "batch",
      "run",
      "detector",
      "conc"
    )) %>%
    make_std_curves()
}

fill_missing_fluxes <- function(df, meta) {
  df_meta <-
    dplyr::left_join(df, meta, by = c("cell_type", "experiment", "id"))

  missing_data <-
    df_meta %>%
    dplyr::group_by(cell_type, experiment, batch) %>%
    tidyr::complete(
      .data$date,
      .data$treatment,
      .data$oxygen,
      tidyr::nesting(
        metabolite,
        type,
        detector,
        time,
        well
      )
    ) %>%
    dplyr::filter(is.na(.data$conc))

  empty_05_t0 <-
    missing_data %>%
    dplyr::filter(experiment == "05" & oxygen == "0.5%" & time == 0 & type == "empty") %>%
    dplyr::select(-c(.data$run, .data$id, .data$conc)) %>%
    dplyr::mutate(oxygen = forcats::fct_recode(oxygen, "21%" = "0.5%")) %>%
    dplyr::left_join(
      df_meta,
      by = c(
        "cell_type",
        "experiment",
        "batch",
        "date",
        "treatment",
        "oxygen",
        "metabolite",
        "type",
        "detector",
        "time",
        "well"
      )
    ) %>%
    dplyr::mutate(oxygen = forcats::fct_recode(oxygen, "0.5%" = "21%"))

  dplyr::bind_rows(df_meta, empty_05_t0)

}

filter_assays <- function(df) {
  mwd <- c("cystine", "glutamine", "isoleucine", "leucine", "lysine", "valine")

  df %>%
    dplyr::mutate(keep = dplyr::case_when(
      metabolite %in% mwd & detector == "mwd" ~ TRUE,
      metabolite %nin% mwd & detector == "fld" ~ TRUE,
      detector %in% c("enzyme", "picogreen") ~ TRUE,
      metabolite == "pyruvate" & experiment == "02" ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    dplyr::filter(.data$keep) %>%
    dplyr::select(-.data$keep) %>%
    dplyr::ungroup()
}

assemble_evap_data <- function(data_list) {
  data_list[["evap"]] %>%
    dplyr::group_by(.data$experiment, .data$oxygen, .data$treatment) %>%
    dplyr::mutate(plate_mass = .data$mass - .data$mass[[1]]) %>%
    dplyr::filter(.data$time != -24) %>%
    dplyr::mutate(volume = 2 * .data$plate_mass / .data$plate_mass[[1]]) %>%
    dplyr::filter(!is.na(.data$volume)) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(.data$data, ~lm(volume ~ time, data = .x))
    ) %>%
    dplyr::mutate(pred_vol = purrr::map2(.data$model, .data$data, predict)) %>%
    dplyr::select(-.data$model) %>%
    tidyr::unnest(c(.data$data, .data$pred_vol)) %>%
    dplyr::select(-c(.data$volume, .data$plate_mass, .data$mass)) %>%
    dplyr::rename(volume = .data$pred_vol) %>%
    tidyr::separate(.data$experiment, c("cell_type", "experiment", "batch", "date"), "_")
}

fill_missing_evap <- function(evap, samples) {
  evap <-
    evap %>%
    dplyr::filter(.data$treatment == "DMSO") %>%
    dplyr::mutate(treatment = "BAY") %>%
    dplyr::bind_rows(evap)

  evap_dup_bay <-
    evap %>%
    dplyr::filter(.data$experiment == "bay") %>%
    dplyr::group_by(.data$oxygen, .data$treatment, .data$time) %>%
    dplyr::summarize(volume = mean(.data$volume, na.rm = TRUE))

  evap_bay_a <-
    samples %>%
    dplyr::filter(experiment == "bay" & batch == "a") %>%
    dplyr::select(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$oxygen,
      .data$treatment
    ) %>%
    dplyr::distinct() %>%
    dplyr::left_join(evap_dup_bay, by = c("oxygen", "treatment"))

  evap_dup_hyp <-
    evap %>%
    dplyr::filter(.data$experiment == "05") %>%
    dplyr::group_by(.data$oxygen, .data$treatment, .data$time) %>%
    dplyr::summarize(volume = mean(.data$volume, na.rm = TRUE)) %>%
    dplyr::mutate(oxygen = replace(.data$oxygen, .data$oxygen == "0.5%", "0.2%"))

  evap_02 <-
    samples %>%
    dplyr::filter(experiment == "02") %>%
    dplyr::select(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$oxygen,
      .data$treatment
    ) %>%
    dplyr::distinct() %>%
    dplyr::left_join(evap_dup_hyp, by = c("oxygen", "treatment"))

  dplyr::bind_rows(evap, evap_02, evap_bay_a)
}

assemble_flux_measurements <- function(conc_clean, evap_clean) {
  abbreviations <-
    tibble::tibble(metabolite = unique(conc_clean$metabolite)) %>%
    dplyr::mutate(abbreviation = dplyr::case_when(
      .data$metabolite == "dna" ~ "cells",
      .data$metabolite == "glucose" ~ "glc",
      .data$metabolite == "asparagine" ~ "asn",
      .data$metabolite == "cystine" ~ "cyx",
      .data$metabolite == "glutamine" ~ "gln",
      .data$metabolite == "isoleucine" ~ "ile",
      .data$metabolite == "tryptophan" ~ "trp",
      TRUE ~ stringr::str_extract(.data$metabolite, "^[a-z]{3}")
    ))

  conc_clean %>%
    dplyr::left_join(evap_clean, by = c(
      "cell_type",
      "experiment",
      "batch",
      "date",
      "oxygen",
      "treatment",
      "time"
    )) %>%
    dplyr::left_join(abbreviations, by = "metabolite") %>%
    dplyr::filter(!is.na(.data$conc)) %>%
    dplyr::select(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$metabolite,
      .data$abbreviation,
      .data$detector,
      .data$type,
      .data$oxygen,
      .data$treatment,
      .data$time,
      .data$well,
      .data$conc,
      .data$volume
    ) %>%
    dplyr::mutate(
      treatment = factor(
        .data$treatment,
        levels = c("none", "DMSO", "BAY"),
        labels = c("None", "DMSO", "BAY")
      ),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%", "0.2%")),
      group = dplyr::case_when(
        experiment %in% c("02", "05", "bay") & treatment == "None" & oxygen == "21%" ~ "21%",
        experiment %in% c("02", "05", "bay") & treatment == "DMSO" ~ "DMSO",
        experiment %in% c("02", "05", "bay") & treatment == "BAY" ~ "BAY",
        experiment %in% c("02", "05", "bay") & oxygen == "0.5%" ~ "0.5%",
        experiment %in% c("02", "05", "bay") & oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY")),
      metabolite = replace(.data$metabolite, .data$metabolite == "dna", "cells"),
      nmol = .data$conc * .data$volume,
      abbreviation = toupper(.data$abbreviation)
    ) %>%
    dplyr::relocate(group, .before = time)
}

calculate_growth_rates <- function(growth_curves) {
  growth_m <- function(df) {
    fit <- MASS::rlm(log(conc) ~ time, data = df)
    names(fit$coefficients) <- c("X0", "mu")
    fit
  }

  growth_curves %>%
    dplyr::select(-c(.data$title, .data$plots)) %>%
    tidyr::unnest(c(.data$data)) %>%
    dplyr::filter(.data$time < 96) %>%
    dplyr::group_by(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$oxygen,
      .data$treatment
    ) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(.data$data, growth_m),
      summary = map(.data$model, broom::tidy)
    ) %>%
    tidyr::unnest(c(.data$summary)) %>%
    dplyr::select(-c(.data$std.error, .data$statistic)) %>%
    tidyr::pivot_wider(names_from = .data$term, values_from = .data$estimate) %>%
    dplyr::mutate(
      X0 = exp(.data$X0),
      group = dplyr::case_when(
        experiment %in% c("02", "05", "bay") & treatment == "None" & oxygen == "21%" ~ "21%",
        experiment %in% c("02", "05", "bay") & treatment == "DMSO" ~ "DMSO",
        experiment %in% c("02", "05", "bay") & treatment == "BAY" ~ "BAY",
        experiment %in% c("02", "05", "bay") & oxygen == "0.5%" ~ "0.5%",
        experiment %in% c("02", "05", "bay") & oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) %>%
    dplyr::arrange(desc(.data$experiment), .data$oxygen, .data$treatment) %>%
    dplyr::select(-c(.data$data, .data$model)) %>%
    dplyr::relocate(group, .before = X0)
}

plot_masses <- function(
  df,
  plot_function,
  ...
) {
  df %>%
    dplyr::group_by(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$metabolite,
      .data$abbreviation
    ) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      title = stringr::str_c(
        .data$metabolite,
        .data$cell_type,
        .data$experiment,
        .data$batch,
        .data$date,
        sep = "_"),
      plots = purrr::map2(
        .data$data,
        .data$title,
        plot_function,
        ...
      ))

}

geom_fit <- function(fit = "line", fo = NULL, method = NULL) {
  if (fit == "reg") {
    ggplot2::geom_smooth(
      method = method,
      formula = fo,
      se = FALSE
    )
  } else if (fit == "line") {
    ggplot2::stat_summary(
      fun = "mean",
      geom = "line"
    )
  }
}

plot_raw_curves <- function(data, title, y, xlab = "Time (h)", ylab, ...) {
  ggplot2::ggplot(data) +
    ggplot2::aes(
      x = time,
      y = {{y}},
      color = interaction(oxygen, treatment, sep = " | ")
    ) +
    ggplot2::geom_point(
      aes(shape = well),
      size = 3,
      alpha = 0.3
    ) +
    geom_fit(...) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8
    ) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      color = "condition"
    ) +
    wmo::theme_wmo()
}

plot_growth_curves <- function(flux_measurements) {

  flux_measurements %>%
    dplyr::filter(metabolite == "cells") %>%
    plot_masses(
      plot_raw_curves,
      y = conc,
      fit = "line",
      ylab = "Cell count"
    )

}

plot_degradation_curves <- function(flux_measurement) {

  flux_measurements %>%
    dplyr::filter(type == "empty") %>%
    plot_masses(
      plot_raw_curves,
      y = log(nmol),
      ylab = "ln(Mass (nmol))",
      fit = "reg",
      method = MASS::rlm,
      fo = y ~ x
    )

}

calculate_degradation_rates <- function(df) {
  degradation_m <- function(df) {
    fit <- MASS::rlm(log(nmol) ~ time, data = df)
    names(fit$coefficients) <- c("intercept", "k")
    fit
  }

  df %>%
    dplyr::select(-c(.data$title, .data$plots)) %>%
    tidyr::unnest(c(.data$data)) %>%
    dplyr::group_by(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$metabolite,
      .data$abbreviation,
      .data$oxygen,
      .data$treatment
    ) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(.data$data, degradation_m),
      summary = map(.data$model, broom::tidy)
    ) %>%
    tidyr::unnest(c(.data$summary)) %>%
    dplyr::filter(term == "k") %>%
    dplyr::select(-c(
      .data$term,
      .data$model,
      .data$data,
      .data$std.error,
      .data$statistic
    )) %>%
    dplyr::rename(k = .data$estimate)
}

clean_degradation_rates <- function(degradation_rates) {
  k <-
    degradation_rates %>%
    dplyr::group_by(.data$metabolite, .data$oxygen, .data$treatment) %>%
    wmo::remove_nested_outliers(k, remove = TRUE) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      ttest = purrr::map(.data$data, ~ t.test(.x$k, mu = 0)),
      summary = purrr::map(.data$ttest, broom::tidy)
    ) %>%
    tidyr::unnest(c(.data$summary)) %>%
    dplyr::filter(p.value < 0.01) %>%
    dplyr::select(
      .data$metabolite,
      .data$oxygen,
      .data$treatment,
      k = .data$estimate
    )

  hyp_02 <-
    k %>%
    dplyr::filter(oxygen == "0.5%") %>%
    dplyr::mutate(oxygen = forcats::fct_recode(oxygen, "0.2%" = "0.5%"))

  bay <-
    k %>%
    dplyr::filter(treatment == "DMSO") %>%
    dplyr::mutate(treatment = forcats::fct_recode(treatment, "BAY" = "DMSO"))

  dplyr::bind_rows(k, hyp_02, bay) %>%
    dplyr::mutate(k = -k) %>%
    dplyr::arrange(metabolite, oxygen, treatment)
}

plot_mass_curves <- function(flux_measurements) {
  flux_measurements %>%
    dplyr::filter(type == "cells" & metabolite != "cells") %>%
    plot_masses(
      plot_raw_curves,
      y = nmol,
      ylab = "ln(Mass (nmol))",
      fit = "line"
    )
}

plot_flux_curves <- function(mass_curves, k, growth_rates) {
  mass_curves %>%
    dplyr::select(-c(.data$title, .data$plots)) %>%
    tidyr::unnest(c(data)) %>%
    dplyr::left_join(k, by = c("metabolite", "oxygen", "treatment")) %>%
    dplyr::left_join(growth_rates, by = c(
      "cell_type",
      "experiment",
      "batch",
      "date",
      "oxygen",
      "treatment"
    )) %>%
    dplyr::mutate(
      k = tidyr::replace_na(.data$k, 0),
      x = exp((.data$mu + .data$k) * .data$time) - 1,
      y = .data$nmol * exp(.data$k * .data$time)
    ) %>%
    plot_masses(
      plot_raw_curves,
      y = y,
      xlab = expression(e^{(mu+k)*t} - 1),
      ylab = expression(M*e^{k*t}),
      fit = "reg",
      method = MASS::rlm,
      fo = y ~ x
    )
}

calculate_fluxes <- function(flux_curves) {
  flux_m <- function(df) {
    fit <- MASS::rlm(y ~ x, data = df)
    names(fit$coefficients) <- c("M0", "m")
    fit
  }

  flux_curves %>%
    dplyr::select(-c(.data$title, .data$plots)) %>%
    tidyr::unnest(c(.data$data)) %>%
    dplyr::filter(.data$time < 96) %>%
    dplyr::group_by(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$metabolite,
      .data$abbreviation,
      .data$oxygen,
      .data$treatment,
      .data$k,
      .data$X0,
      .data$mu
    ) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(.data$data, flux_m),
      summary = purrr::map(.data$model, broom::tidy)
    ) %>%
    tidyr::unnest(c(.data$summary)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(term == "m") %>%
    dplyr::select(-c(
      .data$term,
      .data$model,
      .data$data,
      .data$std.error,
      .data$statistic
    )) %>%
    dplyr::rename(m = .data$estimate) %>%
    dplyr::mutate(
      flux = m * (mu + k) / X0 * 1E6,
      group = dplyr::case_when(
        experiment %in% c("02", "05", "bay") & treatment == "None" & oxygen == "21%" ~ "21%",
        experiment %in% c("02", "05", "bay") & treatment == "DMSO" ~ "DMSO",
        experiment %in% c("02", "05", "bay") & treatment == "BAY" ~ "BAY",
        experiment %in% c("02", "05", "bay") & oxygen == "0.5%" ~ "0.5%",
        experiment %in% c("02", "05", "bay") & oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
      ) %>%
    dplyr::select(-c(.data$k, .data$X0, .data$mu, .data$m)) %>%
    dplyr::relocate(.data$metabolite, .data$abbreviation) %>%
    dplyr::relocate(.data$group, .after = .data$treatment) %>%
    dplyr::arrange(
      .data$metabolite,
      .data$oxygen,
      .data$treatment,
      .data$date
    )
}
