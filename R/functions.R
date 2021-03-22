# functions.R

# %nin% -------------------------------------------------------------------

#' Negate %in%
#'
#' Binary operator that returns a logical vector returning `TRUE` if there is
#' not a match for its left operand.
#'
#' @param x The values to be matched
#' @param y The values to be matched against
#'
"%nin%" <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}

# my_kable ----------------------------------------------------------------

my_kable <- function(data, ...) {
  kableExtra::kable(data, booktabs = TRUE, linesep = "", ...) %>%
    kableExtra::kable_styling(
      latex_options = c("hold_position"),
      font_size = 9
    )
}

# path_to_data ------------------------------------------------------------

path_to_data <- function(nm) {
  list.files(
    system.file(
      "extdata/",
      package = "Copeland.2021.hypoxia.flux"
    ),
    pattern = nm,
    full.names = TRUE,
    recursive = TRUE
  )
}

# path_to_reports ---------------------------------------------------------

path_to_reports <- function(nm) {
  path <-
    stringr::str_c(
      system.file(
        stringr::str_c("analysis"),
        package = "Copeland.2021.hypoxia.flux"
      ),
      "/",
      nm
    )
  path
}

# path_to_plots -----------------------------------------------------------

path_to_plots <- function(nm) {
  path <-
    stringr::str_c(
      system.file(
        stringr::str_c("analysis/"),
        package = "Copeland.2021.hypoxia.flux"
      ),
      "/figures/",
      nm
    )

  if (dir.exists(path)) unlink(path, recursive = TRUE)

  if (!dir.exists(path)) dir.create(path = path, recursive = TRUE)

  path

}

# read_multi_excel --------------------------------------------------------

read_multi_excel <- function(excel_file) {
  sheets <- readxl::excel_sheets(excel_file)
  purrr::map(sheets, ~readxl::read_excel(excel_file, sheet = .x)) %>%
    rlang::set_names(sheets)
}


# clean_technical_replicates ----------------------------------------------

clean_technical_replicates <- function(tbl) {
  tidyr::pivot_longer(
    data = tbl,
    cols = .data$a:.data$c,
    names_to = "replicate",
    values_to = "value"
  ) %>%
    dplyr::group_by(
      dplyr::across(-c(.data$replicate, .data$value))
    ) %>%
    dplyr::mutate(value = replace_outliers(value)) %>%
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE)) %>%
    dplyr::ungroup()
}

# make_std_curves ---------------------------------------------------------

make_std_curves <- function(df) {
  df %>%
    dplyr::filter(!is.na(.data$conc)) %>%
    dplyr::select(where(~all(!is.na(.)))) %>%
    dplyr::group_by(dplyr::across(-c(.data$conc, .data$value))) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      title = stringr::str_c(!!!rlang::syms(dplyr::groups(.)), sep = "_")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      model = furrr::future_map(
        .data$data,
        ~lm(value ~ conc, data = .x, na.action = modelr::na.warn)
      ),
      summary = furrr::future_map(.data$model, broom::glance),
      plots = furrr::future_map2(.data$data, .data$title, make_std_plots)
    ) %>%
    dplyr::group_by(
      dplyr::across(
        -c(.data$data, .data$title, .data$model, .data$summary, .data$plots)
      )
    )
}

# make_std_plots ----------------------------------------------------------

make_std_plots <- function(df, title = NULL) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = .data$conc,
      y = .data$value
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      color = "gray20",
      se = FALSE
    ) +
    ggplot2::geom_point(
      size = 3,
      alpha = 0.3,
      color = "blue"
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8,
      color = "blue"
    ) +
    ggplot2::labs(
      x = "Concentration",
      y = "Value",
      title = title
    )
}

# interp_data -------------------------------------------------------------

interp_data <- function(tbl, std) {
  tbl %>%
    dplyr::filter(is.na(.data$conc)) %>%
    dplyr::select(where(~all(!is.na(.)))) %>%
    dplyr::group_by(dplyr::across(dplyr::group_vars(std))) %>%
    tidyr::nest() %>%
    dplyr::left_join(dplyr::select(std, .data$model)) %>%
    dplyr::mutate(conc = purrr::map2(.data$data, .data$model, wmo::interpolate)) %>%
    tidyr::unnest(c(.data$data, .data$conc)) %>%
    dplyr::select(-c(.data$model, .data$value))
}

# interpolate -------------------------------------------------------------

interpolate <- function(new_df, model) {
  x <- stats::model.frame(model)[[deparse(model$terms[[3]])]]
  p <- polynom::polynomial(stats::coefficients(model))
  new_y <- as.list(new_df[[deparse(model$terms[[2]])]])
  new_x <- unlist(lapply(new_y, function(y) {
    roots <- solve(p, y)
    roots <- round(roots, digits = 8)
    root <- roots[which(Im(roots) == 0 & Re(roots) >= 0 & Re(roots) <= 1.25 * max(x))]
    ifelse(identical(root, numeric(0)), NA, Re(root))
  }))
  new_x
}

# clean_dna_per_cell ------------------------------------------------------

clean_dna_per_cell <- function(filename) {
  filename %>%
    read_multi_excel() %>%
    purrr::map(clean_technical_replicates) %>%
    dplyr::bind_rows(.id = "id") %>%
    tidyr::separate(id, c("cell_type", "volume"), sep = "_", convert = TRUE) %>%
    dplyr::mutate(cells = 1000 * cells)
}

# calculate_cells_per_dna -------------------------------------------------

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

# clean_flux_meta ---------------------------------------------------------

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

# assemble_flux_data ------------------------------------------------------

assemble_flux_data <- function(file_list) {
  nms <- sub("\\..*$", "", basename(file_list))
  sheets <- unique(unlist(purrr::map(file_list, readxl::excel_sheets)))
  data_list <- purrr::map(file_list, read_multi_excel)

  purrr::map(sheets, ~purrr::map(data_list, .x)) %>%
    rlang::set_names(sheets) %>%
    purrr::map(rlang::set_names, nms) %>%
    purrr::map(dplyr::bind_rows, .id = "experiment")
}

# clean_fluxes ------------------------------------------------------------

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

# print_plots -------------------------------------------------------------

print_plots <- function(
  plot_list,
  name_list,
  path_name,
  width = 20,
  height = 15
) {
  path <- path_to_plots(path_name)
  furrr::future_walk2(
    plot_list,
    name_list,
    ~ggplot2::ggsave(
      filename = stringr::str_c(.y, ".pdf"),
      path = path,
      plot = .x,
      device = cairo_pdf,
      width = width,
      height = height,
      units = "cm"
    )
  )
  invisible(path)
}

# clean_flux_std ----------------------------------------------------------

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

# fill_missing_fluxes -----------------------------------------------------

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

# filter_assays -----------------------------------------------------------

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

# assemble_evap_data ------------------------------------------------------

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

# fill_missing_evap -------------------------------------------------------

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

# assemble_flux_measurements ----------------------------------------------

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
      metabolite = replace(.data$metabolite, .data$metabolite == "dna", "cells"),
      nmol = .data$conc * .data$volume,
      abbreviation = toupper(.data$abbreviation)
    )
}

# calculate_growth_rates --------------------------------------------------

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
    dplyr::mutate(X0 = exp(.data$X0)) %>%
    dplyr::arrange(desc(.data$experiment), .data$oxygen, .data$treatment) %>%
    dplyr::select(-c(.data$data, .data$model))
}

# plot_growth_curves ------------------------------------------------------

plot_growth_curves <- function(data, title) {
  ggplot2::ggplot(data) +
    ggplot2::aes(
      x = time,
      y = conc,
      color = interaction(oxygen, treatment, sep = " | ")
    ) +
    ggplot2::geom_point(
      aes(shape = well),
      size = 3,
      alpha = 0.3
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      geom = "line"
    ) +
    ggplot2::labs(
      title = title,
      x = "Time (h)",
      y = "Cell count",
      color = "condition")
}

# plot_degradation_curves -------------------------------------------------

plot_degradation_curves <- function(data, title) {
  ggplot2::ggplot(data) +
    ggplot2::aes(
      x = time,
      y = log(nmol),
      color = interaction(oxygen, treatment)
    ) +
    ggplot2::geom_smooth(
      formula = y ~ x,
      method = MASS::rlm,
      se = FALSE
    ) +
    ggplot2::geom_point(
      aes(shape = well),
      size = 3,
      alpha = 0.3
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8
    ) +
    ggplot2::labs(
      title = title,
      x = "Time (h)",
      y = "ln(Mass (nmol))",
      color = "condition"
    ) +
    ggplot2::theme(legend.title = element_blank())
}

# plot_masses -------------------------------------------------------------

plot_masses <- function(
  df,
  plot_function
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
        plot_function
      ))

}


# calculate_degradation_rates ---------------------------------------------

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


# clean_degradation_rates -------------------------------------------------

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

# plot_mass_curves --------------------------------------------------------

plot_mass_curves <- function(data, title) {
  ggplot2::ggplot(data) +
    ggplot2::aes(
      x = time,
      y = nmol,
      color = interaction(oxygen, treatment)
    ) +
    ggplot2::geom_point(
      aes(shape = well),
      size = 3,
      alpha = 0.3
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 1,
      geom = "line",
    ) +
    ggplot2::labs(
      title = title,
      x = "Time (h)",
      y = "ln(Mass (nmol))",
      color = "condition"
    ) +
    ggplot2::theme(legend.title = element_blank())
}


# plot_flux_curves --------------------------------------------------------

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
    plot_masses(plot_flx_crvs)
}

# plot_flx_crvs -----------------------------------------------------------

plot_flx_crvs <- function(data, title) {
  ggplot2::ggplot(data) +
    ggplot2::aes(
      x = x,
      y = y,
      color = interaction(oxygen, treatment)
    ) +
    ggplot2::geom_smooth(
      formula = y ~ x,
      method = MASS::rlm,
      se = FALSE
    ) +
    ggplot2::geom_point(
      aes(shape = well),
      size = 3,
      alpha = 0.3
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8
    ) +
    ggplot2::labs(
      title = title,
      x = expression(e^{(mu+k)*t} - 1),
      y = expression(M*e^{k*t}),
      color = "condition"
    ) +
    ggplot2::theme(legend.title = element_blank())
}

# calculate_fluxes --------------------------------------------------------

calculate_fluxes <- function(flux_curves) {
  flux_m <- function(df) {
    fit <- MASS::rlm(y ~ x, data = df)
    names(fit$coefficients) <- c("M0", "m")
    fit
  }

  flux_curves %>%
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
    dplyr::mutate(flux = m * (mu + k) / X0 * 1E6) %>%
    dplyr::select(-c(.data$k, .data$X0, .data$mu, .data$m)) %>%
    dplyr::relocate(.data$metabolite, .data$abbreviation) %>%
    dplyr::arrange(
      .data$metabolite,
      .data$oxygen,
      .data$treatment,
      .data$date
    )
}


# calculate_ratios --------------------------------------------------------

calculate_ratios <- function(file_name) {
  readr::read_csv(file_name) %>%
    dplyr::filter(!is.na(area)) %>%
    tidyr::separate(
      .data$filename,
      into = c(NA, "window", "replicate"),
      sep = c(1, 2),
      convert = TRUE
    ) %>%
    tidyr::separate(
      .data$ion,
      into = c("metabolite", "isotope"),
      sep = " ",
      convert = TRUE
    ) %>%
    dplyr::mutate(carbons = case_when(
      .data$metabolite %in% c("citrate") ~ 6,
      .data$metabolite %in% c("2HG", "2OG", "glutamate", "glutamine") ~ 5,
      .data$metabolite %in% c("aspartate", "malate", "succinate") ~ 4,
      .data$metabolite %in% c("lactate", "pyruvate", "alanine", "serine") ~ 3
    )) %>%
    dplyr::filter(.data$window <= .data$carbons) %>%
    tidyr::pivot_wider(names_from = .data$isotope, values_from = .data$area) %>%
    dplyr::mutate(ratio = .data$M1 / .data$M0) %>%
    dplyr::group_by(.data$metabolite)
}


# import_qbias ------------------------------------------------------------

import_qbias <- function(file_list) {
  file_list %>%
    rlang::set_names(
      stringr::str_extract(basename(.), pattern = "(?<=_)\\w(?=\\.csv)")
    ) %>%
    purrr::map_dfr(calculate_ratios, .id = "batch") %>%
    dplyr::group_by(batch, metabolite) %>%
    dplyr::arrange(metabolite)
}


# calculate_predicted_ratios ----------------------------------------------

calculate_predicted_ratios <- function() {
  predicted_ratios <-
    isotope_library %>%
    dplyr::mutate(table = purrr::map2(
      formula,
      polarity,
      ~mzrtools::mz_iso_quant(molecule = .x, polarity = .y)[["prob_matrix"]])
    ) %>%
    dplyr::mutate(pred_ratio = purrr::map_dbl(table, ~ .x[[2, 1]]/.x[[1, 1]])) %>%
    dplyr::select(.data$metabolite, .data$pred_ratio)
}


# calculate_correction_factors --------------------------------------------

calculate_correction_factors <- function(qbias_ratios, pred_ratios) {
  qbias_ratios %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(.data$data, ~MASS::rlm(ratio ~ poly(window, 3), data = .x)),
      predict = purrr::map2(.data$model, .data$data, predict)
    ) %>%
    tidyr::unnest(c(data, predict)) %>%
    dplyr::select(.data$batch, .data$metabolite, .data$window, .data$carbons, .data$predict) %>%
    dplyr::distinct() %>%
    dplyr::left_join(pred_ratios, by = "metabolite") %>%
    dplyr::mutate(cf = .data$predict / .data$pred_ratio) %>%
    dplyr::select(.data$metabolite, "M" = .data$window, .data$carbons, .data$cf) %>%
    dplyr::filter(.data$M < .data$carbons) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(M = M + 1) %>%
    tidyr::pivot_wider(names_from = .data$M, values_from = .data$cf) %>%
    dplyr:: mutate(
      M0 = 1,
      M1 = 1 / .data$`1` * .data$M0,
      M2 = 1 / .data$`2` * .data$M1,
      M3 = 1 / .data$`3` * .data$M2,
      M4 = 1 / .data$`4` * .data$M3,
      M5 = 1 / .data$`5` * .data$M4,
      M6 = 1 / .data$`6` * .data$M5
    ) %>%
    dplyr::select(.data$batch, .data$metabolite, tidyselect::matches("M[0-9]+")) %>%
    tidyr::pivot_longer(
      cols = matches("M[0-9]+"),
      names_to = "M",
      values_to = "cf",
      values_drop_na = TRUE
    ) %>%
    dplyr::arrange(.data$batch, .data$metabolite)
}


# clean_mids --------------------------------------------------------------

clean_mids <- function(mid_files) {
  read_function <- function(filename) {
    readr::read_csv(filename) %>%
      replace(is.na(.), 0)
  }

  mid_files %>%
    rlang::set_names(sub("\\..*$", "", basename(.))) %>%
    purrr::map_dfr(read_function, .id = "experiment") %>%
    tidyr::pivot_longer(
      cols = .data$`2HG M0`:tidyselect::last_col(),
      names_to = "ion",
      values_to = "peak_area",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(
      batch = stringr::str_extract(.data$experiment, "a|b|c"),
      method = stringr::str_extract(.data$experiment, "fs|sim"),
      cell_type = stringr::str_extract(.data$experiment, "lf|pasmc"),
      date = stringr::str_extract(.data$experiment, "\\d{4}-\\d{2}-\\d{2}"),
      .before = .data$experiment
    ) %>%
    dplyr::select(-.data$experiment) %>%
    tidyr::separate(
      .data$id,
      c("tracer", "treatment", "time", "well"),
      sep = " "
    ) %>%
    tidyr::separate(.data$ion, into = c("metabolite", "isotope")) %>%
    dplyr::mutate(
      peak_area = tidyr::replace_na(peak_area, 0),
      oxygen = dplyr::case_when(
        .data$treatment %in% c("21%", "dmso", "bay") ~ "21%",
        .data$treatment == "0.5%" ~ "0.5%",
        .data$treatment == "0.2%" ~ "0.2%"
      ),
      treatment = dplyr::case_when(
        .data$treatment %in% c("21%", "0.5%", "0.2%") ~ "None",
        TRUE ~ toupper(.data$treatment)
      ),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%", "0.2%")),
      treatment = factor(.data$treatment, levels = c("None", "DMSO", "BAY")),
      time = as.double(.data$time),
      time = dplyr::case_when(
        .data$cell_type == "lf" ~ time * 24,
        .data$cell_type == "pasmc" ~ time * 12
      )
    ) %>%
    dplyr::relocate(.data$oxygen, .before = .data$treatment) %>%
    dplyr::select(-filename)
}


# correct_mid -------------------------------------------------------------

correct_mid <- function(mid_clean){
  correction_factors <-
    dplyr::mutate(correction_factors, method = "sim")

  mid_q <-
    mid_clean %>%
    dplyr::left_join(
      correction_factors,
      by = c("batch", "method", "metabolite", "isotope" = "M")
    ) %>%
    dplyr::mutate(
      area_corr = dplyr::case_when(
        .data$method == "sim" ~ peak_area * cf,
        .data$method == "fs" ~ peak_area * 1
      )
    ) %>%
    dplyr::select(-c(.data$peak_area, .data$cf))

  correction_matrices <-
    isotope_library %>%
    dplyr::mutate(
      matrix = purrr::map2(formula, polarity, mzrtools::mz_iso_quant),
      matrix = purrr::map(matrix, purrr::pluck, "prob_matrix")
    )

  mmult <- function(m, df) {
    mid <- df$mid
    if (nrow(m) > length(mid)) {
      m <- m[1:length(mid), 1:length(mid)]
    }
    mid_corr <- mzrtools::mz_iso_correct(m, mid)
    dplyr::bind_cols(df, mid_corr = mid_corr)
  }

  mid_na <-
    mid_q %>%
    dplyr::group_by(
      .data$method,
      .data$cell_type,
      .data$date,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$time,
      .data$well,
      .data$metabolite
    ) %>%
    dplyr::filter(.data$metabolite %nin% c("palmitate", "sedoheptulose")) %>%
    dplyr::mutate(mid = .data$area_corr / sum(.data$area_corr)) %>%
    dplyr::filter(.data$mid != "NaN") %>%
    dplyr::select(-.data$area_corr) %>%
    tidyr::nest() %>%
    dplyr::left_join(
      dplyr::select(correction_matrices, -c(.data$formula, .data$polarity)),
      by = "metabolite"
    ) %>%
    dplyr::mutate(data = purrr::map2(.data$matrix, .data$data, mmult)) %>%
    tidyr::unnest(c(.data$data)) %>%
    dplyr::filter(.data$isotope %in% stringr::str_c("M", 0:6)) %>%
    dplyr::select(
      .data$method:.data$metabolite,
      .data$isotope,
      .data$tracer:.data$well,
      .data$mid,
      .data$mid_corr
    )

}


# remove_mid_outliers -----------------------------------------------------

remove_mid_outliers <- function(mid_correct) {
  mid_correct %>%
    dplyr::group_by(
      dplyr::across(c(.data$method:.data$time, .data$metabolite:.data$isotope))
    ) %>%
    dplyr::mutate(
      outlier = replace_outliers(.data$mid_corr)
    ) %>%
    dplyr::summarise(mid = mean(mid_corr, na.rm = TRUE)) %>%
    dplyr::arrange(
      .data$cell_type,
      .data$metabolite,
      .data$oxygen,
      .data$treatment,
      .data$tracer
    )
}


# replace_outliers --------------------------------------------------------

replace_outliers <- function(vec) {
  if (mad(vec, na.rm = TRUE) == 0) return (vec)
  replace(
    vec,
    abs(vec - median(vec, na.rm = TRUE)) / mad(vec, na.rm = TRUE) > 2,
    NA
  )
}


# plot_mid_curves ---------------------------------------------------------

plot_mid_curves <- function(mids) {
  mids %>%
    dplyr::group_by(
      .data$method,
      .data$cell_type,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$metabolite
    ) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      title = stringr::str_c(
        .data$metabolite,
        .data$method,
        .data$tracer,
        .data$cell_type,
        stringr::str_replace(.data$oxygen, "%", ""),
        .data$treatment,
        sep = "_"
      ),
      plots = purrr::map2(
        .data$data,
        .data$title,
        ~ggplot2::ggplot(.x) +
          ggplot2::aes(
            x = .data$time,
            y = .data$mid,
            color = .data$isotope
          ) +
          ggplot2::geom_point(
            alpha = 0.3
          ) +
          ggplot2::stat_summary(
            fun = "mean",
            geom = "point",
            alpha = 0.8
          ) +
          ggplot2::stat_summary(
            fun = "mean",
            geom = "line",
            alpha = 0.8,
            show.legend = FALSE
          ) +
          ggplot2::scale_x_continuous(
            name = "Time (h)",
            breaks = seq(24, 96, by = 24)
          ) +
          ggplot2::scale_color_manual(values = viridis::viridis(7)) +
          ggplot2::labs(
            y = "Mole fraction",
            title = .y
          )
      )
    )
}


# clean_biomass -----------------------------------------------------------

clean_biomass <- function(biomass_file) {
  readr::read_csv(biomass_file) %>%
    dplyr::mutate(
      diff = post-pre,
      cell_mass = diff / cell_number * 1E12
    )
}


# calculate_biomass -------------------------------------------------------

calculate_biomass <- function(biomass_clean) {
  biomass_clean %>%
    dplyr::group_by(.data$cell_type, .data$date) %>%
    dplyr::mutate(cell_mass = replace_outliers(.data$cell_mass)) %>%
    dplyr::summarise(cell_mass = mean(.data$cell_mass, na.rm = TRUE)) %>%
    dplyr::mutate(cell_mass = replace_outliers(.data$cell_mass)) %>%
    dplyr::summarise(biomass = mean(.data$cell_mass, na.rm = TRUE))
}


#  calculate_biomass_equations --------------------------------------------

calculate_biomass_equations <- function(biomass) {
  biomass_coefs <- function(biomass, metabolites) {
    μmol_per_mass %>%
      dplyr::mutate(coefficient = .data$μmol_per_gDW * biomass) %>%
      dplyr::filter(.data$metabolite %in% metabolites) %>%
      dplyr::select(.data$metabolite, .data$coefficient) %>%
      dplyr::mutate(
        metabolite = dplyr::case_when(
          .data$metabolite %in% c("AcCoA") ~ stringr::str_c(metabolite, ".c"),
          TRUE ~ metabolite
        ))
  }

  biomass_eq <- function(biomass_coefs) {
    biomass_coefs %>%
      dplyr::mutate(
        combined = stringr::str_c(
          as.character(format(coefficient, digits = 0)),
          metabolite,
          sep = " ")
      ) %>%
      dplyr::pull(combined) %>%
      stringr::str_trim() %>%
      stringr::str_c(collapse = " + ") %>%
      stringr::str_c(" -> Biomass")
  }

  metabs <- c(
    "ALA",
    "ASP",
    "GLU",
    "GLN",
    "P5P",
    "AcCoA",
    # "LEU",
    # "ILE",
    # "VAL",
    "SER",
    "CYS",
    "GLY",
    "MEETHF",
    "CO2",
    "DHAP",
    # "ASN",
    # "ARG",
    "G6P"
  )

  biomass %>%
    dplyr::mutate(
      coefs = purrr::map(.data$biomass, biomass_coefs, metabs),
      eq = purrr::map(.data$coefs, biomass_eq)
    )
}


# save_biomass_equations --------------------------------------------------

save_biomass_equations <- function(biomass_equations) {
  path <- path_to_reports("modeling/matlab-input")

  purrr::walk2(
    biomass_equations$coefs,
    biomass_equations$cell_type,
    ~readr::write_csv(
      .x,
      file = file.path(path, stringr::str_c(.y, "_biomass.csv"))
    )
  )

  purrr::map_chr(
    biomass_equations$cell_type,
    ~ file.path(path, stringr::str_c(.x, "_biomass.csv"))
  )
}
