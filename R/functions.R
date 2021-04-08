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
  dir(
    system.file(
      "extdata/",
      package = "Copeland.2021.hypoxia.flux"
    ),
    pattern = nm,
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE
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


# path_to_manuscript ------------------------------------------------------

path_to_manuscript <- function(nm) {
  path <-
    stringr::str_c(
      system.file(
        stringr::str_c("manuscript"),
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

make_std_curves <- function(df, fo = NULL) {

  if (is.null(fo)){
    fo <- ~lm(value ~ conc, data = .x, na.action = modelr::na.warn)
  }

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
      model = furrr::future_map(.data$data, fo),
      summary = furrr::future_map(.data$model, ~broom::glance(.x)),
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
    dplyr::select(-.data$conc) %>%
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
  file_list[stringr::str_detect(file_list, "\\.csv$")] %>%
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

# calculate_biomass_equations --------------------------------------------

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

# write_matlab_input ------------------------------------------------------

write_matlab_input <- function(x, column, suffix) {
  path <- path_to_reports("modeling/matlab-input")
  column <- ensym(column)

  purrr::walk2(
    x[[column]],
    x[["cell_type"]],
    ~readr::write_csv(
      .x,
      file = file.path(path, stringr::str_c(.y, suffix))
    )
  )

  purrr::map_chr(
    x[["cell_type"]],
    ~ file.path(path, stringr::str_c(.x, suffix))
  )

}

# format_reactions --------------------------------------------------------

format_reactions <- function(reactions_file) {
  model_reactions <-
    readr::read_csv(reactions_file) %>%
    dplyr::filter(!is.na(.data$name)) %>%
    dplyr::mutate(
      pathway = factor(pathway, levels = c(
        "transport",
        "glycolysis",
        "pentose phosphate pathway",
        "anaplerosis",
        "tricarboxylic acid cycle",
        "amino acid metabolism",
        "biomass",
        "mixing",
        "dilution"
      ))
    )

  usethis::use_data(model_reactions, overwrite = TRUE)

  model_reactions
}

# format_fluxes -----------------------------------------------------------

format_fluxes <- function(growth_rates, fluxes) {
  growth <-
    growth_rates %>%
    dplyr::mutate(metabolite = "biomass") %>%
    dplyr::select(-.data$X0) %>%
    dplyr::rename(flux = .data$mu)

  fluxes_final <-
    fluxes %>%
    dplyr::filter(.data$experiment %in% c("05", "02", "bay")) %>%
    dplyr::filter(!(.data$experiment == "02" & metabolite == "pyruvate")) %>%
    dplyr::filter(!(.data$experiment == "02" & oxygen == "0.2%")) %>%
    dplyr::select(-.data$abbreviation) %>%
    dplyr::bind_rows(growth) %>%
    dplyr::group_by(.data$metabolite, .data$cell_type, .data$oxygen, .data$treatment) %>%
    wmo::remove_nested_outliers(flux, remove = TRUE)

  flux_names <- tribble(
    ~ metabolite, ~ name,
    "biomass", "BIOMASS",
    "alanine", "ALAR",
    "glucose", "GLUT",
    "glutamine", "GLNR",
    "glutamate", "GLUR",
    "lactate", "MCT",
    "pyruvate", "PYRR",
    "serine", "SERR",
    "cystine", "CYSR",
    "glycine", "GLYR",
    "aspartate", "ASPR"
  )

  model_fluxes <- c(
    "biomass",
    "alanine",
    "glucose",
    "glutamine",
    "glutamate",
    "lactate",
    "pyruvate",
    "serine",
    "cystine",
    "glycine",
    "aspartate"
  )

  fluxes_final %>%
    dplyr::filter(.data$metabolite %in% model_fluxes) %>%
    dplyr::group_by(.data$metabolite, .data$cell_type, .data$oxygen, .data$treatment) %>%
    dplyr::summarise(
      mean = mean(flux, na.rm = TRUE),
      se = sd(flux, na.rm = TRUE)/sqrt(n())
    ) %>%
    dplyr::filter(!is.nan(.data$mean)) %>%
    dplyr::left_join(flux_names, by = "metabolite") %>%
    dplyr::group_by(.data$cell_type) %>%
    tidyr::nest()

}

# format_mids -------------------------------------------------------------

format_mids <- function(mids) {
  model_metabolites <- c(
    "pyruvate",
    "lactate",
    "alanine",
    "2OG",
    "malate",
    "aspartate",
    "glutamate",
    "glutamine",
    "citrate",
    "serine",
    "FBP",
    "3PG"
  )

  mids_filtered <-
    mids %>%
    dplyr::filter(.data$metabolite %in% model_metabolites) %>%
    dplyr::filter(.data$time != 96)

  mids_lf_new <-
    mids_filtered %>%
    dplyr::filter(
      .data$cell_type == "lf" &
        .data$date %in% c(
          "2018-11-15", "2018-11-20", "2018-11-25", "2018-12-16",
          "2018-11-11", "2018-11-16",
          "2019-05-06", "2019-05-10", "2019-05-14"
        )
    ) %>%
    dplyr::filter(
      .data$method == "sim" |
        (.data$method == "fs" & .data$metabolite %in% c("FBP", "3PG"))
    )

  mids_pasmc_05 <-
    mids_filtered %>%
    dplyr::filter(
      .data$cell_type == "pasmc" &
        .data$oxygen %in% c("21%", "0.5%")
    ) %>%
    dplyr::filter(
      .data$method == "sim" |
        (.data$method == "fs" & .data$metabolite %in% c("FBP", "3PG"))
    )

  dplyr::bind_rows(mids_lf_new, mids_pasmc_05) %>%
    dplyr::group_by(
      .data$method,
      .data$cell_type,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$metabolite,
      .data$time,
      .data$isotope
    ) %>%
    wmo::remove_nested_outliers(mid, remove = TRUE) %>%
    # tidyr::nest() %>%
    # dplyr::mutate(data = purrr::map(.data$data, remove_outliers, "mid", remove = TRUE)) %>%
    # unnest(c(data)) %>%
    dplyr::summarise(
      mean = mean(.data$mid, na.rm = TRUE),
      se = sd(.data$mid, na.rm = TRUE)/sqrt(dplyr::n())
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      metabolite = replace(.data$metabolite, .data$metabolite == "glutamine", "GLN"),
      metabolite = replace(.data$metabolite, .data$metabolite == "2OG", "AKG"),
      metabolite = stringr::str_sub(.data$metabolite, 1, 3),
      metabolite = toupper(.data$metabolite)
    ) %>%
    dplyr::select(-.data$method) %>%
    dplyr::arrange(
      .data$metabolite,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$time,
      .data$isotope
    ) %>%
    dplyr::group_by(cell_type) %>%
    tidyr::nest()
}

# plot_growth_curve -------------------------------------------------------

plot_growth_curve <- function(
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  clr = c("oxygen", "treatment")
) {

  flux_measurements %>%
    dplyr::filter(
      cell_type == cell &
        experiment == exp &
        metabolite == "cells" &
        time < 96
    ) %>%
    dplyr::group_by(
      date,
      oxygen,
      treatment,
      time
    ) %>%
    dplyr::summarize(count = mean(conc, na.rm = TRUE)) %>%
    plot_time_lines(
      y = count/1000,
      ylab = expression(paste("Cell count (x", 10^3, ")")),
      clr = clr
    )
}

# plot_growth_rates -------------------------------------------------------

plot_growth_rates <- function(
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  clr = c("oxygen", "treatment"),
  xaxis = clr,
  annot
) {

  if(missing(annot)) {
    annot <-
      tibble::tibble(
        p = t.test(
          mu ~ !!rlang::ensym(clr),
          data = growth_rates,
          subset = c(cell_type == cell & experiment == exp),
          paired = TRUE
        )$p.value,
        y_pos = Inf,
        x_pos = 1.5,
        vjustvar = 1.5) %>%
      dplyr::mutate(lab = annot_p(p))
  }

  growth_rates %>%
    dplyr::filter(cell_type == cell & experiment == exp) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = .data[[xaxis]],
      y = mu
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = .data[[clr]]),
      geom = "col",
      fun = "mean",
      width = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = "mean_se",
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    # ggbeeswarm::geom_beeswarm(
    #   ggplot2::aes(color = .data[[clr]]),
    #   priority = "random",
    #   alpha = 0.4,
    #   cex = 4,
    #   size = 1,
    #   dodge.width = 0.7,
    #   show.legend = FALSE
    # ) +
    # ggplot2::stat_summary(
    #   color = "black",
  #   ggplot2::aes(group = .data[[clr]]),
  #   geom = "crossbar",
  #   fun = "mean",
  #   fatten = 1,
  #   width = 0.4,
  #   position = ggplot2::position_dodge(width = 0.7),
  #   show.legend = FALSE
  # ) +
  ggplot2::geom_text(
    data = annot,
    ggplot2::aes(
      label = lab,
      x = x_pos,
      y = y_pos,
      vjust = vjustvar
    )
  ) +
    ggplot2::labs(
      x = "Treatment",
      y = "Growth rate (/h)"
    ) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::coord_cartesian(ylim = c(0, 0.03)) +
    theme_plots()
}

# annot_p -----------------------------------------------------------------

annot_p <- function(num) dplyr::if_else(num < 0.05, "*", NA_character_)

# clean_viability ---------------------------------------------------------

clean_viability <- function(viability_file) {
  readr::read_csv(viability_file) %>%
    dplyr::mutate(
      viability = 100 * live / (dead + live),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) %>%
    dplyr::group_by(time, oxygen) %>%
    wmo::remove_nested_outliers(viability, remove = TRUE)
}

# plot_blot ---------------------------------------------------------------

plot_blot <- function(blot_image) {

  # blot_image <- magick::image_read_pdf(blot_image)
  blot_image <- magick::image_read(blot_image)

  cowplot::ggdraw() +
    cowplot::draw_image(
      blot_image,
      scale = 1.3,
      hjust = 0.15,
      vjust = 0.1
    ) +
    wmo::theme_wmo(
      base_family = "Calibri",
      base_size = 8
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      panel.border = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold")
    )
}


# plot_densities ----------------------------------------------------------

plot_densities <- function(
  densities,
  exp,
  prot = c("ldha", "hif1a"),
  ylab = prot,
  clr = c("oxygen", "treatment")
) {

  densities %>%
    dplyr::filter(experiment == exp) %>%
    dplyr::filter(protein == prot) %>%
    plot_time_lines(y = fold_change, ylab = ylab, clr = clr)
}

# normalize_qpcr ----------------------------------------------------------

normalize_qpcr <- function(raw_mrna) {
  raw_mrna %>%
    dplyr::mutate(gene = tolower(gene)) %>%
    dplyr::group_by(dplyr::across(c(experiment:gene))) %>%
    dplyr::summarize(ct = mean(ct, na.rm = TRUE)) %>%
    tidyr::pivot_wider(names_from = gene, values_from = ct) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dplyr::across(-c(experiment:time, actin), ~ . - actin)
    ) %>%
    dplyr::select(-actin) %>%
    tidyr::pivot_longer(
      -c(experiment:time),
      names_to = "rna",
      values_to = "dct"
    ) %>%
    dplyr::filter(!is.na(dct)) %>%
    dplyr::group_by(.data$experiment, .data$rna) %>%
    dplyr::mutate(
      ddct = dct - mean(
        dct[oxygen == min(oxygen) & treatment %in% c("None", "DMSO") & time == 0]
      ),
      log2fc = 2 ^ -ddct
    ) %>%
    dplyr::group_by(experiment, oxygen, treatment, time, rna) %>%
    wmo::remove_nested_outliers(log2fc, remove = TRUE)
}

# plot_mrna ---------------------------------------------------------------

plot_mrna <- function(
  mrna,
  exp,
  gene = c("ldha", "hif1a"),
  ylab = prot,
  clr = c("oxygen", "treatment")
) {

  mrna %>%
    dplyr::filter(rna == gene) %>%
    dplyr::filter(experiment == exp) %>%
    plot_time_lines(y = log2fc, ylab = ylab, clr = clr)
}

# plot_high_fluxes --------------------------------------------------------

plot_high_fluxes <- function(
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  clr = c("oxygen", "treatment"),
  xaxis = clr,
  annot
) {

  fluxes %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
        metabolite %in% c("lactate", "glucose")
    ) %>%
    dplyr::left_join(annot, by = c("metabolite", "abbreviation", "cell_type", "experiment")) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = .data[[clr]]),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = .data[[clr]]),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    # ggbeeswarm::geom_beeswarm(
    #   ggplot2::aes(
    #     color = .data[[clr]]),
    #   dodge.width = 0.9,
    #   priority = "random",
    #   alpha = 0.4,
    #   size = 1,
    #   cex = 4,
    #   show.legend = FALSE
    # ) +
    # ggplot2::stat_summary(
  #   color = "black",
  #   aes(group = .data[[clr]]),
  #   geom = "crossbar",
  #   fun = "mean",
  #   fun.min = "mean",
  #   fun.max = "mean",
  #   fatten = 1,
  #   width = 0.7,
  #   position = position_dodge(width = 0.9),
  #   show.legend = FALSE
  # ) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = y_pos,
      vjust = vjust,
      label = pval
    )
  ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Flux (fmol/cell/h)"
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    # ggplot2::coord_cartesian(ylim = c(-750, 1250)) +
    theme_plots()
}

# plot_low_fluxes ---------------------------------------------------------

plot_low_fluxes <- function(
  cell = c("lf", "pasmc"),
  exp = c("02", "05", "bay"),
  clr = c("oxygen", "treatment"),
  xaxis = clr,
  annot
) {

  fluxes %>%
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exp &
        !(metabolite %in% c("lactate", "glucose"))
    ) %>%
    dplyr::left_join(annot, by = c("metabolite", "abbreviation", "cell_type", "experiment")) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = .data[[clr]]),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = .data[[clr]]),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    # ggplot2::stat_summary(
    #   ggplot2::aes(color = .data[[clr]]),
    #   geom = "pointrange",
    #   fun.data = "mean_se",
    #   fatten = 1,
    #   position = position_dodge(width = 0.5),
    #   show.legend = TRUE
    # ) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = y_pos,
        vjust = vjust,
        label = pval
      )
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    # ggplot2::coord_cartesian(ylim = c(-125, 50)) +
    ggplot2::scale_y_continuous(
      trans = ggallin::pseudolog10_trans,
      breaks = c(-100, -10, 0, 10, 100),
      limits = c(-250, 250)
    ) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

# plot_viability ----------------------------------------------------------

plot_viability <- function(viability) {

  plot_time_lines(
    viability,
    y = viability,
    ylab = "Cell viability (%)",
    clr = "oxygen"
  )
}

# read_data ---------------------------------------------------------------

read_data <- function(data_files) {
  data_files[stringr::str_detect(data_files, "\\.csv$")] %>%
    rlang::set_names(stringr::str_extract(., "(lf|pasmc)_(02|05|bay)")) %>%
    purrr::map_dfr(read_csv, .id = "experiment") %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%", "0.2%"), ordered = TRUE),
      treatment = factor(treatment, levels = c("None", "DMSO", "BAY"), ordered = TRUE)
    )
}

# normalize_densities -----------------------------------------------------

normalize_densities <- function(blot_raw) {
  blot_raw %>%
    dplyr::filter(.data$time < 96) %>%
    dplyr::group_by(.data$experiment, .data$gel) %>%
    dplyr::mutate(
      blot_norm = blot / mean(blot, na.rm = TRUE),
      hif1a_norm = hif1a / mean(hif1a, na.rm = TRUE),
      ldha_norm = ldha / mean(ldha, na.rm = TRUE),
      hif1a_ratio = hif1a_norm / blot_norm,
      ldha_ratio = ldha_norm / blot_norm
    ) %>%
    dplyr::select(experiment:time, hif1a_ratio, ldha_ratio) %>%
    tidyr::pivot_longer(
      contains("ratio"),
      names_to = "protein",
      values_to = "density"
    ) %>%
    tidyr::separate(protein, into = c("protein", NA), sep = "_") %>%
    dplyr::group_by(experiment, protein) %>%
    dplyr::mutate(
      fold_change = density /
        mean(density[oxygen == min(oxygen) & treatment %in% c("None", "DMSO") & time == 0])
    ) %>%
    dplyr::group_by(experiment, oxygen, treatment, time, protein) %>%
    wmo::remove_nested_outliers(fold_change, remove = TRUE)
}

# theme_plots -------------------------------------------------------------

theme_plots <- function() {
  wmo::theme_wmo(
    base_family = "Calibri",
    base_size = 8
  ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 1))
    )
}

# plot_time_lines ---------------------------------------------------------

plot_time_lines <- function(
  df,
  y,
  ylab = prot,
  clr = c("oxygen", "treatment")
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
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots()
}

# annot_fluxes ------------------------------------------------------------

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

# arrange_fluxes ----------------------------------------------------------

arrange_fluxes <- function(a, b, c, d, e, f, g, h, i, j) {
  layout <- "
  abc
  def
  ghi
  jjj
"

  a + b + c + d + e + f + g + h + i + j +
    patchwork::plot_layout(design = layout) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(plot.margin = margin(-3, -3, -3, -3))
    )
}

# write_figures -----------------------------------------------------------

write_figures <- function(plot, filename, width, height) {
  path <- system.file("manuscript/figures", package = "Copeland.2021.hypoxia.flux")

  ggplot2::ggsave(
    plot,
    filename = filename,
    path = path,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )

  stringr::str_c(path, "/", filename)

}

# plot_cells_per_dna ------------------------------------------------------

plot_cells_per_dna <- function(dna_per_cell_clean) {
  dna_per_cell_clean %>%
    dplyr::filter(.data$volume == 200 & cells < 400000) %>%
    ggplot2::ggplot() +
    ggplot2::facet_wrap(~cell_type, labeller = ggplot2::as_labeller(toupper)) +
    ggplot2::aes(
      x = cells/1000,
      y = conc
    ) +
    ggplot2::geom_smooth(
      formula = y ~ 0 + x,
      method = "lm",
      color = clrs[[2]],
      size = 0.5,
      se = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
      fun.data = "mean_se",
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      fill = "black",
      size = 2,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = expression(paste("Cell count (x", 10^3, ")")),
      y = "DNA (ng)"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, NA)) +
    theme_plots()
}

# clean_dna_count_hypoxia -------------------------------------------------

clean_dna_count_hypoxia <- function(dna_count_hypoxia_file) {
  df <-
    readr::read_csv(dna_count_hypoxia_file) %>%
    dplyr::mutate(oxygen = factor(oxygen, levels = c("21%", "0.5%"))) %>%
    clean_technical_replicates()

  std <-
    df %>%
    dplyr::filter(!is.na(conc)) %>%
    lm(value ~ conc, data = .)

  df %>%
    dplyr::filter(is.na(conc)) %>%
    dplyr::mutate(conc = interpolate(., std)) %>%
    dplyr::select(-value)
}

# plot_dna_count_hypoxia --------------------------------------------------

plot_dna_count_hypoxia <- function(dna_count_hypoxia) {

  ggplot2::ggplot(dna_count_hypoxia) +
    ggplot2::aes(
      x = count / 1000,
      y = conc,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      size = 0.5,
      formula = y ~ x,
      se = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      pch = 21,
      color = "white",
      # alpha = 0.3,
      size = 2,
      show.legend = FALSE) +
    ggplot2::labs(
      x = expression(paste("Cell count (×", 10^3, ")")),
      y = "DNA (ng)"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 300, 100)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::coord_cartesian(ylim = c(0, NA), xlim = c(0, 300)) +
    theme_plots()
}

# plot_evap_data ----------------------------------------------------------

plot_evap_data <- function(evap_clean) {
  evap_clean %>%
    dplyr::filter(experiment == "05" & cell_type == "lf") %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = volume,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
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
      y = "Volume (mL)"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, 24)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    theme_plots()
}


# plot_k ------------------------------------------------------------------

plot_k <- function(degradation_rates, k) {

  annot <-
    k %>%
    dplyr::select(-k) %>%
    dplyr::mutate(
      label = "*",
      ypos = Inf,
      vjust = 1.5
    )

  degradation_rates %>%
    dplyr::mutate(
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" ~ "0.5%",
        treatment == "DMSO" ~ "DMSO"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO"))
    ) %>%
    dplyr::left_join(annot, by = c("metabolite", "oxygen", "treatment")) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), k),
      y = k
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        color = group,
        label = label,
        y = ypos,
        vjust = vjust
      ),
      position = ggplot2::position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Rate constant (/h)",
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_color_manual(
      values = clrs
    ) +
    ggplot2::scale_fill_manual(
      values = clrs
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
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
    patchwork::plot_layout(design = layout) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(plot.margin = margin(-3, -3, -3, -3))
    )
}

# plot_mids ---------------------------------------------------------------

plot_mids <- function(df) {

  tracer_labels <-
    c(expression(paste("[1,2-"^13, "C"[2], "] glucose")),
      expression(paste("[U-"^13, "C"[6], "] glucose")),
      expression(paste("[U-"^13, "C"[5], "] glutamine")),
      expression(paste("[U-"^13, "C"[3], "] lactate")))

  tracer_levels <-
    c("glc2", "glc6", "q5", "lac3")

  df %>%
    dplyr::mutate(
      tracer = factor(tracer, levels = tracer_levels, labels = tracer_labels),
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" & treatment == "None" ~ "0.5%",
        oxygen == "21%" & treatment == "DMSO" ~ "DMSO",
        oxygen == "21%" & treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY"))
    ) %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = isotope, y = mean, fill = group) +
    ggplot2::facet_grid(
      tracer ~ metabolite,
      labeller = ggplot2::label_parsed,
      # switch = "y",
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge()
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = mean - se,
        ymax = mean + se
      ),
      position = ggplot2::position_dodge(width = 0.9)
    ) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL
    ) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}


# plot_manuscript_mids ----------------------------------------------------

plot_manuscript_mids <- function(model_mids) {
  model_mids %>%
    dplyr::filter(cell_type == "lf") %>%
    tidyr::unnest(c(data)) %>%
    dplyr::mutate(
      metabolite = factor(
        metabolite,
        levels = c("FBP", "3PG", "PYR", "CIT", "AKG", "MAL"),
        labels = c("FBP", "`3PG`", "PYR", "CIT", "AKG", "MAL")
      )
    ) %>%
    dplyr::filter(tracer != "lac3" & time == 72 & !is.na(metabolite)) %>%
    plot_mids()
}

# plot_lf_mids ------------------------------------------------------------

plot_lf_mids <- function(model_mids) {
  model_mids %>%
    dplyr::filter(cell_type == "lf") %>%
    tidyr::unnest(c(data)) %>%
    dplyr::mutate(
      metabolite = factor(
        metabolite,
        levels = c("ALA", "ASP", "GLU", "GLN", "LAC", "SER")
      )
    ) %>%
    dplyr::filter(tracer != "lac3" & time == 72 & !is.na(metabolite)) %>%
    plot_mids()
}


# plot_pasmc_mids ---------------------------------------------------------

plot_pasmc_mids <- function(model_mids) {
  model_mids %>%
    dplyr::filter(cell_type == "pasmc") %>%
    tidyr::unnest(c(data)) %>%
    dplyr::mutate(
      metabolite = factor(
        metabolite,
        levels = c(
          "FBP",
          "3PG",
          "PYR",
          "ALA",
          "SER",
          "LAC",
          "CIT",
          "AKG",
          "GLN",
          "GLU",
          "MAL",
          "ASP"
        )
      ),
      metabolite = forcats::fct_recode(metabolite, "`3PG`" = "3PG")
    ) %>%
    dplyr::filter(tracer != "lac3" & time == 48 & !is.na(metabolite)) %>%
    plot_mids() +
    ggplot2::facet_grid(
      metabolite ~ tracer,
      labeller = ggplot2::label_parsed
    )
}

# clean_model_fluxes ------------------------------------------------------

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


# plot_labeling_rate ------------------------------------------------------

plot_labeling_rate <- function(mids) {
  df <-
    mids %>%
    dplyr::filter(
      .data$method == "sim" &
        .data$cell_type == "lf" &
        .data$tracer == "glc6" &
        .data$metabolite %in% c("pyruvate") &
        .data$isotope == "M0"
    ) %>%
    dplyr::mutate(
      labeled = 1 - mid,
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" ~ "0.5%",
        treatment == "DMSO" ~ "DMSO",
        treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY")),
      idx = dplyr::case_when(
        group %in% c("21%", "DMSO") ~ 1,
        group %in% c("0.5%", "BAY") ~ 2
      )
    )

  a <-
    df %>%
    dplyr::group_by(.data$time, .data$group) %>%
    wmo::remove_nested_outliers(labeled, remove = TRUE) %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = labeled,
      color = group,
      fill = group
    ) +
    ggplot2::geom_smooth(
      method = "nls",
      formula = y ~ SSasympOrig(x, asym, lrc),
      fullrange = TRUE,
      se = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "linerange",
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
    ggplot2::scale_x_continuous(
      breaks = seq(0, 72, 24),
      limits = c(0, 72)
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::labs(
      x = "Time (h)",
      y = "Fraction labeled",
      color = NULL,
      fill = NULL
    ) +
    theme_plots()

  k <-
    df %>%
    dplyr::group_by(.data$group) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      m = purrr::map(data, ~nls(labeled ~ SSasympOrig(time, asym, lrc), data = .x)),
      s = purrr::map(m, broom::tidy)
    ) %>%
    tidyr::unnest(c(s)) %>%
    dplyr::filter(term == "lrc") %>%
    dplyr::mutate(
      k = exp(estimate),
      experiment = dplyr::case_when(
        group %in% c("21%", "0.5%") ~ "hypoxia",
        group %in% c("DMSO", "BAY") ~ "bay"
      ),
      sem = std.error / abs(estimate) * k
    )

  b <-
    k %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = group,
      y = k,
      fill = group
    ) +
    ggplot2::geom_col(
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = k - sem,
        ymax = k + sem
      ),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::labs(
      x = NULL,
      y = "Rate",
      color = NULL,
      fill = NULL
    ) +
    theme_plots()

  list(curve = a, rate = b)
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
    patchwork::plot_layout(design = layout) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(plot.margin = margin(-3, -3, -3, -3))
    )
}

# calculate_flux_differences ----------------------------------------------

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

# normalize_fluxes --------------------------------------------------------

normalize_fluxes <- function(map_fluxes, reference_flux) {

  map_fluxes %>%
    dplyr::group_by(.data$cell_type, .data$treatment) %>%
    dplyr::mutate(
      lb = lb / flux[id == reference_flux],
      ub = ub / flux[id == reference_flux],
      flux = flux / flux[id == reference_flux]
    )

}

# assemble_flux_differences -----------------------------------------------

assemble_flux_differences <- function(map_fluxes) {

  glc_norm <- normalize_fluxes(map_fluxes, "GLUT")
  growth_norm <- normalize_fluxes(map_fluxes, "BIOMASS")

  tibble::tibble(
    map_fluxes = list(map_fluxes, map_fluxes, glc_norm, glc_norm, growth_norm, growth_norm),
    cell = rep(c("lf", "lf"), 3),
    control = rep(c("21%", "DMSO"), 3),
    experiment = rep(c("0.5%", "BAY"), 3)
  ) %>%
    purrr::pmap_dfr(calculate_flux_differences, .id = "normalization") %>%
    dplyr::mutate(normalization = dplyr::case_when(
      normalization %in% 1:2 ~ "none",
      normalization %in% 3:4 ~ "glucose",
      normalization %in% 5:6 ~ "growth"
    ))

}

# make_grid ---------------------------------------------------------------

make_grid <- function(arrow, lhs, rhs) {
  forward <- tidyr::expand_grid(substrate = lhs, product = rhs)
  reverse <- tidyr::expand_grid(substrate = rhs, product = lhs)

  if (arrow == "->") {
    forward
  } else {
    reverse
  }
}


# parse_eq ----------------------------------------------------------------

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

# make_graph --------------------------------------------------------------

make_graph <- function(map_flux_differences, nodes, treat, normalizer) {

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
    dplyr::filter(normalization == normalizer & treatment == treat) %>%
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

# set_graph_style ---------------------------------------------------------

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

# plot_ratio_network ------------------------------------------------------

plot_ratio_network <- function(graph, caption) {
  fold <- c(Inf, 5, 3, 2, 1.5, 1.1, 0.91, 0.67, 0.5, 0.33, 0.2, 0)

  cuts <- log(fold, base = 2)

  lvls <-
    tidygraph::activate(graph, edges) %>%
    dplyr::pull(ratio) %>%
    cut(
      cuts,
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

  set_gr_style()

  ggraph::ggraph(
    graph,
    circular = FALSE,
    x = tidygraph::.N()$x,
    y = tidygraph::.N()$y
  ) +
    ggraph::geom_edge_parallel(
      ggplot2::aes(
        color = cut(ratio, cuts, include.lowest = TRUE),
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
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = caption
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.position = c(0.07, 0.2),
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

# arrange_m4 --------------------------------------------------------------

arrange_m4 <- function(hypoxia_graph_ratio_plot, bay_graph_ratio_plot, m4c) {

  a <- hypoxia_graph_ratio_plot
  b <- bay_graph_ratio_plot

  top <-
    a + b +
    plot_layout(guides = "collect") &
    theme(legend.text.align = 1)

  # layout <- "
  # ab
  # cc
  # "

  top / m4c +
    plot_annotation(tag_levels = "A")
}


# plot_lactate_mids -------------------------------------------------------

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
    plot_mids() + ggplot2::facet_wrap(~ metabolite, scales = "free_x", nrow = 2)
}

# format_time_course_mids -------------------------------------------------

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

# plot_mid_time_course ----------------------------------------------------

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

# plot_normoxia_network ---------------------------------------------------

plot_normoxia_network <- function(hypoxia_graph) {

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

# arrange_s6 --------------------------------------------------------------

arrange_s6 <- function(a, b, c, d) {
  layout <- "
    ab
    cd
    cd"

  a + b + c + d +
    patchwork::plot_layout(design = layout) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(plot.margin = margin(-3, -3, -3, -3))
    )
}


# plot_twoby_fluxes -------------------------------------------------------

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
        color = treatment,
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
      data = dplyr::filter(annot, metabolite == metab),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      color = "Treatment"
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    # ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(NA, 0.2))) +
    # ggplot2::coord_cartesian(ylim = c(-750, 1250)) +
    theme_plots()

}

# analyze_twoby_fluxes ----------------------------------------------------

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
      lab = dplyr::if_else(adj.p.value < 0.05, "*", "")
    )

  list(data = df, annot = annot)

}

# clean_nad ---------------------------------------------------------------

clean_nad <- function(nad_data) {
  nad_data %>%
    dplyr::bind_rows(.id = "metabolite") %>%
    dplyr::rename(rep = replicate) %>%
    clean_technical_replicates() %>%
    tidyr::separate(.data$experiment, c(NA, "date"), "_")
}

# finalize_nad ------------------------------------------------------------

finalize_nad <- function(nad_interp, cells_per_dna) {

  x <-
    dplyr::filter(cells_per_dna, cell_type == "lf" & volume == 200) %>%
    dplyr::pull(slope)

  nad_interp %>%
    dplyr::group_by(metabolite, date, oxygen, treatment, nucleotide) %>%
    dplyr::summarise(conc = mean(conc)) %>%
    dplyr::mutate(
      conc = dplyr::case_when(
        metabolite == "dna" ~ conc * x,
        metabolite == "nad" ~ conc * 720
      ),
      nucleotide = replace(nucleotide, is.na(nucleotide), "Count")
    ) %>%
    tidyr::pivot_wider(-c(metabolite), names_from = nucleotide, values_from = conc) %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("none", "DMSO", "BAY")),
      treatment = forcats::fct_recode(treatment, "None" = "none"),
      Ratio = NADH/NAD,
      dplyr::across(tidyselect::contains("NAD"), ~. / Count * 1000)
    ) %>%
    dplyr::select(-Count) %>%
    tidyr::pivot_longer(c(NAD, NADH, Ratio), names_to = "measurement", values_to = "value") %>%
    dplyr::arrange(measurement, oxygen, treatment)

}

# plot_nad ----------------------------------------------------------------

plot_nad <- function(nad_final) {

  df <-
    nad_final %>%
    dplyr::filter(treatment != "None") %>%
    dplyr::mutate(
      measurement = factor(
        measurement,
        levels = c("NAD", "NADH", "Ratio"),
        labels = c("NAD (nmol/cell)", "NADH (nmol/cell)", "NADH/NAD ratio")
      )
    ) %>%
    dplyr::group_by(measurement, oxygen, treatment)

  annot <-
    df %>%
    dplyr::group_by(measurement) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(value ~ treatment * oxygen + (1 | date), data = .x)),
      s = purrr::map(
        m,
        ~emmeans::emmeans(.x, ~ treatment * oxygen) %>%
          pairs(simple = "each", combine = TRUE) %>%
          broom::tidy()
      )
    ) %>%
    tidyr::unnest(c(s)) %>%
    dplyr::select(measurement, oxygen, treatment, pval = adj.p.value) %>%
    dplyr::mutate(
      lab = dplyr::if_else(pval < 0.05, "*", ""),
      y = Inf,
      vjust = 1.5
    )

  annot1 <-
    dplyr::filter(annot, oxygen != ".") %>%
    dplyr::mutate(
      treatment = "BAY",
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )
  annot2 <-
    dplyr::filter(annot, treatment != ".") %>%
    dplyr::mutate(
      oxygen = "21%",
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )

  ggplot2::ggplot(df) +
    ggplot2::facet_wrap(
      ~ measurement,
      nrow = 1,
      scales = "free_y",
      strip.position = "left"
    ) +
    ggplot2::aes(
      x = treatment,
      y = value,
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
      show.legend = FALSE,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = annot1,
      ggplot2::aes(
        color = oxygen,
        y = y,
        label = lab,
        vjust = vjust
      ),
      position = ggplot2::position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot2,
      ggplot2::aes(
        y = y,
        label = lab,
        vjust = vjust
      ),
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "Oxygen",
      color = "Oxygen"
    ) +
    # ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(NA, 0.2))) +
    # ggplot2::coord_cartesian(ylim = c(-750, 1250)) +
    theme_plots() +
    ggplot2::theme(
      strip.placement = "outside",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )

}
