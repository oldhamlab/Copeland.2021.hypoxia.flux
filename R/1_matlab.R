# matlab.R

clean_biomass <- function(biomass_file) {
  readr::read_csv(biomass_file) %>%
    dplyr::mutate(
      diff = post-pre,
      cell_mass = diff / cell_number * 1E12
    )
}

calculate_biomass <- function(biomass_clean) {
  biomass_clean %>%
    dplyr::group_by(.data$cell_type, .data$date) %>%
    dplyr::mutate(cell_mass = replace_outliers(.data$cell_mass)) %>%
    dplyr::summarise(cell_mass = mean(.data$cell_mass, na.rm = TRUE)) %>%
    dplyr::mutate(cell_mass = replace_outliers(.data$cell_mass)) %>%
    dplyr::summarise(biomass = mean(.data$cell_mass, na.rm = TRUE))
}

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
    dplyr::mutate(
      metabolite = replace(.data$metabolite, .data$metabolite == "glutamine", "GLN"),
      metabolite = replace(.data$metabolite, .data$metabolite == "2OG", "AKG"),
      metabolite = stringr::str_sub(.data$metabolite, 1, 3),
      metabolite = toupper(.data$metabolite)
    ) %>%
    dplyr::select(-.data$method)
}

summarize_mids <- function(df) {
  df %>%
    dplyr::summarise(
      mean = mean(.data$mid, na.rm = TRUE),
      se = sd(.data$mid, na.rm = TRUE)/sqrt(dplyr::n())
    ) %>%
    dplyr::ungroup() %>%
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
