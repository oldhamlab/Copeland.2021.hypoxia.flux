# isotope correction ------------------------------------------------------

isotope_library <-
  tibble::tribble(
    ~ metabolite, ~ formula, ~ polarity,
    "2HG", "C5H8O5", "negative",
    "2OG", "C5H6O5", "negative",
    "alanine", "C3H7NO2", "negative",
    "aspartate", "C4H7NO4", "negative",
    "citrate", "C6H8O7", "negative",
    "glutamate", "C5H9NO4", "negative",
    "glutamine", "C5H10N2O3", "negative",
    "lactate", "C3H6O3", "negative",
    "malate", "C4H6O5", "negative",
    "pyruvate", "C3H4O3", "negative",
    "serine", "C3H7NO3", "negative",
    "succinate", "C4H6O4", "negative",
    "3PG", "C3H7O7P", "negative",
    "aconitate", "C6H6O6", "negative",
    "FBP", "C6H14O12P2", "negative",
    "G3P", "C3H9O6P", "negative",
    "palmitate", "C16H32O2", "negative",
    "PEP", "C3H5O6P", "negative",
    "sedoheptulose", "C7H14O7", "negative",
    "DHAP", "C3H7O6P", "negative",
    "GAP", "C3H7O6P", "negative",
    "G1P", "C6H13O9P", "negative",
    "G6P", "C6H13O9P", "negative",
    "R5P", "C5H11O8P", "negative"
  )


# biomass calculation -----------------------------------------------------

cell_composition <-
  tibble::tribble(
    ~ metabolite, ~ μmol_per_gDW,
    "ALA", 0.6,
    "ARG", 0.377,
    "ASP", 0.359,
    "ASN", 0.288,
    "CYS", 0.145,
    "GLN", 0.322,
    "GLU", 0.386,
    "GLY", 0.538,
    "HIS", 0.143,
    "ILE", 0.324,
    "LEU", 0.564,
    "LYS", 0.57,
    "MET", 0.138,
    "PHE", 0.219,
    "PRO", 0.313,
    "SER", 0.43,
    "THR", 0.386,
    "TRP", 0.044,
    "TYR", 0.182,
    "VAL", 0.416,
    "glycogen", 0.279,
    "dAMP", 0.0148,
    "dGMP", 0.0099,
    "dCMP", 0.0099,
    "dTMP", 0.0148,
    "AMP", 0.033,
    "GMP", 0.0624,
    "CMP", 0.0551,
    "UMP", 0.033,
    "cholesterol", 0.018,
    "phosphatidylcholine", 0.069,
    "phosphatidylethanolamine", 0.026,
    "phosphatidylinositol", 0.01,
    "phosphatidylserine", 0.003,
    "phosphatidylglycerol", 0.001,
    "sphingomyelin", 0.008,
    "cardiolipin", 0.003
  ) %>%
  dplyr::mutate(class = case_when(
    .data$metabolite %in% c("dAMP", "dGMP", "AMP", "GMP") ~ "purine",
    .data$metabolite %in% c("dCMP", "CMP", "UMP") ~ "pyrimidine",
    .data$metabolite == "dTMP" ~ "thymine",
    TRUE ~ metabolite
  ))

metabolite_composition <-
  tibble::tribble(
    ~ class, ~ component, ~ stoichiometry,
    "purine", c("P5P", "GLY", "CO2", "MEETHF"), c(1, 1, 1, 2),
    "pyrimidine", c("P5P", "CO2", "ASP"), c(1, 1, 1),
    "thymine", c("P5P", "CO2", "ASP", "MEETHF"), c(1, 1, 1, 1),
    "glycogen", "G6P", 1,
    "cholesterol", "AcCoA", 18,
    "phosphatidylcholine", c("AcCoA", "DHAP"), c(17.43, 1),
    "phosphatidylethanolamine", c("AcCoA", "DHAP"), c(17.43, 1),
    "phosphatidylserine", c("AcCoA", "DHAP", "SER"), c(17.43, 1, 1),
    "phosphatidylinositol", c("AcCoA", "DHAP", "G6P"), c(17.43, 1, 1),
    "phosphatidylglycerol", c("AcCoA", "DHAP"), c(17.43, 2),
    "sphingomyelin", c("AcCoA", "SER"), c(17.43, 1),
    "cardiolipin", c("AcCoA", "DHAP"), c(34.86, 2)
  ) %>%
  tidyr::unnest(c(component, stoichiometry))

μmol_per_mass <-
  dplyr::left_join(cell_composition, metabolite_composition, by = "class") %>%
  dplyr::mutate(
    component = dplyr::if_else(is.na(.data$component), .data$metabolite, .data$component),
    stoichiometry = tidyr::replace_na(stoichiometry, 1),
    μmol_per_g = .data$μmol_per_gDW * .data$stoichiometry
  ) %>%
  dplyr::group_by(.data$component) %>%
  dplyr::summarise(μmol_per_gDW = sum(.data$μmol_per_g)) %>%
  dplyr::rename(metabolite = .data$component)


# usethis -----------------------------------------------------------------

usethis::use_data(
  isotope_library,
  μmol_per_mass,
  internal = TRUE,
  overwrite = TRUE
)
