# _targets.R

#setup
devtools::load_all()
library(targets)
library(tarchetypes)
source("R/functions.R")

conflicted::conflict_prefer("filter", "dplyr")

options(
  tidyverse.quiet = TRUE,
  usethis.quiet = TRUE
)

future::plan(future::multisession(workers = future::availableCores() - 1))

# target-specific options
tar_option_set(
  packages = c("tidyverse"),
  format = "qs"
)

# list of target objects
list(
  tar_target(dna_per_cell_file, path_to_data("dna-per-cell-number.xlsx")),
  tar_target(dna_per_cell_raw, clean_dna_per_cell(dna_per_cell_file)),
  tar_target(dna_per_cell_std, make_std_curves(dna_per_cell_raw)),
  tar_target(dna_per_cell_clean, interp_data(dna_per_cell_raw, dna_per_cell_std)),
  tar_target(cells_per_dna, calculate_cells_per_dna(dna_per_cell_clean)),
  tar_render(
    dna_per_cell_report,
    path = path_to_reports("dna-per-cell.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  )
)
