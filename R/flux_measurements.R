#' Measurements for extracellular flux calculations
#'
#' A data set containing the interpolated metabolite concentrations, cell counts,
#' and evaporation volumes for extracellular flux determinations.
#'
#' \describe{
#'   \item{cell_type}{
#'   `lf` = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{experiment}{
#'   `02` = 0.2% oxygen for hypoxia \cr
#'   `05` = 0.5% oxygen for hypoxia \cr
#'   `bay` = molidustat treatment \cr
#'   `05-bay` = 0.5% oxygen plus molidustat}
#'   \item{batch}{a group of biological replicates analyzed similarly}
#'   \item{date}{start date of an experiment}
#'   \item{metabolite}{name of measurement}
#'   \item{abbreviation}{abbreviated metabolite name}
#'   \item{detector}{
#'   `fld` = HPLC fluorescence detector \cr
#'   `mwd` = HPLC multi-wavelength detector \cr
#'   `enzyme` = enzymatic assay \cr
#'   `hplc` = OPD-derivatized HPLC detection \cr
#'   `lcms` = liquid chromatography-mass spectrometry assay \cr
#'   `picogreen` = fluorescence dye labeling of DNA}
#'   \item{type}{
#'   `cells` = conditioned medium \cr
#'   `empty` = unconditioned medium}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY` = 10 μM molidustat}
#'   \item{time}{hours}
#'   \item{well}{denotes technical replicates in each experiment}
#'   \item{conc}{measured in number for cells and μM for metabolites}
#'   \item{volume}{measured in mL, extrapolated from evaporation controls}
#'   \item{nmol}{metabolite mass}
#'   }
"flux_measurements"
