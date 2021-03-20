#' Calculated fluxes
#'
#' A data set containing metabolite fluxes.
#'
#'  \describe{
#'   \item{metabolite}{name of measurement}
#'   \item{abbreviation}{abbreviated metabolite name}
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
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY` = 10 μM molidustat}
#'   \item{flux}{fmol / cell / h, positive fluxes indicate secretion, negative fluxes
#'   indicate uptake}
#'   }
"fluxes"
