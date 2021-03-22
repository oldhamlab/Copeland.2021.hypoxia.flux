#' Corrected mass isotopomer distributions
#'
#' A data set containing corrected mass isotope distributions. Technical
#' replicates have been filtered from the data.
#'
#' \describe{
#'   \item{method}{
#'   `fs` = full scan \cr
#'   `sim` = selected ion monitoring}
#'   \item{cell_type}{
#'   `lf` = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{date}{start date of an experiment}
#'   \item{tracer}{
#'   `glc2` = \[1,2-13C2\] glucose \cr
#'   `glc6` = \[U-13C6\] glucose \cr
#'   `q5` = \[U-13C5\] glutamine \cr
#'   `lac3` = \[U-13C3\] lactate}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `none` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY` = 10 μM molidustat}
#'   \item{time}{hours}
#'   \item{metabolite}{name of measurement}
#'   \item{isotope}{mass shift of detected isotope}
#'   \item{mid}{mole fraction of the indicated isotope}
#'   }
"mids"
