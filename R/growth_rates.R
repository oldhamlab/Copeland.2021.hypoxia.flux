#' Growth rate parameters
#'
#' A data set containing the growth rate and X0 values from linear fitting of
#' cell count data.
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
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY` = 10 μM molidustat}
#'   \item{X0}{cell count at time 0}
#'   \item{mu}{growth rate per hour}
#'   }
"growth_rates"
