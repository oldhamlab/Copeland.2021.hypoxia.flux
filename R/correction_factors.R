#' Quadrupole bias correction factors
#'
#' A data set containing experimentally determined correction factors to adjust
#' peak areas for isotopes detected in selected ion monitoring mode. Isotope peak
#' areas are multiplied by the correction factor. The corrected peak area is
#' subsequently used to determine the mass isotope distribution.
#'
#'  \describe{
#'    \item{batch}{correction factors associated with specific experiments}
#'    \item{metabolite}{name of measurement}
#'    \item{M}{isotope}
#'    \item{cf}{correction factor}
#'   }
"correction_factors"
