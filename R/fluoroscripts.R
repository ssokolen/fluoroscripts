# Documentation for provided data-sets

#' Summary information of selected fluorescent proteins.
#'
#' Information was largely taken from the data provided by Talley Lambert
#' and Kurt Thorn under the Creative Commons Attribution 3.0 license
#' (nic.ucsf.edu/FPvisualization/). Less detailed information was added
#' from Addgene's Fluorescent Protein Guide 
#' (https://www.addgene.org/fluorescent-proteins/plasmid-backbones/).
#' Small changes were made to protein names for consistency.
#'
#' @format A data frame with 95 rows and 14 variables summarizing common
#'         fluorescent protein properties. See nic.ucsf.edu/FPvisualization/
#'         for more information.
#' \describe{
#'    \item{protein}{Name of protein}
#'    \item{peak_excitation}{Wavelength (nm) of peak excitation}
#'    \item{peak_emission}{Wavelength (nm) of peak emission}
#'    \item{stokes}{Stokes shift (difference between emission and excitation)}
#'    \item{extinction}{Extinction coefficient}
#'    \item{qy}{Quantum yield}
#'    \item{brightness}{Normalized product of extinction coefficient and 
#'                      quantum yield}
#'    \item{agg}{Aggregation state of protein}
#'    \item{pka}{pKa of protein}
#'    \item{bleach}{Bleaching time in seconds}
#'    \item{mature}{Maturation time in minutes}
#'    \item{lifetime}{Fluorescence lifetime in nanoseconds}
#'    \item{doi}{Literature reference}
#'    \item{addgene}{TRUE if availabel through addgene}
#' }
"protein.summary"

#' Spectral information of selected fluorescent proteins.
#'
#' Information combines spectra available through 
#' http://www.spectra.arizona.edu/ with manual spectra plot digitation from
#' original publications.
#'
#' @format A data frame with 5600 rows and 4 variables. 
#' \describe{
#'    \item{protein}{Name of protein}
#'    \item{wavelength}{Wavelength in nm}
#'    \item{excitation}{Relative excitation curve scaled to range of 0-1}
#'    \item{emission}{Relative emission curve scaled to range of 0-1}
#' }
"protein.spectra"



