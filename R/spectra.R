# Manipulation of emission spectra for area/overlap calculation

#=========================================================================>
# Helper functions

#------------------------------------------------------------------------
#' Calculate flow cytometer channel range.
#'
#' Converts channel specifications into a list of c(min, max) ranges. Parsing
#' script inteprets specifications of the form xxxLP for long pass, xxxSP for
#' short pass, or xxx/yy for bandpass filter.
#
#' @param channels Vector of channels e.g. c('530/30', '670LP').
#' @param min.wavelengths Minimum wavlength for calculating short pass extent.
#' @param max.wavelengths Maximum wavlength for calculating long pass extent.
#
#' @return A list of c(min, max) ranges. 
#' @export
#'
calculate_ranges <- function(channels, min.wavelength, max.wavelength) {

  ranges <- list()
  for (i in 1:length(channels)) {

    if (!is.na(str_match(channels[i], '\\d{3}SP')[1])) {

      value <- as.numeric(str_extract(channels[i], '\\d+'))
      ranges[[i]] <- c(min.wavelength, value)

    } else if (!is.na(str_match(channels[i], '\\d{3}LP')[1])) {

      value <- as.numeric(str_extract(channels[i], '\\d+'))
      ranges[[i]] <- c(value, max.wavelength)

    } else if (!is.na(str_match(channels[i], '\\d{3}/\\d{2}')[1])) {

      value <- as.numeric(str_extract(channels[i], '\\d+(?=/)'))
      range <- as.numeric(str_extract(channels[i], '(?<=/)\\d+'))
      ranges[[i]] <- c(value - range/2, value + range/2)

    } else {

      msg <- 'Filter format must be one of "xxxSP", "xxxLP", or "xxx/yy".'
      stop(msg)

    }
  }

  return(ranges)
}

#=========================================================================>
# Functions relating areas under spectral curves

#------------------------------------------------------------------------
#' Calculate area of spectral curve
#'
#' Calculates the area under an emission spectra within specified channels,
#' excited by given laser wavelength. The output is a dataframe to
#' facilitate use within dplyr's do() function.
#
#' @param wavelengths A vector of wavelengths (nm).
#' @param excitations A vector of excitation intensities.
#' @param emissions A vector of emmision intensities.
#' @param lasers Vector of lasers exciting given proteins e.g. c(488, 642).
#' @param channels Vector of channels for which spectral area is calculated
#'                 e.g. c('530/30', '670LP').
#' @param threshold Fraction of total spectral area observed in a given
#'                  channel that is considered to be practically 0.
#
#' @return A data frame with two columns -- channel and area.
#' @export
#'
#------------------------------------------------------------------------
# Calculating channel area
calculate_area <- function(wavelengths, excitations, emissions, lasers,
                           channels, threshold = 1e-2) {

  # Generating spectral function
  f_excitation <- splinefun(wavelengths, excitations)
  f_emission <- splinefun(wavelengths, emissions)

  scaling.factor <- f_excitation(lasers)/max(excitations)

  # Initializing output
  out <- data.frame(channel = channels, area = 0)

  # Integrating
  w.min <- min(wavelengths)
  w.max <- max(wavelengths)
  ranges <- calculate_ranges(channels, w.min, w.max)

  # A laser functions as a lower limit on emission wavelength
  total.area <- 0
  for (j in 1:length(lasers)) {
    low <- max(w.min, lasers[j])
    high <- max(w.max, lasers[j])

    if (high > low) {
      curve.area <- integrate(f_emission, w.min, w.max)$value
      total.area <- total.area + curve.area * scaling.factor[j]
    }
  }

  if (total.area == 0) return(out)

  # Looping through channels and lasers to increment area 
  for (i in 1:length(channels)) {
    
    for (j in 1:length(lasers)) {
      low <- max(ranges[[i]][1], lasers[j])
      high <- max(ranges[[i]][2], lasers[j])

      if (high > low) {
        curve.area <- integrate(f_emission, low, high)$value 
        out[i, 'area'] <- out[i, 'area'] + curve.area * scaling.factor[j]
      }
    }
  }
  
  # Filtering out values below the threshold 
  out <- mutate(out, area = ifelse(area/total.area < threshold, 0, area))

  return(out)
}
