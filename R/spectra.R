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
