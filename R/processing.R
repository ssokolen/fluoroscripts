# Functions for processing 

#========================================================================>
# Cleaning up acquired spectra

#------------------------------------------------------------------------
#' Identify emission scatter.
#'
#' Identifies sharp spikes in spectra that are composed of relatively few
#' points. Intensity data is filtered with a 3rd order Savitzky–Golay
#' filter using 7 and 21 points. If the first or second derivative
#' at an observed point is 5 times greater when using fewer points, it is
#' flagged as unusually sharp. A continuous cluster of points that
#' occurs within a specified wavelength of the excitation is identified as
#' a spike. If more than one unique excitation is provided, the median spike
#' width is used across all emissions. All integer multiples of the
#' excitation wavelength are also assessed.
#
#' @param excitations Vector of excitation wavelengths or a single excitation
#'                   value.
#' @param emissions Vector of emission wavelengths.
#' @param intensities Vector of emission intensities.
#' @param order Polynomial order to use in Savitzky-Golay filter.
#' @param n1 Number of points to use in first Savitzky–Golay.
#' @param n2 Number of points to use in second Savitzky–Golay.
#' @param r1 Threshold ratio between first derivative of filtered values
#'           (first filter divided by second).
#' @param r2 Threshold ratio between second derivative of filtered values
#'           (first filter divided by second).
#' @param m1 Minimum value of first derivative calculated with the first
#'           filter.
#' @param m2 Minimum value of second derivative calculated with the first
#'           filter.
#
#' @return A vector of booleans corresponding to each intensity where TRUE
#'         indicates a scatter spike.
#' @export
#'
#------------------------------------------------------------------------
flag_scatter <- function(excitations, emissions, intensities, order = 3, 
                         n1 = 7, n2 = 21, r1 = 5, r2 = 5, m1 = 50, m2 = 20) {

  # Generating data frame
  d <- data.frame(excitation = excitations, emission = emissions,
                  intensity = intensities)

  # Defining filters
  filt.n1d1 <- signal::sgolay(order, n1, m = 1)
  filt.n1d2 <- signal::sgolay(order, n1, m = 2)

  filt.n2d1 <- signal::sgolay(order, n2, m = 1)
  filt.n2d2 <- signal::sgolay(order, n2, m = 2)

  # Defining filter function (outputs number of points to filter)
  f_width <- function(excitation, emission, intensity) {

    # Replacing overflow by one point (to preserve peak sharpness)
    # It's safe to assume that all maximum values are overflow
    max.value <- max(intensity)
    filtered <- intensity
    n <- length(intensity)
    indeces <- 1:n

    temp.indeces <- c()
    overflow.indeces <- c()

    for (i in 1:(n - 1)) {
      if (filtered[i] == max.value) {
        temp.indeces <- c(temp.indeces, i)
      } else {
        n.temp <- length(temp.indeces)
        if (n.temp > 0) {
          if (n.temp > 1) {
            overflow.indeces <- c(overflow.indeces, 
                                  temp.indeces[1:(n.temp - 1)])
          }
          temp.indeces <- c()
        }
      }
    }

    if (length(overflow.indeces) > 0) { 
      filtered <- filtered[-overflow.indeces]
      indeces <- indeces[-overflow.indeces]
    }

    # Filtering values
    y.n1d1 <- signal::filter(filt.n1d1, filtered)
    y.n1d2 <- signal::filter(filt.n1d2, filtered)
    y.n2d1 <- signal::filter(filt.n2d1, filtered)
    y.n2d2 <- signal::filter(filt.n2d2, filtered)

    logic1 <- (abs(y.n1d1 / y.n2d1) > r1) & (y.n1d1 > m1)
    logic2 <- (abs(y.n1d2 / y.n2d2) > r2) & (y.n1d2 > m2)

    # Setting all overflow values and all filtered values to NA
    intensity[intensity == max.value] <- NA
    intensity[indeces[logic1 | logic2]] <- NA

    # Counting NA values from center of scatter
    center <- which.min(abs(emission - excitation))

    # Right points
    i.right <- center
    n.right <- 0
    while ((i.right <= n) && (is.na(intensity[i.right]))) {
      i.right <- i.right + 1
      n.right <- n.right + 1
    }

    # Left points
    i.left <- center
    n.left <- 0
    while ((i.left >= 1) && (is.na(intensity[i.left]))) {
      i.left <- i.left - 1
      n.left <- n.left + 1
    }

    # Forcing symmetry
    n.left <- max(n.right, n.left)
    n.right <- max(n.right, n.left)

    i.left <- max(center - n.left, 1)
    i.right <- min(center + n.right, n)

    # Returning half-width
    return(abs(emission[i.right] - emission[i.left])/2)
  }

  # Running f_width to determine width of scatter
  d.width <- d %>%
               group_by(excitation) %>%
               summarize(width = f_width(excitation, emission, intensity))

  # Taking median 
  width <- median(median(d.width$width))

  # Filtering out integer multiples of excitation (using same scatter width)
  multiples <- 1:(max(d$emission)/min(d$excitation))
  logic <- rep(FALSE, length(d$intensity))
  for (i in multiples) {
    logic <- logic | (abs(d$excitation*i - d$emission) < width)
  }

  return(logic)
}
