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
#' Calculate area under spectral curve
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

#------------------------------------------------------------------------
#' Calculate overlap between proteins
#'
#' Calculates the overlaps in area under an emission spectra between all
#' proteins. Areas are calculated by taking into account the specified 
#' channels, and excitation by given laser wavelength. Overlap is reported
#' as a fraction of area of specified protein. As an example the overlap
#' of EGFP by EYFP will likely be different from the overlap of EYFP by EGFP.
#' An overlap of NA means that the area of specified protein is 0 for the
#' given channel (and would result in division by 0). An overlap of 0 means
#' that the relative area of overlapping proteins is practically 0.
#
#' @param wavelengths A vector of wavelengths (nm).
#' @param excitations A vector of excitation intensities.
#' @param emissions A vector of emmision intensities.
#' @param proteins A vector of protein labels that correspond to each
#'                 measurement.
#' @param lasers Vector of lasers exciting given proteins e.g. c(488, 642).
#' @param channels Vector of channels for which spectral area is calculated
#'                 e.g. c('530/30', '670LP').
#' @param threshold Fraction of total spectral area observed in a given
#'                  channel that is considered to be practically 0.
#
#' @return A data frame with four columns -- protein1, protein2, channel 
#'         and overlap. Overlap represents the fraction of protein1 area
#'         occupied by protein2 in a given channel.
#' @export
#'
#------------------------------------------------------------------------
calculate_pairwise_overlap <- function(wavelengths, excitations, emissions, 
                                       proteins, lasers, channels, 
                                       threshold = 1e-2) {

  # Generating data frame
  d <- data.frame(protein = proteins, wavelength = wavelengths,
                  excitation = excitations, emission = emissions)

  # Ordering by wavelength
  d.max <- d %>%
             group_by(protein) %>%
             summarize(wavelength = wavelength[emission == max(emission)]) %>%
             arrange(wavelength)

  d$protein <- factor(as.character(d$protein), 
                    levels = as.character(d.max$protein))
  d <- arrange(d, protein)

  # Calculating areas
  areas <- d %>%
             group_by(protein) %>%
             do(calculate_area(.$wavelength, .$excitation, .$emission,
                               lasers = lasers, channels = channels, 
                               threshold = threshold))

  # Generating a Cartesian product
  areas$id <- 1
  areas2 <- rename(areas, protein2 = protein, area2 = area)
 
  # Calculating overlap 
  combined <- left_join(areas, areas2, by = c('id', 'channel')) %>%
                select(protein1 = protein, protein2, 
                       channel, area1 = area, area2) %>%
                filter(protein1 != protein2) %>%
                mutate(overlap = ifelse(area1 > 0, area2/area1, NA)) %>%
                select(-area1, -area2)

  return(combined)
}

#------------------------------------------------------------------------
#' Assess protein combinations
#'
#' Generates all possible combinations of proteins taken n at time, chooses
#' optimal detection channels by selecting those with least overlap and
#' calculates three metrics -- max total overlap for the selected channels,
#' sum of all spectral areas for the selected channels, and the ratio
#' of largest calculated area to smallest. Areas are calculated by taking into 
#' account excitation by given laser wavelength. As taking n combinations
#' may result in large computation time, it is recommended that some
#' prescreening is performed before running the assessment.
#
#' @param wavelengths A vector of wavelengths (nm).
#' @param excitations A vector of excitation intensities.
#' @param emissions A vector of emmision intensities.
#' @param proteins A vector of protein labels that correspond to each
#'                 measurement.
#' @param n The number of protein combinations to consider.
#' @param lasers Vector of lasers exciting given proteins e.g. c(488, 642).
#' @param channels Vector of channels for which spectral area is calculated
#'                 e.g. c('530/30', '670LP').
#' @param spectra.threshold Fraction of total spectral area observed in a 
#'                          given channel that is considered to be 
#'                          practically 0.
#' @param channel.threshold Fraction of maximum area observed within a given
#'                          channel that is considered to be practically 0
#'                          (compared to all channel areas of the protein).
#' @param intensity.threshold Fraction of maximum area observed across all 
#'                            channels that is considered to be practically 0
#'                            (compared to selected protein subset).
#
#' @return A data frame listing the combination of protein1, protein2,
#'         protein3 individually, combined combination string, maximum
#'         overlap, sum of detection areas, ratio of largest area to
#'         smallest, and the detection strategy string that combines protein,
#'         channel, and channel overlap information.
#'         
#' @export
#'
#------------------------------------------------------------------------
assess_combinations <- function(wavelengths, excitations, emissions, 
                                proteins, n, lasers, channels, 
                                spectra.threshold = 1e-2,
                                channel.threshold = 1e-2,
                                intensity.threshold = 1e-2) {

  # Checking input
  if (n > length(channels)) {
    msg <- 'Number of proteins (n) must not exceed number of channels.'
    stop(msg)
  }

  # Checking input
  if (n <= 1) {
    msg <- 'Number of proteins (n) must be 2 or greater.'
    stop(msg)
  }

  # Generating data frame
  d <- data.frame(protein = proteins, wavelength = wavelengths,
                  excitation = excitations, emission = emissions)

  # Ordering by wavelength
  d.max <- d %>%
             group_by(protein) %>%
             summarize(wavelength = wavelength[emission == max(emission)]) %>%
             arrange(wavelength)

  d$protein <- factor(as.character(d$protein), 
                    levels = as.character(d.max$protein))
  d <- arrange(d, protein)

  # Calculating areas
  areas <- d %>%
             group_by(protein) %>%
             do(calculate_area(.$wavelength, .$excitation, .$emission,
                               lasers = lasers, channels = channels, 
                               threshold = spectra.threshold))

  # Generating a Cartesian product
  areas$id <- 1
  areas2 <- rename(areas, protein2 = protein, area2 = area)
 
  # Calculating overlap 
  overlap <- left_join(areas, areas2, by = c('id', 'channel')) %>%
               select(protein1 = protein, protein2, 
                      channel, area1 = area, area2) %>%
               filter(protein1 != protein2) %>%
               ungroup() %>%
               mutate(overlap = ifelse(area1 > 0, area2/area1, NA)) %>%
               select(-area1, -area2)

  # Combinations
  out <- t(combn(unique(as.character(areas$protein)), n))
  out <- as.data.frame(out)
  colnames(out) <- paste('protein', 1:n, sep = '')
  
  proteins <- paste('protein', 1:n, sep = '')
  out <- cbind(out, unite_(out, 'combination', proteins, sep = '-'))
  out$overlap.max <- NA
  out$area.sum <- NA
  out$area.diff <- NA
  out$detection <- NA

  for (col in proteins) {
    out[, col] <- as.character(out[, col])
  }

  # Looping through each combination
  for (i in 1:nrow(out)) {
    combination <- as.character(out[i, proteins])
    s.areas <- filter(areas, protein %in% combination)

    # Filtering proteins by relative intensity across all channels
    area.sums <- s.areas %>%
                   group_by(protein) %>%
                   summarize(area = sum(area))
    max.area <- max(area.sums$area)

    area.sums <- filter(area.sums, area > intensity.threshold * max.area)

    # Skip iteration if any proteins are dropped
    if (length(area.sums$protein) < n) next

    # Subsetting overlaps
    s.overlap <- filter(overlap, protein1 %in% combination, 
                                 protein2 %in% combination)
    
    # Calculating overlap sums 
    overlap.sums <- s.overlap %>%
                      group_by(protein1, channel) %>%
                      summarize(overlap = sum(overlap)) %>%
                      ungroup() %>%
                      arrange(overlap) %>%
                      mutate(channel = as.character(channel),
                             protein1 = as.character(protein1)) %>%
                      rename(protein = protein1)

    # Calculating fraction of area in each channel
    s.areas <- s.areas %>%
                 group_by(protein) %>%
                 mutate(fraction = area/sum(area)) %>%
                 ungroup() %>%
                 mutate(protein = as.character(protein),
                        channel = as.character(channel))


    # Merging and filtering
    columns <- c('protein', 'channel', 'area', 'fraction')
    overlap.sums <- overlap.sums %>%
                      left_join(s.areas[ , columns],
                                by = c('protein', 'channel')) %>%
                      filter(fraction > channel.threshold)

    overlaps <- c()
    chosen.areas <- c()
    chosen.proteins <- c()
    chosen.channels <- c()

    for (j in 1:n) {
      overlaps <- c(overlaps, overlap.sums$overlap[1])
      chosen.areas <- c(chosen.areas, overlap.sums$area[1])
      chosen.proteins <- c(chosen.proteins, overlap.sums$protein[1])
      chosen.channels <- c(chosen.channels, overlap.sums$channel[1])
      
      overlap.sums <- filter(overlap.sums, 
                             !channel %in% chosen.channels, 
                             !protein %in% chosen.proteins) 

      if (nrow(overlap.sums) == 0) break
    }

    if (length(overlaps) == n) {
      out$overlap.max[i] <- max(overlaps)
      out$area.diff[i] <- max(chosen.areas)/min(chosen.areas)
      out$area.sum[i] <- sum(chosen.areas)
      out$detection[i] <- paste(chosen.proteins, chosen.channels, 
                                round(overlaps, 2),
                                sep = '-', collapse = ', ')
    }
  }

  return(out)
}
