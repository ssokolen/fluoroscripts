# Functions related to interactive plotting


#========================================================================>
# Launcher

#------------------------------------------------------------------------
#' Launch "app" for interactive spectra display.
#'
#' A thin wrapper around shiny's runApp() that provides the path to
#' required server.R and ui.R files.
#
#' @param ... All arguments passed into runApp. See runApp for details.
#
#' @export
#'
#------------------------------------------------------------------------
fluorolive <- function(...) {

  if (!require('opencpu')) {
    msg <- 'The package "opencpu" needs to be installed for interactive plots'
    stop(msg)
  }

  opencpu$browse('/library/fluoroscripts/www', ...)
}

#========================================================================>
# Helper functions

#------------------------------------------------------------------------
#' Finds index of closest element in vector to specified value.
#'
#' Finds index of closest element in vector to specified value.
#'
#' @param values Vector of values for comparison.
#' @param comparison Comparison value.
#'
#' @return An integer vector.
#' @export
#'
#------------------------------------------------------------------------
closest_match <- function(values, comparison) {
  out <- which(abs(values - comparison) == min(abs(values - comparison)))
}

#========================================================================>
# Plotting functions

#------------------------------------------------------------------------
#' Plot a selection of protein spectra over laser and channel data
#' (intended for interactive plotting).
#'
#' Scales protein spectra based on laser excitation and brightness data
#' over top of selected channels.
#'
#' @param proteins Vector of protein names.
#' @param lasers Vector of laser wavelengths.
#' @param channels Vector of channels e.g. c('530/30', '670LP').
#' @param protein.colours Vector of colours as hex or characters.
#' @param laser.colours Vector of colours as hex or characters.
#' @param channel.colours Vector of colours as hex or characters.
#' @param scale.brightness A logical value indicating whether the spectra
#'                         should be scaled by brightness (not all proteins
#'                         have brightness values).
#'
#' @return An R plot object.
#' @export
#'
#------------------------------------------------------------------------
iplot_spectra <- function(proteins, lasers = NULL, channels = NULL, 
                          protein.colours = NA, laser.colours = NA, 
                          channel.colours = NA, scale.brightness = TRUE) {
  
  # Basic input formatting
  proteins <- as.character(proteins)
  lasers <- as.numeric(lasers)
  channel.ranges <- calculate_ranges(channels, 400, 800)

  # Checking if all proteins are valid
  not.valid <- !proteins %in% protein.spectra$protein
  if (any(not.valid)) {
    string <- "No spectra for the following proteins: %s" 
    msg <- sprintf(string, paste(proteins[not.valid], collapse = ', '))
    stop(msg)
  }

  # Checking whether we have all required brightness values
  if (scale.brightness == TRUE) {
    logic <- !is.na(protein.summary$brightness)
    not.valid <- !proteins %in% protein.summary$protein[logic]

    if (any(not.valid)) {
      string <- "Brightness values not available for the following proteins: %s" 
      msg <- sprintf(string, paste(proteins[not.valid], collapse = ', '))
      stop(msg)
    }
  }

  # Setting default line width
  par(lwd = 1)

  # Initializing plot
  plot(0, type = 'n', xlab = 'Emission (nm)', ylab = 'Relative intensity',
       xlim = c(400, 800), ylim = c(0, 1))

  # Adding channel blocks
  if (length(channels) > 0) {
    if (is.na(channel.colours)) channel.colours = rep(NA, length(channels))
      
    for (i in 1:length(channels)) {
      if (is.na(channel.colours[i])) {
        channel.colours[i] <- calculate_colours(mean(channel.ranges[[i]]))
      }
      
      colour.rgb <- col2rgb(channel.colours[i])
      colour <- rgb(colour.rgb[1], colour.rgb[2], colour.rgb[3], 
                    alpha = 50, maxColorValue = 255)

      rect(channel.ranges[[i]][1], -.5, channel.ranges[[i]][2], 1.5, 
           col = colour, border = FALSE)
    }
  }

  # Changing line width
  par(lwd = 3)

  # Adding laser lines
  if (length(lasers) > 0) {
    if (is.na(laser.colours)) laser.colours = rep(NA, length(lasers))
      
    for (i in 1:length(lasers)) {
      if (is.na(laser.colours[i])) {
        laser.colours[i] <- calculate_colours(lasers[i])
      }

      abline(v = lasers[i], col = laser.colours[i])
    }
  }

  # Initializing spectra subset
  d <- protein.spectra %>%
         filter(protein %in% proteins) %>%
         mutate(scaled = 0)

  # Scaling by brightnes if needed
  if (scale.brightness) {
    d <- d %>% 
           left_join(protein.summary[, c('protein', 'brightness')], 
                     by = 'protein') %>%
           mutate(emission = emission * brightness)
  }

  # Looping through and adding the excitation power of each laser on
  # each protein
  if (length(lasers) > 0) {
    for (i in 1:length(lasers)) {
      d.sub <- d %>% 
             filter(wavelength > lasers[i]) %>%
             group_by(protein) %>%
             mutate(scaled = scaled + 
                    emission*excitation[closest_match(wavelength, lasers[i])])
      d[d$wavelength > lasers[i], ] <- d.sub
    }
  }
  else {
    d$scaled <- d$emission
    # If there were no lasers, then original emission values fill scaled
  }

  d$scaled <- d$scaled/max(d$scaled)

  # Plotting lines
  if (is.na(protein.colours)) protein.colours = rep(NA, length(proteins))

  for (i in 1:length(proteins)) {
    if (is.na(protein.colours[i])) {
      peak <- filter(protein.summary, protein == proteins[i])$peak_emission
      protein.colours[i] <- calculate_colours(peak)
    }
    
    s <- filter(d, protein == proteins[i])
    lines(s$wavelength, s$scaled, col = protein.colours[i])
  }

  # Returning
  #invisible()
}
