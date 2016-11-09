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
iplot_spectra <- function(proteins, lasers, channels, protein.colours = NA, 
                          laser.colours = NA, channel.colours = NA,
                          scale.brightness = FALSE) {

  # Basic input formatting
  proteins <- as.character(proteins)
  lasers <- as.numeric(lasers)
  channel.ranges <- calculate_ranges(channels, 500, 800)

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

  # Initializing plot
  plot(0, type = 'n', xlab = 'Emission (nm)', ylab = 'Relative intensity',
       xlim = c(500, 800), ylim = c(0, 1))

  # Adding channel blocks

  # Returning
  invisible()
}
