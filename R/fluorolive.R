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
  path <- system.file('fluorolive', package = 'fluoroscripts')
  
  if (!require('shiny')) {
    message <- 'The package "shiny" needs to be installed for interactive plots'
    error(message)
  }

  runApp(path, ...)
}
