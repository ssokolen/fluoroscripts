# Functions for converting raw protein data stored in csv files into Rdata
# files loaded by the package

#========================================================================
# Spectra

#------------------------------------------------------------------------
#' Import fluorescent protein spectra
#'
#' Imports fluorescent protein data, interpolates the available points,
#' and generates new values using a common basis.
#'
#' @param filename Filename or pattern for matching multiple filenames
#'                 containing csv spectra data. If no filename is
#'                 provided, the current working directory is assumed
#'                 to be the root directory of the package source and
#'                 all csv files are matched in ./inst/protein_spectra.
#' @param destination Filename to save the generated Rdata file. If no
#'                    destination is provided, the current working
#'                    directory is assumed to be the root directory of
#'                    the package source and a protein_spectra.Rdata 
#'                    file is saved to ./data.
#' @param basis A sequence of wavelength values (with the same units
#'              as the data in the csv files) to evaluate all spectra
#'              on.
#' @param overwrite TRUE to overwirte existing Rdata file.
#'
#' @examples
#' # To do
#' @export

generate_spectra <- function(filename = NA, destination = NA, 
                             basis = seq(300, 800, length.out = 200),
                             overwrite = TRUE) {

  # Checking filename
  if (is.na(filename)) {
    msg <- paste('No filename provided, looking for csv files', 
                 'in ./inst/protein_spectra')

    warning(msg)

    filename = './inst/protein_spectra/*.csv'
  }

  spectra.files <- Sys.glob(filename)

  if (length(spectra.files) == 0) {
    msg <- sprintf('No files matching %s', filename)
    stop(msg)
  }

  # Reading in spectra files
  protein.spectra <- list()

  for (file in spectra.files) {
    d <- read.csv(file)

    d.list <- list(d[ , c(1, 2)], d[ , c(3, 4)])

    for (i in c(1, 2)) {

      d <- d.list[[i]]
      colnames(d) <- c('wavelength', 'intensity')

      d <- filter(d, !is.na(wavelength))

      f_interpolate <- splinefun(d$wavelength, d$intensity)
      f_extrapolate <- function(x) {
        low <- min(d$wavelength)
        high <- max(d$wavelength)
        protein.spectra <- ifelse((x > low) & (x < high), f_interpolate(x), 0)
        return(protein.spectra)
      }

      y <- f_extrapolate(basis)
      y[y < 0] <- 0

      # Ensuring a scale of to 0-1
      y <- (y - min(y))/(max(y) - min(y))

      d.list[[i]] <- data.frame(wavelength = basis, intensity = y)
    }

    d.list[[1]]$spectra <- 'excitation'
    d.list[[2]]$spectra <- 'emission'

    d <- rbind(d.list[[1]], d.list[[2]]) %>%
         spread(spectra, intensity) %>%
         select(wavelength, excitation, emission)

    # Trimming everything after the last period
    protein <- gsub('\\.[^.]*$', '', basename(file))
    protein.spectra[[protein]] <- d 
  }

  protein.spectra <- melt(protein.spectra, 
                          id = c('wavelength', 'excitation', 'emission')) %>%
                     select(protein = L1, everything())

  # Checking destination
  if (is.na(destination)) {
    msg <- paste('No destination provided, storing spectra as', 
                 'protein_spectra.RData in ./data')

    warning(msg)

    destination = './data/protein_spectra.RData'
  }

  if (!overwrite && file.exists(destination)) {
    msg <- sprintf('Destination file %s exists, aborting.', destination)
    stop(msg)
  }

  save(protein.spectra, file = destination)
}
