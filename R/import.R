# Functions for importing data

#========================================================================>
# Importing

#------------------------------------------------------------------------
#' Import csv output from Cary Eclipse Spectrophotometer
#'
#' Imports csv output from Cary Eclipse Spectrophotometer and converts it
#' into long format.
#'
#' @param filename Filename or pattern for matching multiple filenames.
#' @param debug TRUE to print out detailed messages and timing (data
#'	        import may take a few minutes).
#' @param cache TRUE to cache processed files. Processed files are
#'		recorded by a hash of their contents, ignoring file names
#'		and ensuring a new import if contents change.
#
#' @return A data.frame with columns of ...
#'
#' @examples
#' # To do
#' @export

import_data <- function(filename, debug = TRUE) {

  # Reading file
  if (debug) cat('Reading file... ')
  t <- Sys.time()
  d <- read.csv(filename, stringsAsFactors = FALSE, 
                blank.lines.skip = FALSE, header = FALSE)
  if (debug) cat(sprintf('done in %.3f min\n', as.numeric(Sys.time() - t)/60))

  # Getting number of data columns 
  # (assuming that all data columns will have an entry in the 3rd row)
  if (debug) cat('Getting number of columns... ')
  t <- Sys.time()
  n.col <- ncol(d)
  while (is.na(d[3, n.col]) || d[3, n.col] == '') {
    n.col <- n.col - 1
  }
  d <- d[, 1:n.col]
  if (debug) cat(sprintf('(%i columns) ', n.col))
  if (debug) cat(sprintf('done in %.3f min\n', as.numeric(Sys.time() - t)/60))
  
  # Dropping second row headers
  d <- d[-2, ] 

  # Getting number of data rows 
  if (debug) cat('Getting number of rows... ')
  t <- Sys.time()
  for (n.row in 1:(nrow(d) - 1)) {
    if (is.na(d[n.row + 1, 1]) || d[n.row + 1, 1] == '') break
  }
  while (! all(is.na(d[n.row + 1, ]) | d[n.row + 1, ] == '')) {
    n.row <- n.row + 1
  }
  d <- d[1:n.row, ]
  if (debug) cat(sprintf('(%i rows) ', n.row))
  if (debug) cat(sprintf('done in %.3f min\n', as.numeric(Sys.time() - t)/60))

  # Modifying column names
  if (debug) cat('Modifying column names... ')
  t <- Sys.time()
  sample.excitation <- rep(d[1, ][seq(1, n.col, by = 2)], each = 2)
  colnames(d) <- paste(sample.excitation, 
                       rep(c('wavelength', 'intensity'), each = 2),
                       sep = '_')

  # Dropping first row headers   
  d <- d[-1, ]
  print(head(d[, 1:5]))
  if (debug) cat(sprintf('done in %.3f min\n', as.numeric(Sys.time() - t)/60))

  # Gathering all columns
  if (debug) cat('Gathering columns... ')
  t <- Sys.time()
  m <- gather(d, id, value)
  if (debug) cat(sprintf('done in %.3f min\n', as.numeric(Sys.time() - t)/60))

  # Adding file name
  return(d)
}