# Functions for importing data

#========================================================================>
# Importing

#------------------------------------------------------------------------
#' Import sample description metadata from csv output from 
#' Cary Eclipse Spectrophotometer
#'
#' Imports sample description metadata from csv output from a Cary Eclipse 
#' Spectrophotometer and converts it into long format.
#'
#' @param filename Filename or pattern for matching multiple filenames.
#' @param aggregate TRUE if metadata should be reported per sample
#'                  (aggregating across the excitation sweep) or per single
#'                  scan.
#' @param cache TRUE to cache processed files. Processed files are
#'		recorded by a hash of their contents, ignoring file names
#'		and ensuring a new import if contents change. [NOT IMPLEMENTED]
#
#' @return A data.frame with columns of scan parameters.
#'
#' @examples
#' # To do
#' @export

import_metadata <- function(filename, aggregate = FALSE, cache = FALSE) {

  # Parsing filename pattern
  files <- Sys.glob(filename)
  if (length(files) == 0) {
    msg = sprintf('No files matching "%s"', filename)
    stop(msg)
  }

  # Initializing output
  options(stringsAsFactors = FALSE)
  out <- data.frame()

  # Looping through each file
  for (file in files) {

    # Reading file as string
    d.string <- readLines(file)
    
    # Finding split between numeric and metadata 
    for (i in 1:length(d.string)) {
      if (substring(d.string[i], 1, 1) == '') break
    }
    split <- i

    # Checking whether data comes from right source
    if (!str_detect(d.string[split + 7], 'Cary Eclipse')) {
     msg <- 'Unexpected data format, output may be corrupted'
     warning(msg)
    }

    # Initializing output data.frame
    n <- length(d.string)
    scan <- 1
    skip <- 32
    skip.extra <- 0
    d <- data.frame()
    dates <- c()

    # Finding rest of breaks
    while (TRUE) {

      # Extracting data
      if (substring(d.string[split], 1, 1) == '') {

        label <- str_trim(str_extract(d.string[split + 2], '.*?(?=_EX_)'))
        date <- str_extract(d.string[split + 3], '\\d.*')
        time <- dmy_hms(date)

        excitation <- str_trim(str_extract(d.string[split + 2], '(?<=_EX_).*'))

        param <- rep(NA, 22)
        for (i in 1:22) {
          param[i] <- d.string[split + i + 6] %>%
                      str_extract('(?<=\\s{3})\\w.*') %>%
                      str_trim()
        }

        # Replacing parameter excitation using column label
        param[8] <- as.numeric(excitation)

        d <- rbind(d, c(list(scan, label, date, time), param))

        scan <- scan + 1
        skip <- skip + skip.extra
        skip.extra <- 0
        split <- split + skip
      } else {
        skip.extra <- skip.extra + 1
        split <- split + 1
        if (split >= n) break
      }
    }

    # Initializing output data.frame
    all.columns <- c('scan', 'label', 'date', 'time', 'instrument', 'serial', 
                     'data.mode', 'scan.mode', 'x.mode', 'start', 'stop', 
                     'excitation.wavelength', 'excitation.slit', 
                     'emission.slit', 'scan.rate', 'data.interval', 
                     'averaging.time', 'excitation.filter',
                     'emission.filter', 'pmt.voltage', 'corrected.spectra',
                     'd.mode', 'excitation.stop', 'excitation.increment', 
                     'multicell.holder', 'multi.zero')

    colnames(d) <- all.columns

    numeric.columns <- c('scan', 'time', 'start', 'stop', 
                         'excitation.wavelength', 'excitation.slit', 
                         'emission.slit', 'scan.rate', 'data.interval', 
                         'averaging.time', 'excitation.stop',
                         'excitation.increment')

    for (column in numeric.columns) {
      d[, column] <- as.numeric(d[, column])
    }

    out <- rbind(out, d)
  }

  # Adding sample and cell numbers
  out <- arrange(out, time)

  row <- 1
  cells <- c()
  samples <- c()
  min.sample <- 0
  
  while (TRUE) {
    if (row >= nrow(out)) break

    n.cells <- 0
    while (TRUE) {
      current <- out$excitation.wavelength[row + n.cells] 
      start <- out$start[row + n.cells]

      if (current == start) {
        n.cells <- n.cells + 1
      } else {
        break
      }
    }
    n.scans <- (out$stop[row] - out$start[row])/
               out$excitation.increment[row] + 1
    
    cells <- c(cells, rep(1:n.cells, n.scans))
    samples <- c(samples, min.sample +  rep(1:n.cells, n.scans))

    min.sample <- min.sample + n.cells

    row <- row + n.cells*n.scans
  }

  out <- cbind(sample = samples, cell = cells, out)

  if (aggregate) {
    agg.columns <- c('cell', 'label', 'date', 'time', 'instrument', 'serial', 
                     'data.mode', 'scan.mode', 'x.mode', 'start', 'stop', 
                     'excitation.slit', 'emission.slit', 'scan.rate', 
                     'data.interval', 'averaging.time', 'excitation.filter',
                     'emission.filter', 'pmt.voltage', 'corrected.spectra',
                     'd.mode', 'excitation.stop', 'excitation.increment', 
                     'multicell.holder', 'multi.zero')

    f_aggregate <- function(d) {
      return(d[1, agg.columns])
    }

    out <- group_by(out, sample) %>% do(f_aggregate(.))
  }

  return(out)
}

#------------------------------------------------------------------------
#' Import numeric data from csv output from Cary Eclipse Spectrophotometer
#'
#' Imports numeric csv output from a Cary Eclipse Spectrophotometer and 
#' converts it into long format.
#'
#' @param filename Filename or pattern for matching multiple filenames.
#' @param cache TRUE to cache processed files. Processed files are
#'		recorded by a hash of their contents, ignoring file names
#'		and ensuring a new import if contents change. [NOT IMPLEMENTED]
#
#' @return A data.frame with the following columns: sample, cell, scan,
#'         label, excitation, emission, intensity
#'
#' @examples
#' # To do
#' @export

import_data <- function(filename, cache = FALSE) {

  # Parsing filename pattern
  files <- Sys.glob(filename)
  if (length(files) == 0) {
    msg = sprintf('No files matching "%s"', filename)
    stop(msg)
  }

  # Getting metadata first
  meta <- import_metadata(filename)

  # Initializing output
  out <- data.frame()

  # Looping through each file
  for (file in files) {

    # Reading file as string
    d.string <- readLines(file)
    
    # Identifying first line that begins with blank
    for (i in 1:length(d.string)) {
      if (substring(d.string[i], 1, 1) == '') break
    }
    n.row <- i - 1

    # Dropping sample metainformation
    d.string <- paste(d.string[1:n.row], sep = '\n')

    # Converting numeric data to data.frame
    d <- read.csv(text = d.string, stringsAsFactors = FALSE, header = FALSE)

    # Stripping extra columns
    # (assuming that all data columns will have an entry in the 3rd row)
    n.col <- ncol(d)
    while (is.na(d[3, n.col]) || d[3, n.col] == '') {
      n.col <- n.col - 1
    }
    d <- d[, 1:n.col]
    
    # Dropping second row headers
    d <- d[-2, ] 

    # Modifying column names
    labels <- d[1, ][seq(1, n.col, by = 2)]
    colnames(d) <- paste(rep(labels, each = 2),
                         rep(c('emission', 'intensity'), n.col/2),
                         sep = '_')

    # Dropping first row headers   
    d <- d[-1, ]

    # Gathering all columns and filtering empty values
    d$row <- 1:(n.row - 2)
    m <- gather(d, id, value, -row) %>%
           filter(value != '') %>%
           separate(id, into = c('id', 'variable'), sep = '(?<=\\d)_(?=[ei])')

    n.meta <- meta %>% 
                group_by(scan) %>%
                mutate(n = ((stop - start)/data.interval + 1)*2)

    qualifiers <- meta[rep(1:nrow(meta), n.meta$n), 
                       c('sample', 'cell', 'scan', 'label',
                         'excitation.wavelength')]

    m <- cbind(qualifiers, m) %>%
           select(-id) %>%
           spread(variable, value) %>%
           select(-row, excitation = excitation.wavelength)

    out <- rbind(out, m)
  }

  return(out)
}
