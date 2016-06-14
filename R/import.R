# Functions for importing data

#========================================================================>
# Misc functions

# Check specified filename glob
check_glob <- function(filename) {

  files <- Sys.glob(filename)
  if (length(files) == 0) {
    msg = sprintf('No files matching "%s"', filename)
    stop(msg)
  }

  return(files)
}

# Match instrument type
match_instrument <- function(instrument) {

  if (length(instrument) > 1) {
    msg = 'Instrument name must be provided as a single string'
    stop(msg)
  }

  # Trimming terminal whitespace
  instrument = str_trim(instrument)

  patterns = c('cary' = '[cC]ary\\s*(?=($|([eE]clipse$)))',
               'synergy4' = '(?<=^|([bB]io[tT]ek))\\s*[sS]ynergy\\s*4')

  matches = str_detect(instrument, patterns)

  valid_names <- c('Cary Eclipse', 'BioTek Synergy 4')
  valid_string <- paste(valid_names, collapse = ', ')

  if (sum(matches) == 2) {
    msg = sprintf('Ambiguous instrument name. Supported instruments:\n\t%s',
                   valid_string)
    stop(msg)
  }
  else if (sum(matches) == 0) {
    msg = sprintf('Invalid instrument name. Supported instruments:\n\t%s',
                   valid_string)
    stop(msg)
  }
  else {
    return(names(patterns)[matches])
  }
}

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
#
#' @return A data.frame with columns of scan parameters.
#'
#' @examples
#' # To do
#' @export

import_cary_meta <- function(filename, aggregate = FALSE) {

  files <- check_glob(filename) 

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
    if (!any(str_detect(d.string[split + c(6, 7)], 'Cary Eclipse'))) {
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
    
    # Row format requires an offset depending on how csv was generated
    if (str_detect(d.string[[1]], 'Z Axis')) {
      k <- 0
    } else {
      k <- 1
    }  

    # Finding rest of breaks
    while (TRUE) {

      # Extracting data
      if (substring(d.string[split], 1, 1) == '') {

        label <- str_trim(str_extract(d.string[split + k + 1], '.*?(?=_EX_)'))
        date <- str_extract(d.string[split + k + 2], '\\d.*')
        time <- dmy_hms(date)

        excitation <- str_trim(str_extract(d.string[split + k + 1], 
                                           '(?<=_EX_).*'))

        param <- rep(NA, 22)
        for (i in 1:22) {
          param[i] <- d.string[split + i + k + 5] %>%
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

    # Modifying PMT voltages to always be numeric
    d <- mutate(d,
                pmt.voltage = gsub('Low', '400', pmt.voltage),
                pmt.voltage = gsub('Medium', '600', pmt.voltage),
                pmt.voltage = gsub('High', '800', pmt.voltage))

    numeric.columns <- c('scan', 'time', 'start', 'stop', 
                         'excitation.wavelength', 'excitation.slit', 
                         'emission.slit', 'scan.rate', 'data.interval', 
                         'pmt.voltage', 'averaging.time', 'excitation.stop',
                         'excitation.increment')

    for (column in numeric.columns) {
      d[, column] <- as.numeric(d[, column])
    }

    # Adding filename
    d$filename <- basename(file)

    out <- rbind(out, d)
  }

  # Adding sample and cell numbers
  out <- arrange(out, time)

  row <- 1
  current.sample <- 0
  
  out$cell <- 0
  out$sample <- 0

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

    cells <- 1:n.cells
    index <- 1

    while (TRUE) {

      out$cell[row] <- cells[index]
      out$sample[row] <- current.sample + cells[index]

      row <- row + 1
     
      if (row >= nrow(out)) break
      if (out$excitation.wavelength[row + 1] < out$excitation.wavelength[row]) {
        break
      }

      if (index == n.cells) {
        index <- 1
      } else {
      index <- index + 1
      }
    }

    row <- row + 1

    current.sample <- current.sample + n.cells
  }

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

  out <- select(out, filename, everything())

  return(out)
}

#------------------------------------------------------------------------
#' Import numeric data from csv output from Cary Eclipse Spectrophotometer
#'
#' Imports numeric csv output from a Cary Eclipse Spectrophotometer and 
#' converts it into long format.
#'
#' @param filename Filename or pattern for matching multiple filenames.
#
#' @return A data.frame with the following columns: filename, 
#'         sample, cell, scan, label, excitation, emission, intensity
#'
#' @examples
#' # To do
#' @export

import_cary_data <- function(filename) {

  files <- check_glob(filename)

  # Initializing output
  out <- data.frame()

  # Looping through each file
  for (file in files) {

    # Getting metadata first
    meta <- import_cary_meta(file)

    # Reading file as string
    d.string <- readLines(file)
    
    # Identifying first line that begins with blank
    for (i in 1:length(d.string)) {
      if (substring(d.string[i], 1, 1) == '') break
    }
    n.row <- i - 1

    # Dropping sample metainformation
    d.string <- paste(d.string[1:n.row], sep = '\n')

    # Conversion of numeric data depends on how the csv was generated
    if (str_detect(d.string[[1]], 'Z Axis')) {

      # Converting numeric data to data.frame
      d <- read.csv(text = d.string, stringsAsFactors = FALSE, header = FALSE)

      # Stripping extra columns
      # (assuming that all data columns will have an entry in the 4th row)
      n.col <- ncol(d)
      while (is.na(d[4, n.col]) || d[4, n.col] == '') {
        n.col <- n.col - 1
      }
      d <- d[, 1:n.col]

      # Dropping top header and column numbers 
      d <- d[-c(1, 3), ] 

      # Renaming 
      colnames(d) <- c('emission', paste(d[1, 2:n.col], 2:n.col, sep = 'X'))

      # Dropping second row headers
      d <- d[-2, ] 

      # Gathering all columns and filtering empty values
      m <- gather(d, id, intensity, -emission) %>%
             mutate(id = gsub('X\\d+', '', id),
                    filename = basename(file)) %>%
             filter(intensity != '')

      n.meta <- meta %>% 
                  group_by(scan) %>%
                  mutate(n = ((stop - start)/data.interval + 1))

      qualifiers <- meta[rep(1:nrow(meta), n.meta$n), 
                         c('sample', 'cell', 'scan', 'label',
                           'excitation.wavelength')]

      m <- cbind(qualifiers, m)

    } else {

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
      colnames(d) <- paste(colnames(d), 1:n.col, sep = 'X')

      # Dropping first row headers   
      d <- d[-1, ]

      # Gathering all columns and filtering empty values
      d$row <- 1:(n.row - 2)
      m <- gather(d, id, value, -row) %>%
             mutate(id = gsub('X\\d+', '', id)) %>%
             filter(value != '') %>%
             separate(id, into = c('id', 'variable'), sep = '(?<=\\d)_(?=[ei])')

      n.meta <- meta %>%
                  group_by(scan) %>%
                  mutate(n = ((stop - start)/data.interval + 1)*2)

      qualifiers <- meta[rep(1:nrow(meta), n.meta$n), 
                         c('sample', 'cell', 'scan', 'label',
                           'excitation.wavelength')]

      print(class(qualifiers))
      print(class(m))

      m <- cbind(qualifiers, m) %>%
             spread(variable, value)

    }

    m <- m %>%
           mutate(emission = as.numeric(emission),
                  intensity = as.numeric(intensity),
                  filename = basename(file)) %>%
           select(filename, sample, cell, scan, label, 
                  excitation = excitation.wavelength, emission, intensity)

    out <- rbind(out, m)
  }

  return(out)
}

#------------------------------------------------------------------------
#' Import numeric data from table output of BioTek Synergy 4 
#'
#' Imports numeric table output from a BioTek Synergy 4 Spectrophotometer and 
#' converts it into long format.
#'
#' @param filename Filename or pattern for matching multiple filenames.
#
#' @return A data.frame with the following columns: filename, read, 
#'         well, row, column, excitation, emission, intensity
#'
#' @examples
#' # To do
#' @export

import_synergy4_data <- function(filename) {

  files <- check_glob(filename)

  # Initializing output
  out <- data.frame()

  # Looping through each file
  for (file in files) {

    # Reading file as string
    d.lines <- readLines(file)
    d.string <- paste(d.lines, collapse = '\n')
    
    # Splitting headers and content
    locations <- str_locate_all(d.string, '(?<=^|(\n\n)).+(?=\n\n)')[[1]]
    headers <- substring(d.string, locations[, 1], locations[, 2])

    content <- str_split(d.string, '(?<=^|(\n\n)).+(?=\n\n)')[[1]][-1]
    content <- str_trim(content)

    # Matching headers of interest
    table_matches <- str_detect(headers, 'Read \\d+:.+ Spectrum')
    table_headers <- headers[table_matches]
    table_content <- content[table_matches]

    # Looping through read tables
    for (i in 1:length(table_headers)) {
      header <- table_headers[i]
      read <- as.numeric(str_extract(header, '(?<=Read )\\d+'))
      d <- read_delim(table_content[i], '\t') %>%
           gather('well', 'intensity', -Wavelength) %>%
           mutate(filename = basename(file),
                  read = read,
                  temp = well,
                  Wavelength = as.numeric(str_replace(
                                 Wavelength, 'nm', ''))) %>%
           separate(temp, c('row', 'column'), sep = '(?<=\\D)') %>%
           mutate(column = as.numeric(column))

      if (str_detect(header, 'EM')) {
        d <- rename(d, emission = Wavelength)
        d$excitation <- NA
      }
      else {
        d <- rename(d, excitation = Wavelength)
        d$emission <- NA
      }

      d <- select(d, filename, read, well, row, column,
                     excitation, emission, intensity)

      out <- rbind(out, d)
    }
  } 

  return(out)
}

#========================================================================>
# Generic functions

#------------------------------------------------------------------------
#' Import sample description metadata from spectrophotometer output.
#'
#' Imports sample description metadata from exported output and converts it 
#' into long format. Import procedure depends on specified instrument. 
#' [This function is a placeholder as only Cary Eclipse import is currently
#' supported]. Calls import_cary_meta().
#'
#' @param filename Filename or pattern for matching multiple filenames.
#' @param instrumet String specifying instrument. Currently limited to
#'                  'Cary Eclipse' or 'cary'.
#' @param ... Arguments passed to specific instrument import functions.
#
#' @return A data.frame with columns of scan parameters.
#'
#' @examples
#' # To do
#' @export

import_meta <- function(filename, instrument = 'cary', ...) {

  instrument = match_instrument(instrument)

  functions <- c('cary' = 'import_cary_meta')

  if (!instrument %in% names(functions)) {
    msg = 'Metadata input for this instrument is not supported'
    stop(msg)
  }

  do.call(functions[instrument], c(list(filename = filename), list(...)))
} 

#------------------------------------------------------------------------
#' Import numeric data from spectrophotometer output.
#'
#' Imports numeric data from exported output and converts it 
#' into long format. Import procedure depends on specified instrument. 
#' Calls import_cary_data() or import_synergy4_data().
#'
#' @param filename Filename or pattern for matching multiple filenames.
#' @param instrumet String specifying instrument. Currently limited to
#'                  'Cary Eclipse' or 'cary', 'BioTek Synergy 4' or 
#'                  'synergy4'.
#' @param ... Arguments passed to specific instrument import functions.
#
#' @return A data.frame with columns of scan parameters.
#'
#' @examples
#' # To do
#' @export

import_data <- function(filename, instrument = 'cary', ...) {

  instrument = match_instrument(instrument)

  functions <- c('cary' = 'import_cary_data',
                 'synergy4' = 'import_synergy4_data')

  if (!instrument %in% names(functions)) {
    msg = 'Numeric data input for this instrument is not supported'
    stop(msg)
  }

  do.call(functions[instrument], c(list(filename = filename), list(...)))
} 
