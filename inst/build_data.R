# Script to generate data supplied with the package.

library(dplyr)
library(tidyr)

# Reading in aggregated data
protein.summary <- read.csv('protein_summary.csv')

# Reading in spectra files
protein.spectra <- c()
spectra.files <- list.files('./protein_spectra', pattern = '*.csv')

for (file in spectra.files) {
  d <- read.csv(file.path('protein_spectra', file))

  d.list <- list(d[ , c(1, 2)], d[ , c(3, 4)])

  for (i in c(1, 2)) {

    d <- d.list[[i]]
    colnames(d) <- c('wavelength', 'intensity')

    d <- filter(d, !is.na(wavelength))

    f_interpolate <- splinefun(d$wavelength, d$intensity)
    f_extrapolate <- function(x) {
      low <- min(d$wavelength)
      high <- max(d$wavelength)
      out <- ifelse((x > low) & (x < high), f_interpolate(x), 0)
      return(out)
    }

    x <- seq(300, 800, length.out = 200)
    y <- f_extrapolate(x)

    # Ensuring a scale of to 0-1
    y <- (y - min(y))/(max(y) - min(y))

    d.list[[i]] <- data.frame(wavelength = x, intensity = y)
  }

  d.list[[1]]$spectra <- 'r.excitation'
  d.list[[2]]$spectra <- 'r.emission'

  d <- rbind(d.list[[1]], d.list[[2]])
  d <- spread(d, spectra, intensity)

  protein <- gsub('.csv$', '', file)
  d$protein <- protein 

  if (!protein %in% protein.summary$protein) {
    msg <- sprintf('No summary data found for "%s"', label)
    stop(msg)
  }

  d <- d[ , c(4, 1, 3, 2)]
  protein.spectra <- rbind(protein.spectra, d)
}

protein.spectra <- mutate(protein.spectra,
                          protein = as.character(protein))

protein.summary <- mutate(protein.summary,
                          protein = as.character(protein))

# Multiplying spectral data by brightness
protein.spectra <- left_join(protein.spectra, 
                             protein.summary[, c('protein', 'brightness')], 
                             by = 'protein')

protein.spectra <- protein.spectra %>%
                     mutate(excitation = brightness * r.excitation,
                     emission = brightness * r.emission) %>%
                     select(-brightness)

save(protein.spectra, file = 'protein_spectra.RData')
save(protein.summary, file = 'protein_summary.RData')

