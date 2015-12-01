# Test of stored data and plotting functions

library(dplyr)
library(ggplot2)
library(fluoroscripts)

d <- filter(protein.spectra, protein %in% c('EGFP', 'EYFP'))

p <- plot_spectra(d$wavelength, d$excitation, d$emission, d$protein,
                  lasers = c(488), channels = c('530/30', '670LP'))
p <- p + coord_cartesian(xlim = c(450, 800))
ggsave('spectra_test.pdf', width = 9, height = 5, units = 'in')

p <- plot_emissions(d$wavelength, d$excitation, d$emission, d$protein,
                    lasers = c(488), channels = c('530/30', '670LP'))
p <- p + coord_cartesian(xlim = c(450, 800))
ggsave('emission_test.pdf', width = 9, height = 5, units = 'in')
