# Test of spectral area functions

library(dplyr)
library(fluoroscripts)

d <- filter(protein.spectra, protein %in% c('EGFP', 'EYFP', 'mOrange'))

areas <- d %>%
           group_by(protein) %>%
           do(calculate_area(.$wavelength, .$excitation, .$emission,
                             lasers = c(488), 
                             channels = c('530/30', '585/42', '670LP')))

#print(areas)

overlaps <- calculate_pairwise_overlap(
  d$wavelength, d$excitation, d$emission, d$protein, 
  lasers = c(488), channels = c('530/30', '585/42', '670LP'))

#print(overlaps)

combinations <- calculate_max_overlap(
  d$wavelength, d$excitation, d$emission, d$protein, 2,
  lasers = c(488), channels = c('530/30', '585/42', '670LP'),
  spectra.threshold = 0.5e-1)

combinations <- arrange(combinations, overlap)

print(combinations)
