# Test of spectral area functions

library(dplyr)
library(fluoroscripts)

d <- filter(protein.spectra, protein %in% c('EGFP', 'EYFP'))
print(head(d))

areas <- d %>%
           group_by(protein) %>%
           do(calculate_area(.$wavelength, .$excitation, .$emission,
                             lasers = c(488), 
                             channels = c('530/30', '585/42', '670LP')))

print(areas)
