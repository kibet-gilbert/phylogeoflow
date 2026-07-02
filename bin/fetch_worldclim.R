install.packages("geodata")
library(geodata)

# Example: Fetch global bioclimatic variables at a 5-minute resolution
bio_data <- worldclim_global(var = "bio", res = 5, path = "path/to/save/directory")

