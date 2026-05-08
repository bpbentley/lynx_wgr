#####################################################
### Extract future BioClim variables from Rasters ###
#####################################################

library(terra)
library(tidyverse)
library(sf)

# 1. Load coordinates
coords <- read.csv("lynx_gps_updated.csv")
points <- vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")

# 2. List multi-band raster files for each time point
# Assumes filenames include the year, e.g., "bio_2040.tif"
raster_files <- list.files(
  "landscape_genomics/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245/",
  pattern = "\\.tif$", full.names = TRUE
)

# 3. Extract all 19 variables from each raster for each time point
extracted_all <- lapply(raster_files, function(f) {
  # Extract year from filename
  year <- str_extract(basename(f), "(?<=-)\\d{4}")
  
  # Load raster and extract values
  r <- rast(f)  # multi-band raster
  vals <- terra::extract(r, points)[, -1]  # drop ID column
  
  # Rename columns to bio_01_2040, ..., bio_019_2100
  colnames(vals) <- paste0("bio_0", 1:ncol(vals), "_", year)
  
  as_tibble(vals)
})

# 4. Combine all extracted values and bind with coordinates
result <- bind_cols(coords, bind_cols(extracted_all))

# 5. Optional: Save to CSV
write.csv(result, "landscape_genomics/extracted_future_bioclim_all.csv", row.names = FALSE, quote = FALSE)
