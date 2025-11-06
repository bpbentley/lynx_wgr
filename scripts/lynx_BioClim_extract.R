##############################################
### Extract BioClim variables from Rasters ###
##############################################

library(terra)
library(tidyverse)
library(geosphere)
library(sf)

# ==== Initial extraction ====

# 1. Define coordinates (see below for fixing the NAs)
coords <- read.csv(file = "lynx_gps_updated.csv")

# 2. Convert to spatial points
points <- vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")

# 3. List raster files (assumes filenames like 'wc2.1_2.5m_bio_1.tif', ..., 'wc2.1_2.5m_bio_19.tif')
rasters <- list.files("landscape_genomics/wc2.1_30s_bio/", full.names = TRUE)

# 4. Loop to extract values
extracted_list <- lapply(rasters, function(rfile) {
  r <- terra::rast(rfile)
  v <- terra::extract(r, points)[, 2]  # skip ID column
  varname <- paste0("bio", gsub("\\D", "", basename(rfile)))  # e.g. "bio1"
  tibble(!!varname := v)
})

# 5. Combine and add coordinates
extracted_data <- bind_cols(extracted_list)
result <- bind_cols(coords, extracted_data)

colnames(result)<-gsub("213","_",colnames(result))
#write.csv(file = "landscape_genomics/extracted_BioClim_variables_updated.csv", result,
#          quote = F, row.names = F)

r <- terra::rast("landscape_genomics/wc2.1_30s_bio/wc2.1_30s_bio_7.tif")

points <- vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")

vals <- terra::extract(r, points)

# Identify NA points
coords$has_NA <- is.na(vals[,2])

plot(r, main = "bio1 with NA points highlighted", xlim = c(-80,-55), ylim = c(44.5,54.5))
points(coords$longitude, coords$latitude, pch = 16, col = ifelse(coords$has_NA, "red", "blue"))
#legend("bottomleft", legend = c("Valid", "NA"), col = c("blue", "red"), pch = 16)

## There are 3 samples with NA points, want to avoid this:

# ==== Troubleshooting NAs ====
r <- terra::rast("landscape_genomics/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")

points <- vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")

vals <- terra::extract(r, points)

# Identify NA points
coords$has_NA <- is.na(vals[,2])

# Plot with NA points highlighted
plot(r, main = "bio1 with NA points highlighted", xlim = c(-80,-55), ylim = c(44.5,54.5))
points(coords$longitude, coords$latitude, pch = 16, col = ifelse(coords$has_NA, "red", "blue"))
legend("bottomleft", legend = c("Valid", "NA"), col = c("blue", "red"), pch = 16)

# All 3 points fall in ocean habitat where BioClim variables aren't available.
## Goal to drag those GPS points to the nearest raster cell:

# ==== Updating GPS points to nearest raster cell ====
fix_oceanic_points_fast <- function(coords, raster_layer, search_radius = 3) {
  # This method is more memory efficient for large rasters
  # search_radius is in number of cells to search around
  
  # Ensure coords is a proper data frame or matrix
  coords <- as.data.frame(coords)
  if(ncol(coords) != 2) {
    stop("Coordinates must have exactly 2 columns (longitude, latitude)")
  }
  
  # Make sure column names are proper
  names(coords) <- c("x", "y")
  
  # Extract values - simplified approach
  extracted_values <- terra::extract(raster_layer, coords)
  
  # Handle different extract return formats
  if(is.data.frame(extracted_values) || is.matrix(extracted_values)) {
    if(ncol(extracted_values) > 1) {
      na_indices <- which(is.na(extracted_values[,2]))  # Skip ID column
    } else {
      na_indices <- which(is.na(extracted_values[,1]))
    }
  } else {
    na_indices <- which(is.na(extracted_values))
  }
  
  if(length(na_indices) == 0) {
    message("No oceanic points found!")
    return(coords)
  }
  
  message(paste("Found", length(na_indices), "oceanic points to fix"))
  fixed_coords <- coords
  
  for(i in na_indices) {
    original_point <- as.numeric(coords[i, ])
    
    # Get cell number for the oceanic point
    cell_num <- cellFromXY(raster_layer, matrix(original_point, ncol = 2))
    
    # Get adjacent cells in expanding rings
    found_replacement <- FALSE
    for(radius in 1:search_radius) {
      
      if(radius == 1) {
        # Get immediate neighbors
        adjacent_cells <- adjacent(raster_layer, cell_num, directions = "queen", 
                                   pairs = FALSE)
      } else {
        # For larger radii, use a different approach
        # Get coordinates of the center cell
        center_xy <- xyFromCell(raster_layer, cell_num)
        
        # Calculate search distance in map units
        res_x <- res(raster_layer)[1]
        res_y <- res(raster_layer)[2]
        search_dist_x <- radius * res_x
        search_dist_y <- radius * res_y
        
        # Define search window
        search_ext <- ext(center_xy[1] - search_dist_x, 
                          center_xy[1] + search_dist_x,
                          center_xy[2] - search_dist_y, 
                          center_xy[2] + search_dist_y)
        
        # Get cells in search window
        adjacent_cells <- cells(raster_layer, search_ext)
      }
      
      # Remove original cell and invalid cells
      adjacent_cells <- adjacent_cells[!is.na(adjacent_cells)]
      adjacent_cells <- adjacent_cells[adjacent_cells != cell_num]
      
      if(length(adjacent_cells) > 0) {
        # Get coordinates for adjacent cells
        adj_coords <- xyFromCell(raster_layer, adjacent_cells)
        
        # Get values for adjacent cells
        adj_values <- terra::extract(raster_layer, adj_coords)
        
        # Handle different extract return formats
        if(is.data.frame(adj_values) || is.matrix(adj_values)) {
          if(ncol(adj_values) > 1) {
            valid_indices <- which(!is.na(adj_values[,2]))  # Skip ID column
          } else {
            valid_indices <- which(!is.na(adj_values[,1]))
          }
        } else {
          valid_indices <- which(!is.na(adj_values))
        }
        
        if(length(valid_indices) > 0) {
          # Calculate distances and find closest
          valid_coords <- adj_coords[valid_indices, , drop = FALSE]
          distances <- sqrt((valid_coords[, 1] - original_point[1])^2 + 
                              (valid_coords[, 2] - original_point[2])^2)
          
          closest_coord <- valid_coords[which.min(distances), ]
          
          # Update coordinates
          fixed_coords[i, 1] <- closest_coord[1]
          fixed_coords[i, 2] <- closest_coord[2]
          
          message(paste("Point", i, "moved from (", 
                        round(original_point[1], 4), ",", round(original_point[2], 4), 
                        ") to (", 
                        round(closest_coord[1], 4), ",", round(closest_coord[2], 4), ")"))
          found_replacement <- TRUE
          break
        }
      }
    }
    
    if(!found_replacement) {
      warning(paste("Could not find land cell within", search_radius, 
                    "cells for point", i, "- keeping original coordinates"))
    }
  }
  
  return(fixed_coords)
}

bioclim_raster <- terra::rast("landscape_genomics/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")
gps_coords <- read.csv(file = "lynx_gps.csv")[,c(2:3)]

# Method 1 (more thorough but slower):
# fixed_coords <- fix_oceanic_points(gps_coords, bioclim_raster, max_distance = 0.05)

# Method 2 (faster, recommended for large datasets):
fixed_coords <- fix_oceanic_points_fast(gps_coords, bioclim_raster, search_radius = 20)

#Found 3 oceanic points to fix
#Point 36 moved from ( -64.2705 , 48.7234 ) to ( -64.2792 , 48.7125 )
#Point 37 moved from ( -66.2877 , 49.9455 ) to ( -66.3958 , 50.0875 )
#Point 39 moved from ( -66.1959 , 49.2173 ) to ( -66.1792 , 49.2042 )

# Updated the GPS file manually with new points.