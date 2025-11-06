# Load libraries
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

setwd("H:/Shared drives/wildlife_genomics_lab/home_BB/project_lynx/")

# Example data frame of collection sites (replace with your own GPS coords)

sites <- read.csv(file = "lynx_gps.csv")
pop <- read.table(file="lynx_clusters_K4.txt")
sites2<-merge(sites, pop, by.x = "sample_id", by.y = "V1")

# Get map data (Canada + USA at medium resolution)
world <- ne_states(country = c("canada", "united states of america"), returnclass = "sf")

# Filter for Quebec, Newfoundland & Labrador, and Maine
regions <- world[world$name %in% c("Québec", "Newfoundland and Labrador", "Maine",
                                   "Prince Edward Island", "Nova Scotia", "New Brunswick",
                                   "Vermont", "New Hampshire", "New York", "Massachusetts"), ]

# Plot the inset
map1<-ggplot() +
  geom_sf(data = regions, fill = "grey90", color = "black") +
  geom_point(data = sites2, aes(x = longitude, y = latitude, color = V3), 
             size = 3) +
  coord_sf(xlim = c(-79.5, -52), ylim = c(43, 56), expand = FALSE) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", title = "") +
  theme(
    panel.grid = element_line(color = "grey80", linetype = "dotted"),
    plot.title = element_text(hjust = 0.5)
  )
map1

ggsave(filename = "plots/lynx_collection_map.svg", height = 8, width = 8)

myshape <- st_read("Lynx_distr/data_0.shp")

# Plot the region
world2 <- world[world$name != "Hawaii",]
map2 <- ggplot() +
  geom_sf(data = world2, fill = "grey90", color = "black") +
  geom_sf(data = myshape, fill = "orange", color = NA, size = 1) +   # overlay shapefile
  coord_sf(xlim = c(-180, -52), ylim = c(20, 90), expand = FALSE) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", title = "") +
  theme(
    panel.grid = element_line(color = "grey80", linetype = "dotted"),
    plot.title = element_text(hjust = 0.5)
  )
map2

ggsave(filename = "plots/lynx_US_map.svg", map2, height = 8, width = 8)

# define gradient from white to red
grad_fun <- colorRamp(c("white", "red"))

# your target values
vals <- c(0.03, 0.05, 0.08, 0.11, 0.16, 0.19)

# rescale to 0–1
scaled_vals <- vals / 0.2

# get RGB values and convert to hex
cols <- rgb(grad_fun(scaled_vals), maxColorValue = 255)

cols
