########################
### Lynx EEMS output ###
########################

#library("devtools")
#install_github("dipetkov/reemsplots2")

library("reemsplots2")
library("ggplot2")
library("sf")
library("dplyr")
library("rnaturalearth")
library("rnaturalearthdata")
library("viridis")

mcmcpath <- "H:/Shared drives/wildlife_genomics_lab/home_BB/project_lynx/EEMS/lynx_EEMS_unfilt_chain2"
plots <- make_eems_plots(mcmcpath, longlat = TRUE, dpi = 600)
names(plots)
#plots$qrates01
#plots$mrates02

# load EEMS mrates & qrates grid
mrates <- plots$mrates01$data
qrates <- plots$qrates01$data

# mrates has columns: x, y, z
# Example data frame of collection sites (replace with your own GPS coords)
sites <- read.csv(file = "lynx_gps.csv")

# Get map data (Canada + USA at medium resolution)
world <- ne_states(country = c("canada", "united states of america"), returnclass = "sf")

# Filter for Quebec, Newfoundland & Labrador, and Maine
regions <- world[world$name %in% c("Québec", "Newfoundland and Labrador", "Maine",
                                   "Prince Edward Island", "Nova Scotia", "New Brunswick",
                                   "Vermont", "New Hampshire", "New York", "Massachusetts"), ]

# Plot the mrates
map1<-ggplot() +
  geom_sf(data = regions, fill = "white", color = "black") +
  geom_point(data = sites, aes(x = longitude, y = latitude), 
             color = "black", size = 3) +
  coord_sf(xlim = c(-79.5, -52), ylim = c(43, 56), expand = FALSE) +
  theme_classic() +
  geom_raster(data = mrates, aes(x, y, fill = z), alpha = 0.7) +
  scale_fill_viridis_c(name = "log(m)", option = "viridis") +
  labs(x = "Longitude", y = "Latitude", title = "") +
  theme(
    panel.grid = element_line(color = "grey80", linetype = "dotted"),
    plot.title = element_text(hjust = 0.5)
  )
map1

ggsave(filename = "plots/EEMS_mrates.svg", map1, height = 8, width = 8)

# Plot the qrates
map2<-ggplot() +
  geom_sf(data = regions, fill = "white", color = "black") +
  geom_point(data = sites, aes(x = longitude, y = latitude), 
             color = "black", size = 3) +
  coord_sf(xlim = c(-79.5, -52), ylim = c(43, 56), expand = FALSE) +
  theme_classic() +
  geom_raster(data = qrates, aes(x, y, fill = z), alpha = 0.7) +
  scale_fill_viridis_c(name = "log(q)", option = "viridis") +
  labs(x = "Longitude", y = "Latitude", title = "") +
  theme(
    panel.grid = element_line(color = "grey80", linetype = "dotted"),
    plot.title = element_text(hjust = 0.5)
  )
map2

ggsave(filename = "plots/EEMS_qrates.svg", map2, height = 8, width = 8)

##############################################################################
install.packages("rEEMSplots",repos=NULL,type="source")
library(rEEMSplots)

#rEEMSplots requires a few other packages in order to work. If these are not already installed, an error message would appear.
library(Rcpp)
library(raster)
library(rgeos)
library(sp)
#install.packages("RcppEigen") #just missing this one
library(RcppEigen) #all set
#install.packages(c("Rcpp", "RcppEigen", "raster", "rgeos", "sp"))
#install.packages("rworldmap") #done
library("rworldmap")
#install.packages("rworldxtra") done
library(rworldxtra)
install.packages(GhostScript) #what is this

mcmcpath = "H:/Shared drives/wildlife_genomics_lab/home_BB/project_lynx/EEMS/data/lynx_EEMS_unfilt_chain1"
plotpath = "H:/Shared drives/wildlife_genomics_lab/home_BB/project_lynx/EEMS/EEMS_plots"

projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

eems.plots(mcmcpath, plotpath, longlat = TRUE, out.png = FALSE,
           plot.width=8,
           plot.height=6,
           add.grid = F,
           add.demes = T,
           projection.in = projection_none,
           projection.out = projection_mercator,
           min.cex.demes = 0.6,
           max.cex.demes = 1,
           col.demes = "red",
           pch.demes = "o",
           add.r.squared=F,
           add.abline = T,
           remove.singletons = F,
           add.map = T,
           col.map="black",
           lwd.map=0.2,
           add.title = F,
           res = 600,
           eems.colors = 
           )
###########################

