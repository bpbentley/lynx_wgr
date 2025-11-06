###############################
### Leveraging BioClim Data ###
###############################

library(adegenet)
library(vcfR)
library(vegan)
library(gdm)
library(raster)
library(dplyr)
library(ape)
library(tidyr)
library(stringr)
library(adespatial)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggpubr)

# ==== Read VCF and convert to genlight ====

vcf <- read.vcfR("popgen/input_files/lynx_WGR_adaptive_snps.vcf")
gen <- vcfR2genlight(vcf)

# ==== Add sample metadata ====

meta1 <- read.csv("landscape_genomics/extracted_BioClim_variables_updated.csv")  # Ensure it includes all required variables
meta2 <- read.csv("landscape_genomics/extracted_future_bioclim_all.csv")
meta1$longitude <- meta1$latitude <- NULL
meta3 <- merge(meta1, meta2, by = "sample_id")
base_vars <- c("bio_010","bio_012","bio_016","bio_019","bio_02","bio_03","bio_08","bio_09") # From correlations (see LFMM script)
pattern <- paste0("^(", paste(base_vars, collapse = "|"), ")")
row.names(meta3)<-meta3$sample_id
filtered_df <- meta3 %>% dplyr::select(matches(pattern))

clean_ids <- ifelse(
  grepl("^(.*)_\\1$", indNames(gen)),       # Match repeated pattern like A202_A202
  sub("_(.*)$", "", indNames(gen)),         # Strip suffix
  indNames(gen)                             # Otherwise keep as-is
)

indNames(gen) <- clean_ids
final_df <- filtered_df[clean_ids, ]
gps <- read.csv(file="lynx_gps_updated.csv")
row.names(gps) <- gps$sample_id
gps <- gps[clean_ids,]

final_df$latitude <- gps$latitude
final_df$longitude <- gps$longitude

#write.csv(file="landscape_genomics/gdm/gdm_all_variable_input.csv", final_df, quote = F)
#current_vars <- read.csv(file = "landscape_genomics/gdm/gdm_input_current_variables.csv")
current_vars <- bind_cols(row.names(final_df), final_df[,c(1:8)])
colnames(current_vars)[1] <- "sample_id"
current_vars1 <- cbind(current_vars, gps$longitude, gps$latitude)
colnames(current_vars1)[c(10,11)] <- c("longitude", "latitude")

# ==== Calculate genetic distance ====
gen_dist <- dist(gen, diag = T, upper = T)  # Euclidean distance between individuals based on adaptive loci
gen_dist1 <- gen_dist/(max(gen_dist))  #rescale by dividing by max value
gen_dist1 <- as.matrix(gen_dist1)
gen_dist1 <- cbind(sample_id = rownames(gen_dist1), gen_dist1)  #make the row names an actual column named "ID"
rownames(gen_dist1) <- NULL  #remove prior row names

# ==== Run the GDM on current conditions ====
## N.B. using the filtered BioClim set of only uncorrelated variables
## Need to convert the site data to numeric and ensure all are in the correct order

# Match the sample ID against a corresponding number in both genetic df and env df
num_id <- as.data.frame(cbind(gen_dist1[,1], seq(1:nrow(current_vars1))))
gen_dist1[,1] <- as.numeric(num_id$V2[match(gen_dist1[,1], num_id$V1)])
colnames(gen_dist1) <- as.numeric(num_id$V2[match(colnames(gen_dist1), num_id$V1)])
current_vars1$sample_id <- as.numeric(num_id$V2[match(current_vars1$sample_id, num_id$V1)])

# change the name to site because ?
colnames(current_vars1)[1] <- "site"
colnames(gen_dist1)[1] <- "site"

# Make sure the matrix is NUMERIC (this is where the major issue was arising)
gen_dist1 <- as.data.frame(gen_dist1)
gen_dist1[] <- lapply(gen_dist1, as.numeric)

# Format for input to the model
gdm.input <- formatsitepair(
  bioData=gen_dist1,
  bioFormat=3,
  predData=current_vars1,
  siteColumn="site",
  XColumn="longitude",
  YColumn="latitude")

# Run the model
gdm <- gdm(gdm.input, geo = T, splines = NULL, knots = NULL)
save(gdm, file = "landscape_genomics/gdm/gdm_model.rda")
summary(gdm)
gdm$explained #70.72% - seems high?

gdm.importance <- gdm.varImp(gdm.input, geo=T, splines=NULL, nPerm=50, parallel=F, cores=1)
imp <- gdm.importance$`Predictor Importance`
imp$predictor <- rownames(imp)
colnames(imp)[1]<-"importance"
imp <- imp[order(imp$importance, decreasing = F),]

bp1<-ggbarplot(data = imp, x = "predictor", y = "importance", orientation = "horiz",
          ylab = "Percent Deviance Explained", xlab = "Predictor Variable") +
  font("xlab", face = "bold") + font("ylab", face = "bold")
bp1

ggsave(filename = "plots/gdm_variance_barplot.svg", bp1, height = 6, width = 8, units = "in")

plotUncertainty(gdm.input, sampleSites=0.70, bsIters=100, geo=T, plot.layout=c(3,4)) #2 pages because there's too many
 
# ==== Extract the splines and plot them ====

splines<-isplineExtract(gdm)
#Geographic bio_010 bio_012 bio_016 bio_019 bio_02 bio_03 bio_08 bio_09

# BIO_09
bio_09<-as.data.frame(cbind(splines$x[,9],splines$y[,9]))
bp2<-ggplot(bio_09, aes(x = V1, y = V2)) +
  ylim(0,1) +
  geom_line(color = "steelblue") +
  labs(x = "Mean Temperature of Driest Quarter",
       y = "Partial ecological distance") + theme_classic()
bp2

# Geography
geo<-as.data.frame(cbind(splines$x[,1],splines$y[,1]))
bp3<-ggplot(geo, aes(x = V1, y = V2)) +
  ylim(0,1) +
  geom_line(color = "steelblue") +
  labs(x = "Geographic",
       y = "Partial ecological distance") + theme_classic()
bp3

# BIO_03
bio_03<-as.data.frame(cbind(splines$x[,7],splines$y[,7]))
bp4<-ggplot(bio_03, aes(x = V1, y = V2)) +
  ylim(0,1) +
  geom_line(color = "steelblue") +
  labs(x = "Isothermality",
       y = "Partial ecological distance") + theme_classic()
bp4

# BIO_10
bio_10<-as.data.frame(cbind(splines$x[,2],splines$y[,2]))
bp5<-ggplot(bio_10, aes(x = V1, y = V2)) +
  ylim(0,1) +
  geom_line(color = "steelblue") +
  labs(x = "Mean Temperature of Warmest Quarter",
       y = "Partial ecological distance") + theme_classic()
bp5

# Combine plot and save to file
combbp<-ggarrange(bp2, bp3, bp4, bp5)

ggsave(filename = "plots/gdm_splines.svg", combbp, height = 6, width = 8, units = "in")

# ==== Trialing the GDM predict function ====
# Current conditions raster stack:
# Updated variable list
sig_variables<- c("bio_010","bio_012","bio_016","bio_019","bio_02","bio_03","bio_08","bio_09")

# Extract numeric codes
codes <- gsub("bio", "", sig_variables)

# Construct filenames (assuming filenames follow "wc2.1_10m_bio_13.tif")
files <- sprintf("landscape_genomics/wc2.1_30s_bio/wc2.1_30s_bio_%s.tif", codes)
files <- gsub("bio__0", "bio_", files)

# Stack them using raster
bio_stack <- stack(files)
names(bio_stack)
new_names <- gsub("wc2.1_30s_bio_", "bio_0", names(bio_stack))
names(bio_stack) <- new_names

future_stack <- stack("landscape_genomics/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2061-2080.tif")
names(future_stack)

layer_numbers <- gsub(".*_", "", names(future_stack))
new_names <- paste0("bio_0", seq(1:19))
names(future_stack) <- new_names
names(future_stack)

future_sig <- subset(future_stack, sig_variables)  # raster

extent_region <- extent(-80, -50, 42, 60)

current_sig_crop <- crop(bio_stack, extent_region)
future_sig_crop <- crop(future_sig, extent_region)
beepr::beep(2)
plot(current_sig_crop[[8]], main = names(current_sig_crop)[8])
plot(future_sig_crop[[1]], main = names(future_sig_crop)[1])

load("landscape_genomics/gdm/gdm_model.rda")

pred.selection <- predict(gdm,
                          current_sig_crop,
                          time=T,
                          predRasts = future_sig_crop) # offset

sppdist  <- st_read("landscape_genomics/gdm/lynx_redlist/data_0.shp")

shape_vect <- vect(sppdist)

r_cropped <- crop(pred.selection, shape_vect)
r_masked <- mask(r_cropped, shape_vect)

# Load raster and convert to dataframe
r_df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE)
names(r_df)[3] <- "value"  # Rename for easier ggplot use

# Load country borders (North America)
na_countries <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")

# Load USA states and Canadian provinces
states <- ne_states(country = "United States of America", returnclass = "sf")
provinces <- ne_states(country = "Canada", returnclass = "sf")

# Ensure CRS match
species_dist <- st_transform(sppdist, crs(r_masked))

# Plot
C<-ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(option = "viridis", na.value = "transparent") +
  geom_sf(data = states, fill = NA, color = "grey80", size = 0.3) +
  geom_sf(data = provinces, fill = NA, color = "grey80", size = 0.3) +
  geom_sf(data = na_countries, fill = NA, color = "grey50", size = 0.5) +
  geom_sf(data = species_dist, fill = NA, color = "black", size = 0.7) +
  coord_sf(xlim = c(-82, -50), ylim = c(43, 60), expand = FALSE) +
  theme_classic() +
  labs(fill = "Genomic offset", title = "2061-2080",
       x = "Longitude", y = "Latitude")
C

E <- ggarrange(A, B, C, D, ncol = 2, nrow = 2, common.legend = T)
ggsave(file="plots/lynx_genomic_offsets.png", E, height = 12, width = 12, units = "in")




# ==== Future climate conditions ====
future_vars <- read.csv(file = "landscape_genomics/gdm/gdm_input_future_variables.csv")
colnames(future_vars)[1] <- "site"
future_vars$site <- as.numeric(num_id$V2[match(future_vars$site, num_id$V1)])

gdm.input.future <- formatsitepair(bioData=gen_dist1, bioFormat=3,
                                   predData=future_vars, siteColumn="site",
                                   XColumn="longitude", YColumn="latitude") 

gdm.future <- gdm(gdm.input.future, geo = T, splines = NULL, knots = NULL)
summary(gdm.input.future)
save(gdm.future, file = "gdm_model_future.rda")
gdm.future$explained

gdm.importance.future <- gdm.varImp(gdm.input.future, geo=T, splines=NULL, nPerm=50, parallel=F, cores=1)

par(mfrow=c(1,1))
plot(gdm.future$observed, gdm.future$predicted, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))

