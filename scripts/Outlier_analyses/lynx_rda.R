########################
### RDA for lynx WGR ###
########################

# ==============================================================================
# Load libraries 
# ==============================================================================

library(vegan)
library(adespatial)
library(spdep)
library(ape)
library(LEA)
library(vcfR)
library(caret)
library(dplyr)
library(UpSetR)
library(viridis)
library(corrplot)

# ==============================================================================
# Read in the interpolated genotype data - RDA can't handle missing data
# ==============================================================================
### This is the same input as LFMM:
Genotypes <- read.lfmm(input.file = "popgen/input_files/lynx_WGR_pruned_snps.lfmm_imputed.lfmm")

vcf_df <- read.vcfR("popgen/input_files/lynx_WGR_pruned_snps.vcf")
locus_names <- vcf_df@fix[, "ID"]  # Vector of SNP IDs

vcf <- readLines("popgen/input_files/lynx_WGR_pruned_snps.vcf")
samples <- vcf[grep("^#CHROM", vcf)] %>%
  strsplit("\t") %>%
  unlist() %>%
  tail(-9)  # first 9 columns are fixed VCF fields

# Clean the IDs for the repeated nature produced by PLINK
clean_ids <- ifelse(
  grepl("^(.*)_\\1$", samples),       # Match repeated pattern like A202_A202
  sub("_(.*)$", "", samples),         # Strip suffix
  samples                             # Otherwise keep as-is
)

row.names(Genotypes)<-clean_ids
colnames(Genotypes)<-locus_names

# Load metadata & order to match input
meta <- read.table(file="populations.txt")
samples<-as.data.frame(row.names(Genotypes))
colnames(samples) <- "Samples"
meta_order <- as.data.frame(merge(samples, meta, by.x = "Samples", by.y = "V1"))
meta_order <- meta_order[match(clean_ids, meta_order$Samples), ]

# Load environmental data
env <- read.csv("landscape_genomics/extracted_BioClim_variables_updated_nolatlon.csv", row.names = 1)
env_data <-  env[match(clean_ids, row.names(env)), ]

# Subset the environmental variables based on the analyses in the LFMM script:
cor_matrix <- cor(env_data, use = "pairwise.complete.obs")
high_corr <- findCorrelation(cor_matrix, cutoff = 0.8, names = TRUE, verbose = FALSE)
high_corr <- high_corr[high_corr != "bio_019"]
high_corr[11]<-"bio_07"
env_filtered <- env_data %>% select(-all_of(high_corr))
vir_colors <- viridis(200)

png(filename = "Plots/corr_plot_bioclim_variables.png", height = 6, width = 6, units = "in", res = 300)
corrplot(cor_matrix, type = "upper", tl.cex = 0.8,
         method = "square", order = "hclust", col = vir_colors)
dev.off()

env_scaled <- as.data.frame(scale(env_filtered))

env_scaled$sample_id <- row.names(env_scaled)

# Load GPS coordinates
coords <- read.csv(file = "lynx_gps_updated.csv")
coords <-  coords[match(clean_ids,coords$sample_id), ]

# Merge coords with env
expl <- as.data.frame(merge(coords, env_scaled, by = "sample_id"))
expl_order <- expl[match(clean_ids, expl$sample_id), ]

# ==============================================================================
# Run the RDA
# ==============================================================================
env_scaled$sample_id <- NULL
rdaout <- rda(Genotypes ~ ., data=env_scaled)

rdaout
RsquareAdj(rdaout)

summary(eigenvals(rdaout, model = "constrained"))

png(filename = "Plots/RDA_scree.png", height = 6, width = 8, units = "in", res = 300)
par(mfrow=c(1,1))
screeplot(rdaout)
dev.off()
# Retain 2 axes

#signif.full <- anova.cca(rdaout, parallel=getOption("mc.cores")) # default is permutation=999
#signif.full

#signif.axis <- anova.cca(rdaout, by="axis", parallel=getOption("mc.cores")) # took too long - skipping
#signif.axis

vif.cca(rdaout)

# Set levels BEFORE creating the pop factor
meta_order$V2 <- factor(meta_order$V2,
                        levels = c("Labrador", "Maine", "New_Brunswick", "Newfoundland", "Quebec_NSLR", "Quebec_SSLR"))

# Merge and make sure rows are matched to RDA sites
env_scaled$sample_id <- rownames(env_scaled)
env_pop <- merge(env_scaled, meta_order, by.x = "sample_id", by.y = "Samples")

# Ensure ordering matches rows of `rdaout`
env_pop <- env_pop[match(clean_ids, env_pop$sample_id), ]
pop <- env_pop$V2  # factor

# Color palette
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c")

# Plot with colors
png(filename = "plots/RDA_biplot.png", height = 6, width = 8, units = "in", res = 300)
plot(rdaout, type = "n", scaling = 3)
#points(rdaout, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
points(rdaout, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg[as.numeric(pop)])
text(rdaout, scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomright", legend = levels(pop), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
dev.off()

eig_vals <- rdaout$CCA$eig

# Calculate % variance explained
var_explained <- eig_vals / sum(eig_vals) * 100

# Print % for each axis
var_explained

# ==============================================================================
# Identify candidate SNPs involved in local adaptation
# ==============================================================================

load.rda <- scores(rdaout, choices=c(1:2), display="species")  # Species scores for the first four constrained axes

png(filename = "Plots/RDA_loadings_per_axis.png", height = 6, width = 10, units = "in", res = 300)
par(mfrow = c(1, 2))
hist(load.rda[,1], xlab="Loadings on RDA1", main = NA)
hist(load.rda[,2], xlab="Loadings on RDA2", main = NA)
dev.off()

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 1206
cand2 <- outliers(load.rda[,2],3) # 735


ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))

colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")

all_snps <- unique(c(cand1$snp, cand2$snp))

# Create binary presence/absence matrix
upset_data <- data.frame(
  SNP = all_snps,
  RDA1 = all_snps %in% cand1$snp,
  RDA2 = all_snps %in% cand2$snp
)

# UpSetR requires rownames
rownames(upset_data) <- upset_data$snp
upset_data$snp <- NULL

# Convert logicals to 1/0 (numeric)
upset_data[, -1] <- lapply(upset_data[, -1], as.integer)

# Create UpSet plot
png(filename = "Plots/RDA_UpSet_plot.png", width = 8, height = 6, units = "in", res = 300)
upset(upset_data, sets = c("RDA1", "RDA2"))
dev.off()

cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=8)  # 8 columns for 8 predictors
colnames(foo) <- c("bio_010","bio_012","bio_016","bio_019","bio_02","bio_03","bio_08","bio_09")

env_scaled$sample_id <- NULL
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- Genotypes[,nam]
  foo[i,] <- apply(env_scaled,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  # 9 duplicate detections

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1

table(foo[foo[,1]==2,2]) #  18 duplicates on axis 2

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

colnames(cand)[12] <- "predictor"
colnames(cand)[13] <- "correlation"

table(cand$predictor) 

write.csv(file = "landscape_genomics/RDA_loci_just_climate.csv", cand, quote = F, row.names = F)

