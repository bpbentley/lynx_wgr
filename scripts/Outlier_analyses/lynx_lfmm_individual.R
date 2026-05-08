#########################
### LFMM For Lynx WGR ###
#########################

# ==== Load libraries ====

library(LEA)
library(tidyverse)
library(corrr)
library(corrplot)
library(mice)
library(viridis)
library(caret)
library(ggpubr)
library(lfmm)
library(factoextra)
library(vcfR)
library(ggvenn)
library(UpSetR)
library(qvalue)

# ==== 1. Read in data and format for downstream use ====

# vcf2geno doesn't retain sample or locus IDs, so pull from VCF
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

# Read in the geno file produced from the vcf with vcf2geno
# Apply the sample names to the rows.
genofile <- read.geno("popgen/input_files/lynx_WGR_pruned_snps.geno")
rownames(genofile) <- clean_ids

# Load environmental data
env <- read.csv("landscape_genomics/extracted_BioClim_variables_updated_nolatlon.csv", row.names = 1)

# ==== 2. Process the environment data to remove highly correlated variables ====

# Remove highly correlated predictors
# Visualize correlation structure
cor_matrix <- cor(env, use = "pairwise.complete.obs")
vir_colors <- viridis(200)
corrplot(cor_matrix, type = "upper", tl.cex = 0.8,
         method = "square", order = "hclust", col = vir_colors)

# Remove highly correlated predictors (|r| > 0.8)
high_corr <- findCorrelation(cor_matrix, cutoff = 0.8, names = TRUE, verbose = FALSE)
high_corr <- high_corr[high_corr != "bio_019"]
high_corr[11]<-"bio_07"
env_filtered <- env %>% select(-all_of(high_corr))

# Scale the environmental variables
env_scaled <- scale(env_filtered)

# Retain only IDs in the genotype file (2 extra in env)
common_ids <- intersect(rownames(genofile), rownames(env_scaled))
env_scaled <- env_scaled[common_ids, ]

# ==== 3. Impute missing genotypes with snmf ====

# Convert the missing data "9" to NA and count how many missing genotypes
genofile[genofile == 9] <- NA
sum(is.na(genofile)) ## 1,011,057 to be imputed
sum(genofile, na.rm = T)

miss_geno <- as.data.frame(rowSums(is.na(genofile))/rowSums(!is.na(genofile)) * 100)
colnames(miss_geno)<-"Percent_missing"
write.table(file="missing_genotypes.txt", miss_geno, quote = F, sep = "\t")

# Use snmf:
setwd("popgen/input_files/")
lynx_WGR_pruned_snps.snmfProject<-load.snmfProject("lynx_WGR_pruned_snps.snmfProject")
impute(lynx_WGR_pruned_snps.snmfProject, "lynx_WGR_pruned_snps.geno", method = 'mode', K = 4, run = 1)
setwd("H:/Shared drives/wildlife_genomics_lab/home_BB/project_lynx")

# ==== 4. Create PCA composites of BioClim variables as an alternative [optional] ====

# Pull out the 55 samples from the total sample list
env2 <- env[rownames(env) %in% common_ids,]
env_matrix <- as.matrix(env2)

# Scale the data
env2_scaled <- scale(env_matrix)

# Run the PCA
env_pca <- prcomp(env2_scaled, center = FALSE, scale. = FALSE)
summary(env_pca)
plot(env_pca)
env_pca_vals<-as.data.frame(env_pca$x)
ggscatter(data = env_pca_vals, x = "PC1", y = "PC2")

# Scree-plot for choosing value of K:
png(filename = "plots/Scree_plot_PCA_composites_for_LFMM.png",
    width = 8, height = 6, units = "in", res = 300)
plot(env_pca, type = "lines", main = "Scree Plot of BioClim PCA")
dev.off()

# Total contributions
var_contr<-as.data.frame(summary(env_pca)$importance[3, ]) # Based on above, retain 4 PCs
var_contr$PC<-seq(1:nrow(var_contr))
colnames(var_contr)<-c("Cumulative_variance", "PC")
ggline(var_contr, x = "PC", y = "Cumulative_variance")

env_pcs <- env_pca$x[, 1:4]

# ==== 5a. Run the LFMM analysis on the 8 remaining BioClim variables ====
# Pull in the imputed file:
geno_imputed <- read.lfmm(input.file = "popgen/input_files/lynx_WGR_pruned_snps.lfmm_imputed.lfmm")
row.names(geno_imputed)<-clean_ids
colnames(geno_imputed)<-locus_names

# Using the eight uncorrelated BioClim variables:
project_1 <- lfmm_ridge(Y = geno_imputed, X = env_scaled, K = 3)

# Run the LFMM program:
pv1 <- lfmm_test(Y = geno_imputed, 
                X = env_scaled, 
                lfmm = project_1,
                calibrate = "gif")

# Extract the p-values
pvalues1 <- as.data.frame(pv1$calibrated.pvalue)

# FDR correction:
qvals1 <- as.data.frame(apply(pvalues1, 2, function(p) qvalue(p)$qvalues))

qvals1$bio_010_P<--log10(qvals1$bio_010)
qvals1$signif <- ifelse(qvals1$bio_010_P > 2, "significant", "non-significant")
table(qvals1$signif)

outlier_pal<-c("black", "red")
p1<-ggscatter(data = qvals1, x = "SNP", y = "bio_010_P", color = "signif",
              pal = outlier_pal,
              title = "LFMM on BioClim variables: bio_010",
              subtitle = paste0("Significance at p < 0.001. Sig. SNPs = ",table(qvals1$signif)[2]),
              xlab = "SNP", ylab = "-log10(p-value)") +
  font("xlab", face = "bold") + font("ylab", face = "bold")
p1
ggsave(file = "plots/bio_010_LFMM_composite.png", p1, height = 6, width = 8, units = "in")

sig_loci_list <- lapply(1:8, function(i) {
  sig_idx <- which(qvals1[, i] < 0.05)
  data.frame(
    BIO = colnames(qvals1)[i],
    Locus = rownames(qvals1)[sig_idx],
    Index = sig_idx
  )
})

sig_loci_bio_df <- do.call(rbind, sig_loci_list)

sig_loci_bio_unique <- sig_loci_bio_df %>%
  group_by(Locus) %>%
  summarise(BIOs = paste(unique(BIO), collapse = ","), .groups = "drop")

table(sig_loci_bio_unique$BIOs)

# Number of loci:
nrow(sig_loci_bio_unique)

# Plot the overlap with UpSetR:
# Get unique PCs and unique loci
all_bio <- unique(sig_loci_bio_df$BIO)

# Pivot to wide format
binary_df <- sig_loci_bio_df %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(names_from = BIO, values_from = value, values_fill = 0)

# Set Locus as row names
binary_mat <- as.data.frame(binary_df)
rownames(binary_mat) <- binary_mat$Locus
binary_mat$Locus <- NULL

png(file = "plots/upset_plot_BioClim_LFMM.png",
    height = 6, width = 8, units = "in", res = 300)
upset(binary_mat, 
      sets = sort(all_bio), 
      nintersects = 30,
      order.by = "freq")
dev.off()

# Write loci to file
write.csv(file="landscape_genomics/LFMM_loci_BioClim_variables_FDR0.05.csv", sig_loci_bio_unique,
          quote = F, row.names = F)

# ==== 5b. Run the LFMM analysis on the PCA composites ====

# Using the PCAs (and the lfmm file from 5a)
project_2 <- lfmm_ridge(Y = geno_imputed, X = env_pcs, K = 3)
pv2 <- lfmm_test(Y = geno_imputed, 
                X = env_pcs, 
                lfmm = project_2, 
                calibrate = "gif")
pvalues2 <- as.data.frame(pv2$calibrated.pvalue)

# FDR correction:
qvals2 <- as.data.frame(apply(pvalues2, 2, function(p) qvalue(p)$qvalues))

qvals2$PC1_P<--log10(qvals2$PC1)
qvals2$signif <- ifelse(qvals2$PC1 < 0.05, "significant", "non-significant")
table(qvals2$signif)

outlier_pal<-c("black", "red")
p1<-ggscatter(data = qvals2, x = "SNP", y = "PC1_P", color = "signif",
          pal = outlier_pal,
          title = "LFMM on composite environmental variables: PC1",
          subtitle = paste0("Significance at p < 0.01. Sig. SNPs = ",table(qvals2$signif)[2]),
          xlab = "SNP", ylab = "-log10(p-value)") +
  font("xlab", face = "bold") + font("ylab", face = "bold")
ggsave(file = "plots/PC1_LFMM_composite.png", p1, height = 6, width = 8, units = "in")

fviz_pca_biplot(env_pca,
                repel = TRUE,     # avoid label overlap
                col.var = "contrib",  # color arrows by contribution
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                col.ind = "gray40",   # color of individuals
                label = "var",        # only show variable labels
                pointsize = 2)

sig_loci_list_comp <- lapply(1:4, function(i) {
  sig_idx <- which(qvals2[, i] < 0.05)
  data.frame(
    PC = paste0("PC", i),
    Locus = rownames(qvals2)[sig_idx],
    Index = sig_idx
  )
})

names(sig_loci_list_comp) <- paste0("PC", 1:4)
sig_loci_df <- do.call(rbind, sig_loci_list_comp)

sig_loci_unique <- sig_loci_df %>%
  group_by(Locus) %>%
  summarise(PCs = paste(unique(PC), collapse = ","), .groups = "drop")

venn_list <- sig_loci_df %>%
  group_by(PC) %>%
  summarise(Loci = list(unique(Locus))) %>%
  deframe()

v1<-ggvenn(venn_list,
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#D55E00"),
       stroke_size = 0.3, set_name_size = 4)
v1
ggsave(filename = "plots/venn_loci_PCs_LFMM.png", v1,
       height = 8, width = 8, units = "in")

write.csv(file="landscape_genomics/LFMM_loci_BioClim_PCA_composites_FD0.05.csv",
          sig_loci_unique,
          quote = F, row.names = F)

# ==============================================================================
# 6. Check for overlapping loci 
# ==============================================================================
bioclim_loci <- as.data.frame(read.csv(file="landscape_genomics/LFMM_loci_BioClim_variables_FDR0.05.csv", header = F)[,1])
composite_loci <- as.data.frame(read.csv(file = "landscape_genomics/LFMM_loci_BioClim_PCA_composites_FD0.05.csv", header = F)[,1])

bioclim_loci$method <- "variables"
composite_loci$method <- "composites"

colnames(bioclim_loci)[1]<-colnames(composite_loci)[1]<-"locus"

comb <- as.data.frame(rbind(bioclim_loci, composite_loci))

comb2 <- comb[comb$locus != "Locus",]

venn_list2 <- comb2 %>%
  group_by(method) %>%
  summarise(Loci = list(unique(locus))) %>%
  deframe()

ggvenn(venn_list2,
       fill_color = c("#E69F00", "#56B4E9"),
       stroke_size = 0.3, set_name_size = 4)
