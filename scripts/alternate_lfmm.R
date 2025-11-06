# LFMM from Brenna Forester
# https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html

library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(caret)

gen <- read.lfmm("popgen/input_files/lynx_WGR_pruned_snps.lfmm_imputed.lfmm")
dim(gen)

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
rownames(gen) <- clean_ids
colnames(gen) <- locus_names

env <- read.csv("landscape_genomics/extracted_BioClim_variables_updated_nolatlon.csv", row.names = 1)
common_ids <- intersect(rownames(gen), rownames(env))
env2 <- env[common_ids, ]

# Remove highly correlated predictors
# Visualize correlation structure
cor_matrix <- cor(env2, use = "pairwise.complete.obs")

# Remove highly correlated predictors (|r| > 0.8)
high_corr <- findCorrelation(cor_matrix, cutoff = 0.8, names = TRUE, verbose = FALSE)
high_corr <- high_corr[high_corr != "bio_019"]
high_corr[11]<-"bio_07"
env_filtered <- env2 %>% select(-all_of(high_corr))

# Scale the environmental variables
pred.pca <- rda(env_filtered, scale=T)
summary(pred.pca)$cont

screeplot(pred.pca, main = "Screeplot: Eigenvalues of Lynx Predictor Variables")
round(scores(pred.pca, choices=1:8, display="species", scaling=0), digits=3)


pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)

screeplot(pred.pca, main = "Screeplot of Lynx Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")


gen.pca <- rda(gen, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")

K <- 5

lynx.lfmm <- lfmm_ridge(Y=gen, X=pred.PC1, K=5)
lynx.pv <- lfmm_test(Y=gen, X=pred.PC1, lfmm=lynx.lfmm, calibrate="gif")
lynx.pv$gif

lynx.qv <- qvalue(lynx.pv$calibrated.pvalue)$qvalues

length(which(lynx.qv < 0.001))
