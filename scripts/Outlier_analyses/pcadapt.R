#############################
### PCAdapt for Lynx WGR  ###
#############################

#==== Load libraries & change directory ====

library(pcadapt)
library(qvalue)
library(ggpubr)
library(viridis)
library(vcfR)

setwd("popgen")

# ==== Run PCAdapt ====
# Note - make sure to convert the VCF to PLINK format ahead of this step.

# NB requires all PLINK outputs in the same directory (.bed .bim .fam)
path_to_file <- paste0("input_files/", "lynx_WGR_pruned_snps.bed")
filename <- read.pcadapt(path_to_file, type = "bed")

# Run the function with a large K-value initially
x <- pcadapt(input = filename, K = 20)

# Check the outputs
plot(x, option = "screeplot")

# Bottoms out at ~K=8, plot the first 10 values of K.
p1<-plot(x, option = "screeplot", K = 10) + theme_classic()
p1
ggsave("../Plots/PCAdapt_scree.png", p1, height = 6, width = 8, units = "in")

# Based on Scree Plot, K=3 with the elbow method
# Re-run with the selected K value (in this case, K=3).
x <- pcadapt(filename, K = 3)
summary(x)

# ==== Visualize PCAdapt outputs ====

# Plot the outliers as a Manhattan plot
p2<-plot(x , option = "manhattan") + theme_classic()
p2
ggsave("../Plots/PCAdapt_Manhattan.png", p2, height = 6, width = 8, units = "in")

# Check the QQ-plot:
plot(x, option = "qqplot") + theme_classic()

# Check the p-value distribution:
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# Check the statistic distribution (not informative here).
plot(x, option = "stat.distribution")

# Correct the p-value for FDR using the q-value function.
qval <- qvalue(x$pvalues)$qvalues

# Extract the outliers and determine how many are present
## Adjust alpha value as needed.
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers) #1,847 SNPs

# ==== Add the population information ====

popmap<-read.table(file="popmap.txt", header = F, col.names = c("sample","pop"))
fam<-read.table(file="input_files/lynx_WGR_pruned_snps.fam", header = F)
fam$ID<-paste0(fam$V1,"_", fam$V2)

# Issues with vcf2plink conversion corrupting sample names
## This command reduced the sample names back to original.
fam$clean_ids <- ifelse(
  grepl("^(.*)_\\1$", fam$ID),       # Match repeated pattern like A202_A202
  sub("_(.*)$", "", fam$ID),         # Strip suffix
  fam$ID                             # Otherwise keep as-is
)

# Order the samples so they match the data:
poporder<-fam$clean_ids
popmap2<-popmap[match(poporder, popmap$sample),]

# Plot a PCA based on all SNPs
pca_out<-as.data.frame(x$scores)
pca_out$sample<-popmap2$sample
pca_out$pop<-popmap2$pop

pal<-viridis(n = 6, alpha = 0.7)
p3<-ggscatter(pca_out, x = "V1", y = "V2", fill = "pop", shape = 21, size = 5, palette = pal,
          xlab = paste0("PC1 (",round((x$singular.values[1]^2)*100,2),"%)"),
          ylab = paste0("PC2 (",round((x$singular.values[2]^2)*100,2),"%)"),
          legend.title = "Collection site")
p3
ggsave("../Plots/PCAdapt_PCA.png", p3, height = 6, width = 8, units = "in")

# Check in matches the PCAdapt outs:
plot(x, option = "scores", pop = popmap2$pop)

# ==== Generate a matrix of locus IDs which can be removed from the neutral dataset ====
bim<-read.table(file="input_files/lynx_WGR_pruned_snps.bim", header = F)
all_loci<-bim$V2
outlier_loci<-all_loci[outliers]

write.table(file="../landscape_genomics/PCAdapt_outlier_loci_names_K3_FDR0.05.txt", outlier_loci, quote = F, row.names = F,
            col.names = F, sep = "\t")

# ==== Plot the outlier loci on the PCA ====
loads<-as.data.frame(x$loadings[,c(1,2)])
rownames(loads)<-all_loci

loads$outlier <- ifelse(rownames(loads) %in% outlier_loci, "outlier", "non_significant")

outpal<-c("black", "red")
p4<-ggscatter(data = loads, x = "V1", y = "V2", col = "outlier", palette = outpal,
          xlab = "PC1", ylab = "PC2")
p4
ggsave("../Plots/PCAdapt_Outlier_PCA_K3_0.05.png", p4, height = 6, width = 8, units = "in")
