############################################
### Population genomics of Canadian lynx ###
####### Whole-genome resequence data #######
############# Neutral SNPs only ############
############################################

## Input data: VCF generated through bwa-mem alignments & GATK hard calls
## VCF Filters: (see GitHub repo) -- depth, mapping quality, MAF, LD, etc.
## Input file contains only neutral SNPs, outliers removed
## Outliers detected through: (1) PCAdapt; (2) LFMM
## Author: Blair Bentley

# ==== Load required libraries ====

library(hierfstat)
library(mapdata)
library(viridis)
library(tidyr)
library(adegenet)
library(vcfR)
library(ggplot2)
library(poppr)
library(ggpubr)
library(viridis)
library(reshape2)
library(ape)
library(dartR)
library(strataG)
library(StAMPP)
library(LEA)

# ==== Set working directory to "popgen" subdirectory ====
setwd("popgen")

# ==== Load genotype and population data ====
vcf <- vcfR::read.vcfR("input_files/lynx_WGR_neutral_snps.vcf")
pop <- read.table("../popgen/popmap.txt", header=F, sep="\t", stringsAsFactors = TRUE,
                  col.names = c("Indiv","STRATA"))
# ==== Convert to needed file, rename samples, add pop data ====
genind <- vcfR::vcfR2genind(vcf)
genl<-vcfR2genlight(vcf)

# The pipeline duplicated SOME of the sample names in the VCF, so need to address that:
indNames(genind)<- ifelse(
  grepl("^(.*)_\\1$", indNames(genind)),       # Match repeated pattern like A202_A202
  sub("_(.*)$", "", indNames(genind)),         # Strip suffix
  indNames(genind)                             # Otherwise keep as-is
)

indNames(genl)<- ifelse(
  grepl("^(.*)_\\1$", indNames(genl)),       # Match repeated pattern like A202_A202
  sub("_(.*)$", "", indNames(genl)),         # Strip suffix
  indNames(genind)                             # Otherwise keep as-is
)

# Need to re-order the popmap to make sure it's the same as the genotypes:
sample_order <- indNames(genind)
pop2<-pop[match(sample_order, pop$Indiv),]


# Now add the data to the genotypes:
genind@pop <- pop2$STRATA
genl@pop <- pop2$STRATA

# Run some summary metrics:
tab1<-as.data.frame(table(genind@pop)) # Print the total number of samples used in the analyses per population.
tab1

# Get heterozygosity
#hf <- genind2hierfstat(genind)
#ho_individual1 <- 1 - rowSums(hf[, -1] == 0 | hf[, -1] == 2, na.rm=TRUE) / rowSums(!is.na(hf[, -1]))
#ho_individual2 <- summary(genind)$Hobs

# Generate a palette for plotting:
pop_pal <- viridis(n=length(levels(factor(pop2$STRATA))), alpha = 0.8)
poporder<-c("Newfoundland", "Labrador", "Quebec_NSLR",
            "Quebec_SSLR", "New_Brunswick", "Maine")

# -----------------------------------------------------------------------------

# ==== Principal Component Analysis (PCA) ====
# Run first with nf = 20 to determine how many PCs to retain
pca1 <- glPca(genl,center = T, scale = T, nf = 4) # Took ~30 mins to run on 90,000 SNPs
barplot(100*pca1$eig/sum(pca1$eig), col = heat.colors(50), main="PCA Eigenvalues") # retain first 4 axes, incremental decrease after 2
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#proportion of explained variance by first four axes
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
a4<-pca1$eig[4]/sum(pca1$eig) # proportion of variation explained by 3rd axis
pcvar <- data.frame(Axis = c(1:4), Proportion = c(a1,a2,a3,a4))
pcvar

#Extract PC scores to color by location
pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$pop <- pop(genl)
pca1.scores$ind <- genl@ind.names

# Plot the PCA
set.seed(9)
num_pops <- length(levels(factor(pop2$pop)))

# plot PC 1 and 2
pca1.p<-ggscatter(pca1.scores, x = "PC1", y = "PC2", shape = 21, fill = "pop",
                  palette = pop_pal, size = 5, legend.title = "Collection site",
                  title = "PCA: Canadian lynx - Neutral loci",
                  subtitle = paste0(genl@n.loc," SNPs"),
                  xlab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)"),
                  ylab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"))
pca1.p

# #plot  PC 2 and 3: not informative for this dataset
pca2.p<-ggscatter(pca1.scores, x = "PC2", y = "PC3", shape = "Coast", color = "pop",
                  palette = pop_pal, ellipse = F, ellipse.level = 0.95,
                  xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  ylab = paste0("PC3 (",round(pcvar[3,2]*100,2),"%)"))
pca2.p

ggsave("../plots/lynx_neutral_PCA_FINAL.svg",pca1.p,width=8,height=6,dpi=600,units="in")
# ------------------------------------------------------------------------------

# ==== Discriminant Analysis of PCA (DAPC) ====
n_individuals<-nrow(pop2)
n_pops<-length(levels(factor(pop2$STRATA)))

grp <- find.clusters(genind, max.n.clust = 12) # Choose the number of clusters
table(pop(genind), grp$grp)

grp_all <- find.clusters(genind, max.n.clust=12, n.pca=40,
                         choose.n.clust = FALSE)
BIC<-as.data.frame(cbind(seq(1,12,1), grp_all$Kstat))
p1<-ggline(BIC, x = "V1", y = "V2", plot_type = "b",
       col = "navy",
       xlab = "Number of clusters (K)",
       ylab = "BIC Value",
       title = "Selection of optimum number of clusters (K; DAPC)") +
  font("xlab", face = "bold") + font("ylab", face = "bold")
p1
ggsave("../Plots/lynx_DAPC_K_selection.png", p1, height = 6, width = 8, units = "in", dpi = 300)

grp_all$Kstat

grp <- find.clusters(genind, n.clust = 3, choose.n.clust = F, n.pca = 40) # Choose the number of clusters
table(pop(genind), grp$grp)

png(filename = "../plots/lynx_DAPC_cluster_assign.png", width = 8, height = 6, units = "in", res = 300)
table.value(table(pop(genind), grp$grp), col.labels = paste0("Cluster ", 1:3))
dev.off()

# DAPC
dapc <- dapc(genind, grp$grp, n.pca=40, n.da = 3, var.contrib = TRUE)
dapc_pal <- magma(n = 5)[c(1,3,4)]
dpca_result <- scatter(dapc, col=dapc_pal, scree.pca = TRUE, posi.da = "topleft",
                       posi.pca = "topright",
                       pch = 21, cell = 0, cstar = 1,
                       solid = 1, cex = 3, clab = 1)
set.seed(4)
compoplot(dapc, posi="topleft", txt.leg=paste("Cluster",1:3),
          xlab="Individuals", col=dapc_pal, lab=genind@pop, show.lab = T)

# Run secondary DAPC with first as prior
set.seed(5)
dapc_a_score <- dapc(genind,n.pca=40,n.da=3)

png(filename = "../plots/lynx_DAPC_ascore.png", width = 8, height = 6, units = "in", res = 300)
temp_score <- optim.a.score(dapc_a_score)
dev.off()

dapc2 <-dapc(genind, grp$grp, n.pca=6, n.da=3)
png(filename = "../plots/lynx_DAPC_scatter.png", width = 8, height = 6, units = "in", res = 300)
dpca_result <- scatter(dapc2, col=dapc_pal, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1, posi.da = "topright",
                       posi.pca = "bottomright")
dev.off()
load_dpca2 <- as.data.frame(dapc2$var.contr)
percent= dapc2$eig/sum(dapc2$eig)*100
dapc_prior=as.data.frame(dapc2$ind.coord)
dapc_prior$IND <- row.names(dapc_prior)
dapc_prior$SITE <- pop2[match(dapc_prior$IND, pop2$Indiv),2]

dapc1<-ggscatter(dapc_prior, x = "LD1", y = "LD2", shape = 21, fill = "SITE", size = 5, xlab = paste0("DA1 (", round(percent[1],2),"%)"),
                 ylab = paste0("DA2 (", round(percent[2],2),"%)"), palette = pop_pal, alpha = 0.67)
dapc1
ggsave("../plots/lynx_DAPC_scatt_ggpubr.png", dapc1, height = 6, width = 8, units = "in", dpi = 300)

# ==== Structure/Admixture ====
geno <- ped2geno(input.file = "input_files/lynx_WGR_pruned_snps.ped")
geno <- paste0("input_files/lynx_WGR_pruned_snps.geno")
project <- snmf(geno, K = 1:6, repetitions = 100, entropy = TRUE, alpha = 100, project = "new")
project <- load.snmfProject("./input_files/lynx_WGR_pruned_snps.snmfProject")

K_values <- 2:6
ce_means <- sapply(K_values, function(k) {
  mean(cross.entropy(project, K = k))
})

ce_sds <- sapply(K_values, function(k) {
  sd(cross.entropy(project, K = k))
})

ce_mean_df<-as.data.frame(cbind(K_values, ce_means, ce_sds))

ggline(ce_mean_df, x = "K_values", y = "ce_means",
       xlab = "K-value", ylab = "Cross-entropy mean (±SD ;100 reps)") +
  geom_errorbar(aes(ymin = ce_means - ce_sds, ymax = ce_means + ce_sds, group = K_values),
                width = 0.1, position = position_dodge(0.2)) +
  font("xlab", face = "bold") + font("ylab", face = "bold")

for(k in 2:5){
#best_K <- K_values[which.min(ce_means)]
best_K = k
cat("Best K =", best_K, "\n")

# Get the run with lowest cross-entropy for best K
best_run <- which.min(cross.entropy(project, K = k))
#best_run = 14

cat("Best run =", best_run, "\n")
Q_best <- Q(project, K = best_K, run = best_run)

# Plot barplot of ancestry proportions
barplot(t(Q_best), col = rainbow(best_K), border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions", main = paste("sNMF: K =", best_K))
# Again, this suggests that K=3 but will run with K=2 through K=6

admix<-read.table(file=paste0("./input_files/lynx_WGR_pruned_snps.snmf/K",k,"/run",best_run,"/lynx_WGR_pruned_snps_r",best_run,".",k,".Q"))

sample_order <- read.table(file="./input_files/lynx_WGR_pruned_snps.fam")[,c(1,2)]
sample_order <- paste0(sample_order[,1],"_", sample_order[,2])
sample_order<- ifelse(
  grepl("^(.*)_\\1$", sample_order),       # Match repeated pattern like A202_A202
  sub("_(.*)$", "", sample_order),         # Strip suffix
  sample_order                             # Otherwise keep as-is
)

admix2 <- as.data.frame(cbind(sample_order, admix))
colnames(admix2)[c(2:ncol(admix2))] <- paste0("Cluster",rep(1:k))
meta <- read.table(file = "../populations.txt")
admix3 <- as.data.frame(merge(admix2, meta, by.x = "sample_order", by.y = "V1"))
colnames(admix3)[k+2] <- "Pop"
admix3 <- admix3[match(sample_order, admix3$sample_order), ]

admix4<-list()
for(q in 1:(ncol(admix3)-2)){
  df<-as.data.frame(admix3[,c(q+1,1,ncol(admix3))])
  df$K<-q
  colnames(df)<-paste0("V",seq(1:4))
  admix4[[q]]<-df
}
admix4<-do.call(rbind, admix4)
colnames(admix4)<-c("Prop","Sample","Pop","K")
admix4$K<-as.factor(admix4$K)
admix4$Pop <- factor(admix4$Pop, levels = poporder)

admix4 <- admix4 %>%
  arrange(Pop, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

admix_pal<-viridis(n=length(levels(admix4$K)))

p1<-ggbarplot(admix4,
          x = "Sample", y = "Prop", fill = "K",
          color = NA,              # no borders
          palette = admix_pal,        # or choose your own
          x.text.angle = 90,
          xlab = FALSE, ylab = "Ancestry proportion",
          legend.title = "Cluster") +
  facet_grid(~Pop, scales = "free_x", space = "free_x") +
  theme(strip.text = element_text(size = 6, face = "bold"),
        axis.text.x = element_text(size = 6),
        legend.position = "top")
p1
ggsave(filename = paste0("../plots/Admixture_snmf/lynx_WGR_sNMF_K",k,".png"),
       p1, height = 4, width = 10, units = "in", dpi = 300)
}
