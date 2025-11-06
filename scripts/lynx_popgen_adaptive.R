############################################
### Population genomics of Canadian lynx ###
####### Whole-genome resequence data #######
############ Adaptive SNPs only ############
############################################

## Input data: VCF generated through bwa-mem alignments & GATK hard calls
## VCF Filters: (see GitHub repo) -- depth, mapping quality, MAF, LD, etc.
## Input file contains only outliers detected through: (1) PCAdapt; (2) LFMM
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
vcf <- vcfR::read.vcfR("input_files/lynx_WGR_adaptive_snps.vcf")
pop <- read.table("../popgen/popmap.txt", header=F, sep="\t", stringsAsFactors = TRUE,
                  col.names = c("Indiv","STRATA"))
# ==== Convert to needed file, rename samples, add pop data ====
genind <- vcfR::vcfR2genind(vcf)
genl<-vcfR2genlight(vcf)

# The pipeline duplicated sample names in the VCF, so need to address that:
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
                  title = "PCA: Canadian lynx - Adaptive loci",
                  subtitle = paste0(genl@n.loc," SNPs"),
                  xlab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)"),
                  ylab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"))
pca1.p

ggsave("../Plots/adaptive_pca.png",pca1.p,width=8,height=8,units="in")

# ==== Discriminant Analysis of PCA (DAPC) ====
n_individuals<-nrow(pop2)
n_pops<-length(levels(factor(pop2$STRATA)))

grp <- find.clusters(genind, max.n.clust = 12) # Choose the number of clusters
table(pop(genind), grp$grp)

grp_all <- find.clusters(genind, max.n.clust=12, n.pca=200,
                         choose.n.clust = FALSE)
BIC<-as.data.frame(cbind(seq(1,12,1), grp_all$Kstat))
ggline(BIC, x = "V1", y = "V2", plot_type = "b",
       col = "navy",
       xlab = "Number of clusters (K)",
       ylab = "BIC Value",
       title = "Selection of optimum number of clusters (K)") +
  font("xlab", face = "bold") + font("ylab", face = "bold")
grp_all$Kstat

grp <- find.clusters(genind, n.clust = 5, choose.n.clust = F, n.pca = 200) # Choose the number of clusters
table(pop(genind), grp$grp)

table.value(table(pop(genind), grp$grp), col.labels = paste0("Cluster ", 1:3))

# DAPC
dapc <- dapc(genind, grp$grp, n.pca=50, n.da = 5, var.contrib = TRUE)
dapc_pal <- magma(n = 5)
dpca_result <- scatter(dapc, col=dapc_pal, scree.pca = TRUE, posi.da = "topleft",
                       posi.pca = "topright",
                       pch = 21, cell = 0, cstar = 1,
                       solid = 1, cex = 3, clab = 1)
set.seed(4)
compoplot(dapc, posi="topleft", txt.leg=paste("Cluster",1:5),
          xlab="Individuals", col=dapc_pal, lab=genind@pop, show.lab = T)

# Run secondary DAPC with first as prior
set.seed(5)
dapc_a_score <- dapc(genind,n.pca=50,n.da=5)
temp_score <- optim.a.score(dapc_a_score)

dapc2 <-dapc(genind, grp$grp, n.pca=2, n.da=5)
dpca_result <- scatter(dapc2, col=dapc_pal, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1, posi.da = "topright",
                       posi.pca = "bottomright")
load_dpca2 <- as.data.frame(dapc2$var.contr)
percent= dapc2$eig/sum(dapc2$eig)*100
dapc_prior=as.data.frame(dapc2$ind.coord)
dapc_prior$IND <- row.names(dapc_prior)
dapc_prior$SITE <- pop2[match(dapc_prior$IND, pop2$Indiv),2]

dapc1<-ggscatter(dapc_prior, x = "LD1", y = "LD2", shape = 21, fill = "SITE", size = 5, xlab = paste0("DA1 (", round(percent[1],2),")%"),
                 ylab = paste0("DA2 (", round(percent[2],2),")%"), palette = pop_pal, alpha = 0.67)
dapc1

# ==== Structure/Admixture ====
geno <- ped2geno(input.file = "input_files/lynx_WGR_forAdmix.ped")
project <- snmf(geno, K = 1:15, repetitions = 5, entropy = TRUE, project = "new")

K_values <- 1:15
ce_means <- sapply(K_values, function(k) {
  mean(cross.entropy(project, K = k))
})

ce_sds <- sapply(K_values, function(k) {
  sd(cross.entropy(project, K = k))
})

ce_mean_df<-as.data.frame(cbind(K_values, ce_means, ce_sds))

ggline(ce_mean_df, x = "K_values", y = "ce_means",
       xlab = "K-value", ylab = "Cross-entropy mean (5 reps)")

best_K <- K_values[which.min(ce_means)]
cat("Best K =", best_K, "\n")

# Get the run with lowest cross-entropy for best K
best_run <- which.min(cross.entropy(project, K = best_K))
Q_best <- Q(project, K = best_K, run = best_run)

# Plot barplot of ancestry proportions
barplot(t(Q_best), col = rainbow(best_K), border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions", main = paste("sNMF: K =", best_K))
# Again, this suggests that K=2 but will run with K=2 through K=6

# Below pulled from Uro work with Admix

admix<-read.table(file="input_files/lynx_WGR_forAdmix.snmf/K6/run1/lynx_WGR_forAdmix_r1.6.Q")
admix<-cbind(admix,pop2)

admix2<-list()
for(q in 1:(ncol(admix)-2)){
  df<-as.data.frame(admix[,c(q,(ncol(admix)-1),ncol(admix))])
  df$K<-q
  colnames(df)<-c("V1","V2","V3","V4")
  admix2[[q]]<-df
}
admix2<-do.call(rbind, admix2)
colnames(admix2)<-c("Prop","Sample","Pop","K")
admix2$K<-as.factor(admix2$K)
admix2$Pop <- factor(admix2$Pop, levels = poporder)

admix2 <- admix2 %>%
  arrange(Pop, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

admix_pal<-viridis(n=length(levels(admix2$K)))

ggbarplot(admix2,
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
