# Population Genetics and GEA Analysis in R
library(tidyverse)
library(vcfR)
library(adegenet)
library(LEA)
library(pcadapt)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(RcppCNPy)
library(viridis)
library(ggpubr)

# Set working directory
setwd("./popgen")

# ============================================================================
# ANGSD Results Analysis - PCAngsd & NGSAdmix
# ============================================================================
# Read in the sample and population data and make sure they're in the right order

samples<-read.table(file="../all_samples_bam.txt")
samples$V1<-gsub("/scratch3/workspace/bbentley_smith_edu-lynx/aligned/","",gsub("_dedup.bam","",samples$V1))

pop<-read.table(file="../populations.txt")

sample_data<-merge(samples, pop, by="V1")
colnames(sample_data)<-c("Sample","Population")

# Read the covariance matrix
cov <- as.matrix(read.table(file = "angsd/pca.cov"), header = F)
mme.pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as.data.frame(cbind(sample_data$Sample, sample_data$Population, eigenvectors)) #combine with our population assignments
df = type_convert(pca.vectors)

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
pca.eigenval.sum

varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
varPC4

pal<-viridis(n = length(levels(factor(sample_data$Population))), alpha = 0.5)


colnames(df)[c(1,2)]<-c("Sample","Population")
pcap<-ggscatter(data = df,
          x="V3", y="V4", col = "black", size = 5,
          fill = "Population", palette = pal, shape = 21,
          xlab = paste0("PC1 (",round(varPC1,2),"%)"), ylab = paste0("PC2 (",round(varPC2,2),"%)"),
          title = "ANGSD likelihoods with samtools model" ,
          subtitle = paste0("PCA explains ",round(pca.eigenval.sum,2),"% of the variance"))
pcap
ggsave(filename = "../plots/lynx_ANGSD_PCA_PC1_2.svg", pcap, height = 6, width = 8)

pcap2<-ggscatter(data = df,
                x="V3", y="V5", col = "black", size = 5,
                fill = "Population", palette = pal, shape = 21,
                xlab = paste0("PC1 (",round(varPC1,2),"%)"), ylab = paste0("PC3 (",round(varPC3,2),"%)"),
                title = "ANGSD likelihoods with samtools model" ,
                subtitle = paste0("PCA explains ",round(pca.eigenval.sum,2),"% of the variance"))
pcap2
ggsave(filename = "../plots/lynx_ANGSD_PCA_PC1_3.svg", pcap2, height = 6, width = 8)

### NGSAdmix

# Choose the most appropriate level of K:
data<-list.files("ngsadmix/reps/", pattern = ".log", full.names = T)
klist<-list()
for(q in 1:length(data)){
  df<-read.table(file=data[q], sep = "\t")
  K=as.numeric(gsub("_.*","",gsub(".*_K","",df[1,1])))
  like=as.numeric(gsub("best like=","",gsub(" after.*","",df[11,1])))
  rep=gsub(".*_","",gsub(".*_K","",df[1,1]))
  df2<-as.data.frame(cbind(K,rep,like))
  klist[[q]]<-df2
}
reps<-do.call(rbind,klist)
reps$K<-as.numeric(reps$K)
reps$like<-as.numeric(reps$like)

ggboxplot(data = reps, x = "K", y = "like")


bigData<-lapply(1:120, FUN = function(i) readLines(data[i]))
library(stringr)
foundset<-sapply(1:120, FUN= function(x) bigData[[x]][which(str_sub(bigData[[x]], 1, 1) == 'b')])
as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )
logs<-data.frame(K = rep(c(1,10,11,12,2:9), each=10))
logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))
tapply(logs$like, logs$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

#########################

logs <- as.data.frame(read.table("ngsadmix/logfile"))
logs$K <- c(rep("1", 10), rep("2", 10), rep("3", 10),
            rep("4", 10), rep("5", 10), rep("6", 10))
write.table(logs[, c(2, 1)], "logfile_formatted", row.names = F, 
            col.names = F, quote = F)



admix<-read.table(file="ngsadmix/reps/run_K3_rep1.qopt")
admix<-cbind(admix,sample_data)

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

admix_pal<-viridis(n=3)
#korder<-admix2[c(1:length(levels(factor(admix2$Sample)))),2]
gps<-read.csv(file="../lynx_gps.csv")
gps2<-gps[,c(1:3)]
admix3<-merge(admix2, gps2, by.x = "Sample", by.y = "sample_id")
population_order <- c("Newfoundland", "Labrador", "Quebec_NSLR", "Quebec_SSLR", "New_Brunswick", "Maine")
admix3$Pop <- factor(admix3$Pop, levels = population_order)
admix3 <- admix3[order(admix3$Pop, decreasing = T), ]

ggbarplot(admix3, x = "Sample", y = "Prop", fill = "K", palette = admix_pal, orientation = "horiz") +
  rotate_x_text(90)
ggsave(filename = "../plots/lynx_admix_bar_K4.svg", plot = admix_bp, height = 8, width = 6)


colnames(admix)<-c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Sample", "Population")
admix_gps<-merge(admix, gps2, by.x = "Sample", by.y = "sample_id")
write.table(file="ngsadmix/admix_out_meta.txt", x = admix_gps, quote = F, row.names = F, sep = "\t")

sample_data$Population <- factor(sample_data$Population, levels = population_order)
sample_data2<-sample_data[order(sample_data$Population, decreasing = F), ]

#remotes::install_github('royfrancis/pophelper')
library(pophelper)
slist1<-readQ("ngsadmix/run_K3_rep1.qopt",filetype="basic")
slist1 <- alignK(slist1[1])
labset_order = sample_data[,2,drop=FALSE]
colnames(labset_order)<-"STRATA"
labset_order$STRATA <- as.character(labset_order$STRATA)
plotQ(slist1,  clustercol= admix_pal, grplab = labset_order, grplabsize=3,
      showsp=FALSE, ordergrp=T, imgtype="pdf", exportpath = getwd(),
      showlegend=TRUE, legendpos="right", legendkeysize = 6, legendtextsize = 6,
      legendmargin=c(2,2,2,0), width=20, height=5, sortind="all",
      subsetgrp = c("Newfoundland", "Labrador", "Quebec_NSLR", "Quebec_SSLR", "New_Brunswick", "Maine"),
      outputfilename = "../plots/lynx_admix_pophelper_K3")



# ============================================================================
# GATK Results Analysis
# ============================================================================

# Load VCF data
vcf <- read.vcfR("gatk/popgen_filtered.recode.vcf.gz")
genind <- vcfR2genind(vcf)

# Subset to high coverage samples
high_cov_samples <- read.table("../high_coverage_samples.txt", header = FALSE)$V1
high_cov_samples<-as.data.frame(gsub("/scratch3/workspace/bbentley_smith_edu-lynx/aligned/","",gsub("_dedup.bam","",high_cov_samples)))
colnames(high_cov_samples)<-"Sample"
high_cov_samples2<-merge(high_cov_samples, sample_data, by = "Sample")

genind_hc <- genind[rownames(genind@tab) %in% high_cov_samples, ]

# PCA on hard calls
pca_gatk <- dudi.pca(genind_hc, center = TRUE, scale = FALSE, nf = 10, scannf = FALSE)

pca_gatk_df <- data.frame(
  sample = rownames(pca_gatk$li),
  PC1 = pca_gatk$li[, 1],
  PC2 = pca_gatk$li[, 2],
  PC3 = pca_gatk$li[, 3]
) %>%
  left_join(sample_data, by = "sample")

var_explained_gatk <- pca_gatk$eig / sum(pca_gatk$eig) * 100

p2 <- ggplot(pca_gatk_df, aes(x = PC1, y = PC2, color = population)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = paste0("PC1 (", round(var_explained_gatk[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained_gatk[2], 2), "%)"),
       title = "PCA from GATK (High Coverage Samples)") +
  theme_minimal()

ggsave("popgen/pca_gatk.png", p2, width = 10, height = 8, dpi = 300)

# ============================================================================
# Population Structure Analysis
# ============================================================================

# STRUCTURE-like analysis using LEA
# Convert to lfmm format
vcf2lfmm("gatk/popgen_filtered.recode.vcf.gz", 
         output.file = "popgen/genotypes.lfmm")

# Run snmf for K = 1 to 10
project <- NULL
project <- snmf("popgen/genotypes.lfmm", K = 1:10, 
                entropy = TRUE, repetitions = 5,
                project = "new")

# Plot cross-entropy criterion
plot(project, col = "blue", pch = 19, cex = 1.2)

# Extract Q-matrix for best K
best_k <- which.min(cross.entropy(project, K = 1:10))
qmatrix <- Q(project, K = best_k)

# Plot admixture
structure_df <- data.frame(
  sample = sample_data$sample[sample_data$sample %in% high_cov_samples],
  qmatrix
) %>%
  left_join(sample_data, by = "sample") %>%
  arrange(population)

# Reshape for plotting
structure_long <- structure_df %>%
  select(sample, population, starts_with("X")) %>%
  pivot_longer(cols = starts_with("X"), names_to = "cluster", values_to = "proportion")

p3 <- ggplot(structure_long, aes(x = factor(sample, levels = structure_df$sample), 
                                 y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ population, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = paste("STRUCTURE-like plot (K =", best_k, ")"),
       x = "Samples", y = "Admixture Proportion")

ggsave("popgen/structure_plot.png", p3, width = 15, height = 8, dpi = 300)

# ============================================================================
# Genome-Environment Association (GEA)
# ============================================================================

# Prepare environmental data
env_data <- sample_data %>%
  filter(sample %in% high_cov_samples) %>%
  select(-sample, -population) %>%
  select_if(is.numeric) %>%
  scale()

# Run RDA (Redundancy Analysis)
# Load genetic data as allele frequencies
gen_freq <- scaleGen(genind_hc, NA.method = "mean")

# RDA
rda_result <- rda(gen_freq ~ ., data = as.data.frame(env_data))

# Summary
rda_summary <- summary(rda_result)
print(rda_summary)

# Plot RDA
plot(rda_result, scaling = 3, main = "RDA - Genome-Environment Association")

# Extract candidate SNPs (outliers)
load_rda <- scores(rda_result, choices = c(1:2), display = "species")
hist(load_rda[, 1], main = "Loadings on RDA1")
hist(load_rda[, 2], main = "Loadings on RDA2")

# Define outliers (top 1% of loadings)
cutoff <- quantile(abs(load_rda), 0.99)
outliers <- which(abs(load_rda) >= cutoff, arr.ind = TRUE)
candidate_snps <- unique(rownames(load_rda)[outliers[, 1]])

cat("Number of candidate SNPs:", length(candidate_snps), "\n")

# Save candidate SNPs
write.table(candidate_snps, "gea/candidate_snps.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ============================================================================
# pcadapt analysis for local adaptation
# ============================================================================

# Run pcadapt
bed_file <- "popgen/pruned_snps.bed"
pcadapt_result <- pcadapt(input = bed_file, K = 3)

# Plot results
plot(pcadapt_result, option = "screeplot")
plot(pcadapt_result, option = "scores", pop = sample_data$population)
plot(pcadapt_result, option = "manhattan")
plot(pcadapt_result, option = "qqplot")

# Get p-values and adjust for multiple testing
pvals <- pcadapt_result$pvalues
qvals <- qvalue(pvals)$qvalues

# Identify outliers
alpha <- 0.01
outliers_pcadapt <- which(qvals < alpha)

cat("Number of outlier SNPs (pcadapt):", length(outliers_pcadapt), "\n")

# ============================================================================
# Summary Statistics
# ============================================================================

# Calculate basic population genetics statistics
summary_stats <- list(
  n_samples_total = nrow(sample_data),
  n_samples_high_cov = length(high_cov_samples),
  n_populations = length(unique(sample_data$population)),
  n_snps_total = nrow(vcf@fix),
  n_snps_filtered = nrow(genind_hc@tab),
  candidate_snps_rda = length(candidate_snps),
  outlier_snps_pcadapt = length(outliers_pcadapt)
)

# Save summary
capture.output(summary_stats, file = "analysis_summary.txt")

# Create final summary plot
library(gridExtra)
final_plot <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)
ggsave("final_summary_plot.png", final_plot, width = 20, height = 16, dpi = 300)

cat("Analysis complete! Check the following files:\n")
cat("- popgen/pca_angsd.png: PCA from ANGSD\n")
cat("- popgen/pca_gatk.png: PCA from GATK\n")
cat("- popgen/structure_plot.png: Population structure\n")
cat("- gea/candidate_snps.txt: Candidate SNPs from GEA\n")
cat("- analysis_summary.txt: Summary statistics\n")
cat("- final_summary_plot.png: Combined plots\n")
