##################################################
### Combine PCAdapt with LFMM for neutral loci ###
##################################################

library(ggvenn)
library(dplyr)
library(tibble)
library(vcfR)

pcadapt <- read.csv(file="landscape_genomics/PCAdapt_outlier_loci_names_FD0.05.txt", 
                    header = F, col.names = "Locus")
pcadapt$Method <- "PCAdapt"


lfmm_bioclim <- as.data.frame(read.table(file="landscape_genomics/LFMM_loci_BioClim_variables_FDR0.05.csv",
                         header = T))
colnames(lfmm_bioclim) <- "Locus"
lfmm_bioclim$Locus <- gsub(",.*","",lfmm_bioclim$Locus)
lfmm_bioclim$Method <- "LFMM"

rda <- as.data.frame(read.csv(file="landscape_genomics/RDA_loci_just_climate.csv",
                              header = T)[,2])
colnames(rda) <- "Locus"
rda$Method <- "RDA"

comb_methods <- as.data.frame(rbind(pcadapt, lfmm_bioclim, rda))
venn_list <- comb_methods %>%
  group_by(Method) %>%
  summarise(Loci = list(unique(Locus))) %>%
  deframe()


v1<-ggvenn(venn_list,
       fill_color = c("#E69F00", "#56B4E9", "dark red"),
       stroke_size = 0.3, set_name_size = 4)
v1

ggsave(filename = "plots/venn_lfmm_pcadapt_rda.png", v1,
       height = 6, width = 6, units = "in")

locus_counts <- table(comb_methods$Locus)
overlapping_loci <- names(locus_counts[locus_counts >= 2])

# ==== Extract the outlier loci from the VCF ====
vcf <- read.vcfR("popgen/input_files/lynx_WGR_pruned_snps.vcf")
vcf_clean <- vcf[!(vcf@fix[,"ID"] %in% overlapping_loci), ]
write.vcf(vcf_clean, file = "popgen/input_files/lynx_WGR_neutral_snps.vcf")

vcf_outliers <- vcf[(vcf@fix[,"ID"] %in% overlapping_loci), ]
write.vcf(vcf_outliers, file = "popgen/input_files/lynx_WGR_adaptive_snps.vcf")

# ==== Draw Manhattan's I guess? ====
vcf_info <- as.data.frame(getFIX(vcf))[, c("CHROM", "POS", "ID")]
colnames(vcf_info) <- c("CHR", "BP", "SNP")


vcf <- read.vcfR("popgen/input_files/lynx_WGR_adaptive_snps.vcf")
