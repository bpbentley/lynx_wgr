############################
### OutFLANK - Lynx WGR ####
### 2025-06-24 B Bentley ###
####### From Vignette ######
############################
#devtools::install_github("whitlock/OutFLANK")


library(devtools)
library(qvalue)
library(vcfR)
library(adegenet)
library(poppr)
library(tidyverse)
library(hierfstat)
library(OutFLANK)

# ==== Read in VCF and convert to OutFLANK ====
vcf <- read.vcfR("popgen/input_files/lynx_WGR_pruned_snps.vcf")

geno <- extract.gt(vcf) # Character matrix containing the genotypes
position <- getPOS(vcf) # Positions in bp
chromosome <- getCHROM(vcf) # Chromosome information

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

table(as.vector(G))
G[is.na(G)] <- 9

# ==== Load popmap and extract ====
popmap<-read.table(file="popmap.txt", header = F, col.names = c("sample","pop"))
ids<-colnames(vcf@gt)[c(2:ncol(vcf@gt))]
clean_ids <- ifelse(
  grepl("^(.*)_\\1$", ids),       # Match repeated pattern like A202_A202
  sub("_(.*)$", "", ids),         # Strip suffix
  ids                             # Otherwise keep as-is
)

popmap2<-popmap[order(popmap$sample),]
pop<-popmap2$pop
id<-popmap2$sample

my_fst <- MakeDiploidFSTMat(t(G), locusNames = position, popNames = pop)
head(my_fst)
plot(my_fst$He, my_fst$FST)

plot(my_fst$FST, my_fst$FSTNoCorr)
abline(0,1)

# Already LD-trimmed with PLINK upstream so no further trimming needed
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
head(P1)
tail(P1)

my_out <- P1$OutlierFlag==TRUE
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")
