###################
#### Lynx GONE2 ###
###################

library(ggpubr)
library(dplyr)

# ==== Pull cluster assignment from NGSadmix ====
### To be used for subsetting VCF
ngs <- read.csv(file="popgen/ngsadmix/lynx_K4_props.csv")
ngs2 <- ngs[,c(1:5)]

### Classify into clusters. For the purposes on Ne, if >30% secondary ancestry
### consider the sample to be "admixed" and removed from this analysis

df_classified <- ngs2 %>%
  rowwise() %>%
  mutate(
    maxK = which.max(c_across(Cluster1:Cluster4)),  # cluster with highest proportion
    maxProp = max(c_across(Cluster1:Cluster4)),     # value of that cluster
    secondProp = sort(c_across(Cluster1:Cluster4), decreasing = TRUE)[2],  # second-highest
    assignment = case_when(
      secondProp >= 0.3 ~ "Admixed",                      # admixed if >30% in another cluster
      TRUE ~ paste0("Cluster", maxK)                      # otherwise assign to main cluster
    )
  ) %>%
  ungroup() %>%
  select(sample_order2, Cluster1:Cluster4, assignment)

meta <- read.table(file="populations.txt")
df_classified<-merge(df_classified, meta, by.x = "sample_order2", by.y = "V1")

### Write to file (for cluster so tsv):
# Cluster 1 (i.e. Newfoundland, N=3)
cluster1<-df_classified[df_classified$assignment == "Cluster1",]
write.table(file = "popgen/ngsadmix/cluster1_assignment_K4.tsv", cluster1, quote = F, row.names = F, sep = "\t")

# Cluster 2 (Quebec NSLR, Labrador, 1x New Brunswick, N=18) New Brunswick not in VCF
cluster2<-df_classified[df_classified$assignment == "Cluster2",]
write.table(file = "popgen/ngsadmix/cluster2_assignment_K4.tsv", cluster2, quote = F, row.names = F, sep = "\t")

# Cluster 3 (Maine; N=16)
cluster3<-df_classified[df_classified$assignment == "Cluster3",]
write.table(file = "popgen/ngsadmix/cluster3_assignment_K4.tsv", cluster3, quote = F, row.names = F, sep = "\t")

# Cluster 4 (Quebec SSLR, 2x New Brunswick, 1x Maine, N=10)
cluster4<-df_classified[df_classified$assignment == "Cluster4",]
write.table(file = "popgen/ngsadmix/cluster4_assignment_K4.tsv", cluster4, quote = F, row.names = F, sep = "\t")

# Admixed (Maine + 1x Quebec SSLR; N=10)
admixed<-df_classified[df_classified$assignment == "Admixed",]
write.table(file = "popgen/ngsadmix/admixed_assignment_K4.tsv", admixed, quote = F, row.names = F, sep = "\t")

# ==== Use the above to subset the VCF ====

Ne <- read.table(file="GONE2/cluster2_GONE2.ped_GONE2_Ne", header = T)

Ne$logNe <- log10(Ne$Ne_diploids)

ggline(data = Ne, x = "Generation", y = "logNe", plot_type = "l")
ggline(data = Ne, x = "Generation", y = "Ne_diploids", plot_type = "l")
