######################
### Het from ANGSD ###
######################

library(ggpubr)
library(viridis)

pop_pal <- viridis(n=6, alpha = 1)
poporder<-c("Newfoundland", "Labrador", "Quebec_NSLR",
            "Quebec_SSLR", "New_Brunswick", "Maine")

all_files<-list.files(path = "pop_stats/ANGSD_het", full.names = T)

het_list<-list()
for(q in 1:length(all_files)){
  df<-scan(all_files[q])
  sample = gsub(".est.ml","",gsub(".*/","",all_files[q]))
  het_prop = df[2]/sum(df)
  het_perc = het_prop * 100
  df2<-as.data.frame(cbind(sample, het_prop, het_perc))
  het_list[[q]]<-df2
}
het_df<-do.call(rbind, het_list)
het_df[,c(2:3)] <- lapply(het_df[,c(2:3)], as.numeric)

meta <- read.table(file="lynx_clusters_K4.txt")

het_meta <- merge(het_df, meta, by.x = "sample", by.y = "V1")
colnames(het_meta)[c(4,5)]<-c("population", "cluster")
dep<-read.csv(file="metrics/all_depth.csv")
het_meta_dep<-merge(het_meta, dep, by.x = "sample", by.y = "Sample_ID", all = T)
het_meta_dep<-het_meta_dep[!is.na(het_meta_dep$het_perc),]

ggboxplot(data = het_meta, x = "population", y = "het_perc", add = "jitter") +
  rotate_x_text(45)

het_meta$cluster <- factor(het_meta$cluster,
                           levels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed"),
                           ordered = T)
angsd_het<-ggboxplot(data = het_meta, x = "cluster", y = "het_perc", add = "jitter",
          add.params = list(color = "population", size = 3), palette = pop_pal,
          xlab = "Cluster assigned by NGSAdmix (K=4)",
          ylab = "Heterozygosity (%; ANGSD)",
          legend.title = "Sampling Location") +
  rotate_x_text(45) + font("xlab", face = "bold") + font("ylab", face = "bold")
angsd_het
ggsave("plots/lynx_Het_ANGSD.svg", angsd_het, height = 6, width = 8)

ggscatter(data = het_meta_dep, x = "Average_depth", y = "het_perc")
