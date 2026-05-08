##################################
### Plotting NGSAdmix outputs ###
#################################

library("ggpubr")
library("tidyverse")
library("viridis")

samples <- read.table(file="all_samples_bam.txt", header = F)
samples2 <- as.data.frame(gsub("_dedup.bam","",gsub(".*\\/","",samples$V1)))
colnames(samples2)<-"Sample"

meta<-read.table("populations.txt")

# ==== Look at the likelihoods of K ====
data<-list.files("popgen/ngsadmix/reps/", pattern = ".log", full.names = T)
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

best_reps <- reps %>%
  group_by(K) %>%
  slice_max(order_by = like, n = 1) %>%
  ungroup()

poporder <- c("Newfoundland", "Labrador", "Quebec_NSLR", "Quebec_SSLR", "New_Brunswick", "Maine")

# ==== Plot K=2-6 for pub ====
## Loop

for(k in 2:6){
 
# Best loglikelihood:
best_run = best_reps[best_reps$K == k,2]

K2 <- read_delim(file = paste0("popgen/ngsadmix/qopt/run_K",k,"_",best_run,".qopt"),delim = " ", col_names = F)
K2[,k+1]<-NULL

K2_admix<-cbind(K2,samples2)
K2_admix<-merge(meta, K2_admix, by.y = "Sample", by.x = "V1")
K2_admix_pal<-viridis(n=k)

K2_admix2<-list()
for(q in 1:(ncol(K2_admix)-2)){
  df <- cbind(K2_admix[c(1,2,q+2)])
  df$K = q
  colnames(df)<-c("Sample","Pop","Prop","K")
  K2_admix2[[q]]<-df
}
K2_admix2<-do.call(rbind, K2_admix2)
K2_admix2$K<-as.factor(K2_admix2$K)
K2_admix2$Pop <- factor(K2_admix2$Pop, levels = poporder)

K2_admix2 <- K2_admix2 %>%
  arrange(Pop, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

p1<-ggbarplot(K2_admix2,
              x = "Sample", y = "Prop", fill = "K",
              color = NA,              # no borders
              palette = K2_admix_pal,        # or choose your own
              x.text.angle = 90,
              xlab = FALSE, ylab = "Ancestry proportion",
              legend.title = "Cluster") +
  facet_grid(~Pop, scales = "free_x", space = "free_x") +
  theme(strip.text = element_text(size = 6, face = "bold"),
        axis.text.x = element_text(size = 6),
        legend.position = "top")

p1

ggsave(paste0("plots/NGSAdmix/NGSAdmix_K",k,".svg"), p1, height = 4, width = 10)
ggsave(paste0("plots/NGSAdmix/NGSAdmix_K",k,".png"), p1, height = 5, width = 10)
}
