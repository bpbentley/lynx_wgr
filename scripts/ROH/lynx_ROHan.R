#####################
### ROHan outputs ###
#####################

library(ggpubr)

sfiles<-list.files(path = "ROHan/100kb/rhomu_1e4", pattern = "summary", full.names = T)

rlist<-list()
for(q in 1:length(sfiles)){
  df<-read.table(sfiles[q], sep = ":", skip = 1)
  df$V2 <- trimws(gsub("\\\\t", " ", df$V2))
  df$V2 <- gsub("\\\t", " ", df$V2)
  
  sample<-gsub("_100kb_autosome_1e4.summary.txt","",basename(sfiles[q]))
  
  # 1. Extract the two middle theta values (line 2)
  theta_noRoh <- df$V2[2]
  theta_noRoh_vals <- as.numeric(unlist(strsplit(theta_noRoh, " ")))
  theta_noRoh_middle <- theta_noRoh_vals[2]  # the two middle values
  
  theta_Roh <- df$V2[3]
  theta_Roh_vals <- as.numeric(unlist(strsplit(theta_Roh, " ")))
  theta_Roh_middle <- theta_Roh_vals[2]  # the two middle values
  
  # 2. Extract the first value in Segments in ROH(%)
  roh_line <- df$V2[7]
  roh_first <- as.numeric(sub(" .*", "", roh_line))
  
  # 3. Extract the Avg. length of ROH
  avg_line <- df$V2[10]
  avg_val <- as.numeric(sub(" .*", "", avg_line))
  
  # Combine results
  results <- data.frame(
    Sample_Name = sample,
    Theta_Outside_ROH = theta_noRoh_middle,
    Theta_Include_ROH = theta_Roh_middle,
    ROH_Percent = roh_first,
    Avg_ROH_Length = avg_val
  )
  rlist[[q]]<-results
  
}

rdf<-do.call(rbind, rlist)

meta<-read.table(file="lynx_clusters_K4.txt")
colnames(meta)<-c("Sample_Name", "Population", "Cluster")

rdf2<-merge(rdf, meta, by = "Sample_Name")
rdf2$log_Avg_ROH_Length <- log(rdf2$Avg_ROH_Length)

rdf2$Cluster<-factor(rdf2$Cluster, ordered = T,
                     levels = c("Cluster1", "Cluster2",
                                "Cluster3", "Cluster4",
                                "Admixed"))

ROH1<-ggboxplot(rdf2, x = "Cluster", y = "ROH_Percent", add = "jitter",
          add.params = list(color = "Population"))
ROH1
ggsave(filename = "plots/ROHan_percROH.svg", ROH1, height = 6, width = 8)

ROH2<-ggboxplot(rdf2, x = "Cluster", y = "Avg_ROH_Length", add = "jitter",
                add.params = list(color = "Population"))
ROH2
ggsave(filename = "plots/ROHan_ROHlen.svg", ROH2, height = 6, width = 8)

ROH3<-ggscatter(rdf2, x = "ROH_Percent", y = "log_Avg_ROH_Length", col = "Cluster")
ROH3

THETA1<-ggboxplot(rdf2, x = "Cluster", y = "Theta_Outside_ROH", add = "jitter",
                  add.params = list(color = "Population"))
THETA1

THETA2<-ggboxplot(rdf2, x = "Cluster", y = "Theta_Include_ROH", add = "jitter",
                  add.params = list(color = "Population"))
THETA2
ggsave(filename = "plots/ROHan_thetaIncROH.svg", THETA2, height = 6, width = 8)


# ==== Individual ROHs ====
alist<-list.files(path = "ROHan/100kb/rhomu_8e4/ind", full.names = T)

for(a in 1:length(alist)){
  df<-read.table(file=alist[a], header = TRUE, comment.char = "")
  
}