# ==============================================================================
# GProfiler 2 for lynx adaptive genes
# ==============================================================================

library(gprofiler2)
library(ggpubr)
library(dplyr)
library(viridis)

# Pull gene list for snpEff output
genes<-read.table(file="landscape_genomics/snpEff/lynx_adaptive.genes.txt", header = T, sep = "\t")

# Summarize to unique genes
gene_list<-levels(factor(genes$GeneId)) #379 genes

# Run gProfiler against humand gene DB
gprofiler_out <- gost(query = gene_list, organism = "hsapiens",
                      user_threshold = 0.05, correction_method = "fdr",
                      evcodes = TRUE,)

gost<-gostplot(gprofiler_out, capped = FALSE, interactive = F)
gost
ggsave(filename = "plots/lynx_gProfiler_gost.svg", gost)

gprofiler_results <- gprofiler_out$result

#write.csv(file = "landscape_genomics/snpEff/adaptive_genes_lynx_Aug2025.csv", gene_list,
#          quote = F, row.names = F)

# ==============================================================================
# Biological Process
# ==============================================================================

bp_terms <-  gprofiler_results[ gprofiler_results$source == "GO:BP",]

#write.table(file = "landscape_genomics/gProfiler2/BP_enrichment.txt", bp_terms2,
#            quote = F, row.names = F, sep = "\t")

top_bp_terms <- head(bp_terms[order(bp_terms$p_value, decreasing = FALSE), ], 20)

top_bp_terms$term_name <- factor(top_bp_terms$term_name,
                                 levels = top_bp_terms$term_name[order(top_bp_terms$p_value)])
order_vec <- top_bp_terms$term_name[order(top_bp_terms$p_value, decreasing = TRUE)]
order_vec2 <- top_bp_terms$term_name[order(top_bp_terms$intersection_size, decreasing = F)]

p1<-ggbarplot(data = top_bp_terms,
          x = "term_name",
          y = "intersection_size",
          fill = "p_value",
          order = order_vec,
          ylab = "Number of genes",
          xlab = "GO term name",
          legend.title = "p-value") +
  coord_flip() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  font("xlab", face = "bold") + font("ylab", face = "bold") +
  theme(legend.text = element_text(angle = 45, hjust = 1))
p1

ggsave(file="plots/lynx_snpEff_adaptive_BP_Pvalue.svg", p1, height = 8, width = 10)

p2<-ggbarplot(data = top_bp_terms,
              x = "term_name",
              y = "intersection_size",
              fill = "p_value",
              order = order_vec2,
              ylab = "Number of genes",
              xlab = "GO term name",
              legend.title = "p-value") +
  coord_flip() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  font("xlab", face = "bold") + font("ylab", face = "bold") +
  theme(legend.text = element_text(angle = 45, hjust = 1))
p2

ggsave(file="plots/lynx_snpEff_adaptive_BP_intersect.svg", p1, height = 8, width = 10)

# ==============================================================================
# Cellular Component
# ==============================================================================
cc_terms <-  gprofiler_results[ gprofiler_results$source == "GO:CC",]

#write.table(file = "landscape_genomics/gProfiler2/cc_enrichment.txt", cc_terms2,
#            quote = F, row.names = F, sep = "\t")

top_cc_terms <- head(cc_terms[order(cc_terms$p_value, decreasing = FALSE), ], 20)

top_cc_terms$term_name <- factor(top_cc_terms$term_name,
                                 levels = top_cc_terms$term_name[order(top_cc_terms$p_value)])
order_vec <- top_cc_terms$term_name[order(top_cc_terms$p_value, decreasing = TRUE)]
order_vec2 <- top_cc_terms$term_name[order(top_cc_terms$intersection_size, decreasing = F)]

p3<-ggbarplot(data = top_cc_terms,
              x = "term_name",
              y = "intersection_size",
              fill = "p_value",
              order = order_vec,
              ylab = "Number of genes",
              xlab = "GO term name",
              legend.title = "p-value") +
  coord_flip() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  font("xlab", face = "bold") + font("ylab", face = "bold") +
  theme(legend.text = element_text(angle = 45, hjust = 1))
p3

ggsave(file="plots/lynx_snpEff_adaptive_CC_Pvalue.svg", p3, height = 8, width = 10)

# ==============================================================================
# Molecular Function
# ==============================================================================
mf_terms <-  gprofiler_results[ gprofiler_results$source == "GO:MF",]

#write.table(file = "landscape_genomics/gProfiler2/mf_enrichment.txt", mf_terms2,
#            quote = F, row.names = F, sep = "\t")

top_mf_terms <- head(mf_terms[order(mf_terms$p_value, decreasing = FALSE), ], 20)

top_mf_terms$term_name <- factor(top_mf_terms$term_name,
                                 levels = top_mf_terms$term_name[order(top_mf_terms$p_value)])
order_vec <- top_mf_terms$term_name[order(top_mf_terms$p_value, decreasing = TRUE)]
order_vec2 <- top_mf_terms$term_name[order(top_mf_terms$intersection_size, decreasing = F)]

p4<-ggbarplot(data = top_mf_terms,
              x = "term_name",
              y = "intersection_size",
              fill = "p_value",
              order = order_vec,
              ylab = "Number of genes",
              xlab = "GO term name",
              legend.title = "p-value") +
  coord_flip() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  font("xlab", face = "bold") + font("ylab", face = "bold") +
  theme(legend.text = element_text(angle = 45, hjust = 1))
p4

ggsave(file="plots/lynx_snpEff_adaptive_MF_Pvalue.svg", p4, height = 8, width = 10)

# ==============================================================================
# Write the list to file:
# ==============================================================================
write.table(file="landscape_genomics/adaptive_gene_list.txt", gene_list,
            sep = "\t", quote = F, row.names = F, col.names = F)
