############################################################
# Cluster 29 Proliferation Signature Analysis
# Description:
#   - Filters significant genes (p_val_adj < 0.05)
#   - Computes proliferation score
#   - Visualizes in UMAP 2D and 3D
############################################################

library(Seurat)
library(dplyr)
library(plotly)
library(ggplot2)


cluster_of_interest <- "29"
pval_threshold <- 0.05

############################################################
# DEFINE HUMAN PROLIFERATION GENE SET (UPPERCASE)
############################################################

prolif_genes <- c(
  "BUB1","BUB1B","CCNB2","CCNA2","CDKN3","CDCA3","PLK1","FOXM1","TOP2A",
  "CENPH","CENPW","CENPK","NUF2","SKA3","KNSTRN","NCAPG","NCAPD2",
  "KIF15","KIF2C","HMMR","CEP55","UBE2C",
  "MKI67","BIRC5","SHCBP1","TPX2","RRM2","UHRF1"
)

############################################################
# FILTER SIGNIFICANT GENES FROM YOUR MARKER TABLE
# Assumes you already have a data.frame called `markers_all`
# with columns: cluster, gene, p_val_adj
############################################################

sig_genes <- markers_all %>%
  filter(cluster == cluster_of_interest,
         p_val_adj < pval_threshold,
         gene %in% prolif_genes) %>%
  pull(gene)

cat("Significant proliferation genes in cluster 29:\n")
print(sig_genes)

############################################################
# CHECK GENE PRESENCE IN SEURAT OBJECT
############################################################

sig_genes_present <- sig_genes[sig_genes %in% rownames(data)]

cat("Genes present in Seurat object:\n")
print(sig_genes_present)

############################################################
# COMPUTE MEAN EXPRESSION SCORE PER CELL
############################################################

expr_df <- FetchData(data, vars = sig_genes_present)
expr_df$cluster <- data$seurat_clusters

expr_df$ProlifScore <- rowMeans(expr_df[, sig_genes_present])

data$ProlifScore_cluster29 <- expr_df$ProlifScore

############################################################
# SUMMARY STATISTICS
############################################################

summary_table <- data@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_score = mean(ProlifScore_cluster29, na.rm = TRUE)) %>%
  arrange(desc(mean_score))

print(summary_table)

############################################################
# UMAP 2D VISUALIZATION
############################################################

p2d <- FeaturePlot(
  data,
  features = "ProlifScore_cluster29",
  reduction = "umap"
)

print(p2d)

############################################################
# UMAP 3D VISUALIZATION 
############################################################

umap_coords <- Embeddings(data, "umap")
umap_df <- as.data.frame(umap_coords)
umap_df$ProlifScore <- data$ProlifScore_cluster29
umap_df$Cluster <- data$seurat_clusters

p3d <- plot_ly(
  data = umap_df,
  x = ~umap_1,
  y = ~umap_2,
  z = ~umap_3,
  type = "scatter3d",
  mode = "markers",
  color = ~ProlifScore,
  marker = list(size = 3),
  text = ~paste("Cluster:", Cluster,
                "<br>Score:", round(ProlifScore, 2))
) %>%
  layout(title = "UMAP 3D - Cluster 29 Proliferation Signature")

p3d

############################################################
# VIOLIN PLOT BY CLUSTER
############################################################

p_violin <- VlnPlot(
  data,
  features = "ProlifScore_cluster29",
  group.by = "seurat_clusters"
)

print(p_violin)