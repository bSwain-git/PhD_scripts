library(data.table)
set.seed(1811)

#use all corees available
future::plan(future::multisession,
             workers = 4)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(scran)
library(scater)

sce <- readRDS('sce_QC_norm_HVG_dimred_corrected.rds')

##repeated for both corrected and uncorrected values:


SNNgraph10 <-
  buildSNNGraph(altExp(sce, 'MNN_corrected'),
                use.dimred = 'corrected',
                k = 10)
SNNgraph20 <-
  buildSNNGraph(altExp(sce, 'MNN_corrected'),
                use.dimred = 'corrected',
                k = 20)
SNNgraph30 <-
  buildSNNGraph(altExp(sce, 'MNN_corrected'),
                use.dimred = 'corrected',
                k = 30)

library(igraph)
clust10 <- cluster_walktrap(SNNgraph10)$membership
clust20 <- cluster_walktrap(SNNgraph20)$membership
clust30 <- cluster_walktrap(SNNgraph30)$membership

#assess the separation of clusters; higher modularity means better sep.
ratios10 <- clusterModularity(SNNgraph10, clust10, as.ratio = T)
ratios20 <- clusterModularity(SNNgraph20, clust20, as.ratio = T)
ratios30 <- clusterModularity(SNNgraph30, clust30, as.ratio = T)

#most higher values (observed/expected weights) should therefore fall on the diagonal
#(same cluster.)
library(pheatmap)
hm10 <- pheatmap(
  log2(ratios10 + 1),
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(c("white", "blue"))(100),
  main = 'SNN k=10'
)
hm20 <- pheatmap(
  log2(ratios20 + 1),
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(c("white", "blue"))(100),
  main = 'SNN k=20'
)
hm30 <- pheatmap(
  log2(ratios30 + 1),
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(c("white", "blue"))(100),
  main = 'SNN k=30'
)

library(gridExtra)
grid.arrange(arrangeGrob(grobs = list(hm10[[4]], hm20[[4]], hm30[[4]]), ncol =
                           3))

##adjacency graphs
cluster.gr10 <-
  graph_from_adjacency_matrix(ratios10,
                              mode = "upper",
                              weighted = TRUE,
                              diag = FALSE)
cluster.gr20 <-
  graph_from_adjacency_matrix(ratios20,
                              mode = "upper",
                              weighted = TRUE,
                              diag = FALSE)
cluster.gr30 <-
  graph_from_adjacency_matrix(ratios30,
                              mode = "upper",
                              weighted = TRUE,
                              diag = FALSE)

par(mfrow = c(1, 3))

plot(cluster.gr10,
     edge.width = E(cluster.gr10)$weight * 2,
     main = 'k=10')
plot(cluster.gr20,
     edge.width = E(cluster.gr20)$weight * 2,
     main = 'k=20')
plot(cluster.gr30,
     edge.width = E(cluster.gr30)$weight * 2,
     main = 'k=30')


sce$clust10 <- clust10 %>%
  as.factor()
sce$clust20 <- clust20 %>%
  as.factor()
sce$clust30 <- clust30 %>%
  as.factor()

altExp(sce, 'MNN_corrected')$clust10 <- clust10 %>%
  as.factor()
altExp(sce, 'MNN_corrected')$clust20 <- clust20 %>%
  as.factor()
altExp(sce, 'MNN_corrected')$clust30 <- clust30 %>%
  as.factor()

#plot the clusters on a dimRed plot (UMAP):
clustPlot10 <-
  plotUMAP(altExp(sce, 'MNN_corrected'), colour_by = 'clust10')
clustPlot20 <-
  plotUMAP(altExp(sce, 'MNN_corrected'), colour_by = 'clust20')
clustPlot30 <-
  plotUMAP(altExp(sce, 'MNN_corrected'), colour_by = 'clust30')

multiplot(clustPlot10, clustPlot20, clustPlot30)

#k = 30 seems best:
sce$clustWalktrap <- sce$clust30
sce$clust10 <- NULL
sce$clust20 <- NULL
sce$clust30 <- NULL

altExp(sce, 'MNN_corrected')$clustWalktrap <- sce$clust30
altExp(sce, 'MNN_corrected')$clust10 <- NULL
altExp(sce, 'MNN_corrected')$clust20 <- NULL
altExp(sce, 'MNN_corrected')$clust30 <- NULL


saveRDS(sce, 'sce_clustered.rds')
