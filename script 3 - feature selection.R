library(data.table)
set.seed(1811)
#use all corees available
future::plan(future::multisession,
             workers = 4)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(scran)
sce <- readRDS('sce_QC_norm.rds')

#use the spike-ins to model the technical noise:
gene.var <-
  modelGeneVarWithSpikes(sce,
                         parametric = T,
                         spikes = 'ERCC',
                         block = sce$batch)

#per batch, plot the mean variance trend for all genes, highlighting the spike-ins;
# plot a trendline using the spike-ins.
# These form figure 5.3:
par(mfrow = c(1, 3))
blocked.stats <- gene.var$per.block
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(
    current$mean,
    current$total,
    main = paste('Batch', i),
    pch = 16,
    cex = 0.5,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
  )
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col = "red", pch = 16)
  curve(curfit$trend(x),
        col = 'dodgerblue',
        add = TRUE,
        lwd = 2)
}

# how many genes are 'significantly' (loose definition) high variance?
length(gene.var$FDR[gene.var$FDR < 0.05])
# Let's arbitrarily use the top 2000 genes for our highly variable subset:
HVG <- getTopHVGs(gene.var, n = 2000)
rowData(sce)$HVG <- rownames(rowData(sce)) %in% HVG

saveRDS(sce, 'sce_QC_norm_HVG.rds')
