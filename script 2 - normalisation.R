library(data.table)
set.seed(1811)

#use all corees available
future::plan(future::multisession,
             workers = 4)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(scran)
# read in QC data
sce <- readRDS('sce_QC.rds')
# now remove the cells marked for discard:
sce <- sce[, !sce$discard]
# annotate spike-in genes:
rowData(sce)$is.spike <- rownames(sce) %like% '^ERCC'
rowData(sce)$is.spike %>% table

# Don't want to accidentally include spike-ins with endogenous genes, so move them to a
# separate slot in the SCE object:
trans.type <- vector(length = nrow(sce))
trans.type[rowData(sce)$is.spike] <- 'ERCC'
trans.type[!rowData(sce)$is.spike] <- 'endogenous'
trans.type <- as.factor(trans.type)
sce <- splitAltExps(sce, trans.type, ref = 'endogenous')

#do a rough clustering with a graph-based method, doing it per batch:
quickclust <- quickCluster(sce, block = sce$batch)
#generate deconvolution based normalisation factors:
deconv <- calculateSumFactors(sce, cluster = quickclust)
summary(deconv)

##calculate size factors based on library size:
libSize <- librarySizeFactors(sce)
# and on spike-in abindance:
spikeSF <- computeSpikeFactors(sce, 'ERCC') %>%
  sizeFactors()

#plot comparisons of all of the methods:
# This is figure 5.2:

png(file = 'Figs/Comparing size factors.png',
    width = 900,
    height = 500)
par(mfrow = c(1, 2))
plot(libSize ~ deconv, xlab = 'SF(deconvolution)', ylab = 'SF(library size)')
abline(0, 1, col = 'red')
plot(spikeSF ~ deconv, xlab = 'SF(deconvolution)', ylab = 'SF(spikes)')
abline(0, 1, col = 'red')
dev.off()

#choose deconvolution method for now:
sizeFactors(sce) <- deconv
# apply normalisation:
sce <- logNormCounts(sce)

#save SCE:
saveRDS(sce, 'sce_QC_norm.rds')
