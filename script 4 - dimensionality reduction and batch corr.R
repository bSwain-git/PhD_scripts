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
sce <- readRDS('sce_QC_norm_HVG.rds')

# run principal components analsis using the HVGs identified earlier:
sce <- runPCA(sce, subset_row = rowData(sce)$HVG)
#plot PCA and colour by phenotype (naive or primed):
plotPCA(sce, colour_by = 'Characteristics.phenotype.', ncomponents = 2)
# and batch:
plotPCA(sce, colour_by = 'batch', ncomponents = 2)
##Choosing PCs:

# by elbow point:
par(mfrow = c(1, 1))
percent.var <- attr(reducedDim(sce), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")

# Factoring in the total technical noise using the total variance of spikes:
denoised <-
  denoisePCA(
    sce,
    technical = modelGeneVarWithSpikes(
      sce,
      spikes = 'ERCC',
      block = sce$Comment.experiment.batch.
    ),
    subset.row = rowData(sce)$HVG
  )
chosen.denoised <- ncol(reducedDim(denoised))

# This is figure 5.4:
cumsum(percent.var) %>% plot(.,
                             main = 'Choosing number of PCs',
                             xlab = 'No. of PCs',
                             ylab = '% variance explained')
abline(v = c(chosen.elbow, chosen.denoised),
       col = c("red", 'blue'))

#let's not use the number of clusters to choose PCs, as there is v. little structure to the data:
# pcs <- reducedDim(sce,'PCA')
# choices <- getClusteredPCs(pcs)
# chosen.clust<-metadata(choices)$chosen

# abline(v=c(chosen.clust), col=c("green"))

legend(
  'bottom',
  legend = c("By elbow point", "By biological variation"),
  col = c("red", "blue"),
  lty = 1
)


#take forward the method that considers the technical variance:
reducedDim(sce, 'PCA_trimmed') <-
  reducedDim(sce, 'PCA')[, c(1:chosen.denoised)]

#looking at UMAP options:
# n_neighbors

sce <- runUMAP(sce, dimred = "PCA_trimmed", n_neighbors = 5)
out5 <- plotReducedDim(sce, dimred = "UMAP",
                       colour_by = "Characteristics.phenotype.") + ggtitle("n_neighbors = 5")


sce <- runUMAP(sce, dimred = "PCA_trimmed", n_neighbors = 20)
out20 <- plotReducedDim(sce, dimred = "UMAP",
                        colour_by = "Characteristics.phenotype.") + ggtitle("n_neighbors = 20")


sce <- runUMAP(sce, dimred = "PCA_trimmed", n_neighbors = 80)
out80 <- plotReducedDim(sce, dimred = "UMAP",
                        colour_by = "Characteristics.phenotype.") + ggtitle("n_neighbors = 80")

multiplot(out5, out20, out80, cols = 3)
#20 seems nicest looking again:


sce <- runUMAP(sce, dimred = "PCA_trimmed", n_neighbors = 20)


#min_dist

sce <-
  runUMAP(sce,
          dimred = "PCA_trimmed",
          min_dist = 0.05,
          n_neighbors = 20)
out005 <- plotReducedDim(sce, dimred = "UMAP",
                         colour_by = "Characteristics.phenotype.") + ggtitle("min_dist = 0.05")

sce <-
  runUMAP(sce,
          dimred = "PCA_trimmed",
          min_dist = .1,
          n_neighbors = 20)
out01 <- plotReducedDim(sce, dimred = "UMAP",
                        colour_by = "Characteristics.phenotype.") + ggtitle("min_dist = 0.1")


sce <-
  runUMAP(sce,
          dimred = "PCA_trimmed",
          min_dist = .5,
          n_neighbors = 20)
out05 <- plotReducedDim(sce, dimred = "UMAP",
                        colour_by = "Characteristics.phenotype.") + ggtitle("min_dist = .5")

multiplot(out005, out01, out05, cols = 3)
#default of 0.1 seems good:

sce <-
  runUMAP(sce,
          dimred = "PCA_trimmed",
          min_dist = .1,
          n_neighbors = 20)


## plot the UMAP coloured by batch and phenotype:
UMAP_pheno <- plotReducedDim(sce, dimred = "UMAP",
                             colour_by = "Characteristics.phenotype.") + ggtitle("Coloured by induction")
UMAP_batch <- plotReducedDim(sce, dimred = "UMAP",
                             colour_by = I(as.factor(sce$batch))) + ggtitle("Coloured by batch")

multiplot(UMAP_pheno, UMAP_batch)

#there is a slight batch effect, as batch 3 is separated robustly. Let's use MNN correction:

sce1 <- sce[, sce$batch == 1]
sce2 <- sce[, sce$batch == 2]
sce3 <- sce[, sce$batch == 3]

library(batchelor)

# First, normalise and equalise coverage between batches:
norm.list <-
  multiBatchNorm(sce1,
                 sce2,
                 sce3,
                 subset.row = rownames(sce)[rowData(sce)$HVG],
                 normalize.all = T)

### use fastMNN to perform PCA-based MNN detection and correction.
sceMNN <- fastMNN(
  norm.list,
  d = chosen.denoised,
  k = 20,
  BSPARAM = BiocSingular::RandomParam(deferred = TRUE),
  correct.all = T,
  subset.row = rownames(sce)[rowData(sce)$HVG]
)
#This produces a corrected PCA matrix:
dim(reducedDim(sceMNN, 'corrected'))

#on which we can run UMAP again:
sceMNN <-
  runUMAP(sceMNN,
          dimred = "corrected",
          min_dist = .1,
          n_neighbors = 20)

sceMNN$batch <- as.factor(sceMNN$batch)

##plotting UMAP of 'corrected' data shows no difference between batches now:
UMAP_batch_MNN <- plotReducedDim(sceMNN, dimred = "UMAP",
                                 colour_by = "batch") + ggtitle("After MNN batch correction")
#How much variance have we lost with this transformation?
totVarLost <- 1 - metadata(sceMNN)$merge.info$lost.var
totVarLost <- (1 - (totVarLost[1, ] * totVarLost[2, ])) * 100

sceMNN$phenotype <- sce$Characteristics.phenotype.
plotUMAP(sceMNN, colour_by = 'batch')
plotUMAP(sceMNN, colour_by = 'phenotype')

# store the MNN corrected data in a slot of the original SCE:
altExp(sce, 'MNN_corrected') <- sceMNN


saveRDS(sce, 'sce_QC_norm_HVG_dimred_corrected.rds')
