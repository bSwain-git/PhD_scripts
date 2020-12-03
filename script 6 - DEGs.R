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

sce <- readRDS('sce_clustered.rds')

##getting markers that are upregulated only - use uncorrected expression,
#and define batch effect. Use wilcox because of the very differnt group sizes. ####

markers <- findMarkers(
  sce,
  sce$clustWalktrap,
  direction = 'up',
  test.type = 'wilcox',
  pval.type = 'any',
  assay.type = 'logcounts',
  block = sce$batch
)

library(AnnotationHub)
## annotate gene symbols ####
ah <- AnnotationHub()
ahDb <- query(ah,
              pattern = c("Homo sapiens", "EnsDb"),
              ignore.case = TRUE)
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)
edb <- ah[[id]]
annotations <- genes(edb,
                     return.type = "data.frame")
annotations <- dplyr::select(annotations, c(gene_id, gene_name))
colnames(annotations) <- c('accession', 'symbol')
rowData(sce) <-
  dplyr::left_join(as.data.frame(rowData(sce)), annotations)


markerSymbols <-
  lapply(markers, function(clusterMarkers0)
    clusterMarkers0$geneSymbol <-
      rowData(sce)$symbol[match(rownames(clusterMarkers0), rowData(sce)$accession)])

for (i in 1:length(markers)) {
  currentClust0 <- markers[[i]]
  currentClust0$symbol <- markerSymbols[[i]]
  markers[[i]] <- currentClust0
}

# annotate sce with marker genes per cluster:
markerInfo <- data.frame('accession' = rownames(sce))
for (annotating0 in 1:length(markers)) {
  currentMarkers <- markers[[annotating0]] %>%
    subset(FDR < 0.05)
  markerInfoCol <- rownames(sce) %in% rownames(currentMarkers) %>%
    data.frame()
  colnames(markerInfoCol) <- paste0('marker_', annotating0)
  markerInfo <- cbind(markerInfo, markerInfoCol)
}
stopifnot(all(rownames(sce) == markerInfo$accession))
rowData(sce) <- cbind(rowData(sce), markerInfo[, -1])

##end####
## get a list of top markers by cluster:####
topMarkerList <- list()
logFCList <- list()

for (i in 1:length(unique(sce$clustWalktrap))) {
  markers_i <- markers[[i]]
  topMarkers_i <- markers_i[markers_i$Top <= 10, ]
  rownames(topMarkers_i) <- topMarkers_i$symbol
  
  topMarkerList[[i]] <- topMarkers_i
  
  logFCs_i <-
    as.matrix(dplyr::select(as.data.frame(topMarkers_i), -c(1:3, 'symbol')))
  colnames(logFCs_i) <- sub("logFC.", "", colnames(logFCs_i))
  
  logFCList[[i]] <- logFCs_i
}
##end####

##getting markers that are downregulated only now####

markersDown <- findMarkers(
  sce,
  sce$clustWalktrap,
  direction = 'down',
  test.type = 'wilcox',
  pval.type = 'any',
  assay.type = 'logcounts',
  block = sce$batch
)

markerSymbolsDown <-
  lapply(markersDown, function(clusterMarkers0)
    clusterMarkers0$geneSymbol <-
      rowData(sce)$symbol[match(rownames(clusterMarkers0), rowData(sce)$accession)])

for (i in 1:length(markersDown)) {
  currentClust0 <- markersDown[[i]]
  currentClust0$symbol <- markerSymbolsDown[[i]]
  markersDown[[i]] <- currentClust0
}

# annotate sce with marker genes per cluster:
markerInfoDown <- data.frame('accession' = rownames(sce))
for (annotating0 in 1:length(markersDown)) {
  currentMarkersDown <- markersDown[[annotating0]] %>%
    subset(FDR < 0.05)
  markerInfoCol <-
    rownames(sce) %in% rownames(currentMarkersDown) %>%
    data.frame()
  colnames(markerInfoCol) <- paste0('marker_', annotating0)
  markerInfoDown <- cbind(markerInfoDown, markerInfoCol)
}
stopifnot(all(rownames(sce) == markerInfoDown$accession))
rowData(sce) <- cbind(rowData(sce), markerInfoDown[, -1])

##end####
## get a list of top markers by cluster:####
topMarkerListDown <- list()
logFCListDown <- list()

for (i in 1:length(unique(sce$clustWalktrap))) {
  markers_i <- markersDown[[i]]
  topMarkers_i <- markers_i[markers_i$Top <= 10, ]
  rownames(topMarkers_i) <- topMarkers_i$symbol
  
  topMarkerListDown[[i]] <- topMarkers_i
  
  logFCs_i <-
    as.matrix(dplyr::select(as.data.frame(topMarkers_i), -c(1:3, 'symbol')))
  colnames(logFCs_i) <- sub("logFC.", "", colnames(logFCs_i))
  
  logFCListDown[[i]] <- logFCs_i
}
##end####


saveRDS(list(markersUp = markers, markersDown = markersDown),
        'markerFiles.rds')
saveRDS(list(logFCListUp = logFCList, logFCListDown = logFCListDown),
        'logFCFiles.rds')
saveRDS(
  list(TopmarkersUp = topMarkerList, TopmarkersDown = topMarkerListDown),
  'TopMarkerFiles.rds'
)


saveRDS(sce, 'sce_final.rds')
# sce <- readRDS('sce_final.rds')

##do markers for naive vs. primed:

markers <- findMarkers(
  sce,
  sce$Characteristics.phenotype.,
  direction = 'up',
  # up regulation
  test.type = 't',
  # have many in each group so no need for non-parametric
  pval.type = 'all',
  # needs to be up or down strictly vs the other groups (kind of loses releavnce when only two groups)
  assay.type = 'logcounts',
  log.p = T,
  # natural log the p-values to maintain accuracy on very small p vals.
  block = sce$batch
)


markerSymbols <-
  lapply(markers, function(clusterMarkers0)
    clusterMarkers0$geneSymbol <-
      rowData(sce)$symbol[match(rownames(clusterMarkers0), rowData(sce)$accession)])

for (i in 1:length(markers)) {
  currentClust0 <- markers[[i]]
  currentClust0$symbol <- markerSymbols[[i]]
  markers[[i]] <- currentClust0
}

saveRDS(markers, 'primedVsNaiveMarkers.rds')

markers.naive <- as.data.frame(markers$naive) %>%
  arrange(log.FDR) %>%
  mutate('FDR' = exp(log.FDR))

markers.primed <- as.data.frame(markers$primed) %>%
  arrange(log.FDR) %>%
  mutate('FDR' = exp(log.FDR))

# markersClust<-readRDS('markerFiles.rds')
markers1 <- markersClust$markersUp$`1` %>%
  as.data.frame() %>%
  arrange(FDR)
