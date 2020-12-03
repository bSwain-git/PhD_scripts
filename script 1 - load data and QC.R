## dataset comparing the expression in H9 hESCs
## cultured to select naive or primed phenotypes.

library(data.table)
set.seed(1811)
#phenotypic data from database
pheno <- readRDS("data/pheno_E-MTAB-6819.rds")
#use all corees available
future::plan(future::multisession,
             workers = 4)

# read in list of raw data
filesList <- list.files('data/', pattern = '.tsv')
# read in first counts matrix
counts <- read.table(paste0('data/', filesList[1]), header = T)

# append all counts matrices to form a large overall counts matrix:
for (i in 2:length(filesList)) {
  fileI <- read.table(paste0('data/', filesList[i]), header = T)
  counts <- dplyr::inner_join(counts, fileI)
  
}

###

library(SingleCellExperiment)
library(tidyverse)
library(data.table)
# make a unique identifier in the counts matrix that matches phenotypic data:
pheno$sampleID <-
  make.unique(str_sub(pheno$Comment.SUBMITTED_FILE_NAME., end = -13))

# remove the extraneous columns in counts matrix:
dim(counts)
countsFiltered <-
  cbind(counts$GeneID, counts[, colnames(counts) %in% pheno$sampleID]) %>%
  column_to_rownames('counts$GeneID')
dim(countsFiltered)

# remove duplicate phenotypic info:
pheno0 <-
  subset(pheno, !duplicated(dplyr::select(pheno, Assay.Name)))

# check that all cells in counts matric have a phenotypic entry:
stopifnot(sum(pheno0$Assay.Name %in% colnames(countsFiltered)) == ncol(countsFiltered))

# collate the info into a single object of SingleCellExperiment class:
sce <-
  SingleCellExperiment(assays = list(counts = as.matrix(countsFiltered)))
colData(sce) <- DataFrame(pheno0)
colnames(sce) <- pheno0$sampleID
rowData(sce)$accession <- rownames(sce)
sce$batch <- sce$Comment.experiment.batch.


## QC steps:
dim(sce)
# remove genes that were undetected in all cells:
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(AnnotationHub)
# annotate mitochondrial genes from accession numbers:
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
annotations <-
  dplyr::select(annotations, c(gene_id, gene_name, seq_name))
mitochondrialGenes <-
  annotations$gene_id[annotations$seq_name %like% 'MT']



library(scater)
# which genes are mitochondrial:
is.mito0 <- which(rownames(sce) %in% mitochondrialGenes)
# which are spike-in transcripts
is.ERCC0 <- which(rownames(sce) %like% '^ERCC-')
# Calculate QC metrics per cell, split by batch:
df <- perCellQCMetrics(sce, subsets = list('Mito' = is.mito0,
                                           'ERCCspike' = is.ERCC0))
# Detect library size outliers:
qc.lib <-
  isOutlier(df$sum,
            log = TRUE,
            type = "lower",
            batch = sce$batch)
# Features per cell outliers:
qc.nexprs <-
  isOutlier(df$detected,
            log = TRUE,
            type = "lower",
            batch = sce$batch)
# Outliers by proportion of spike counts:
qc.spikes <-
  isOutlier(df$subsets_ERCCspike_percent,
            type = "higher",
            batch = sce$batch)
# And finally by proportion of mitochondrial counts:
qc.mito.batch <-
  isOutlier(df$subsets_Mito_percent,
            type = "higher",
            batch = sce$batch)
# Print htresholds for reference:
attr(qc.lib, 'thresholds')
attr(qc.nexprs, 'thresholds')
attr(qc.spikes, 'thresholds')
attr(qc.mito.batch, 'thresholds')
# discard cells that fail any step:
sce$discard <- qc.lib | qc.nexprs | qc.mito.batch | qc.spikes
# add QC info to SCE object:
colData(sce) <- cbind(colData(sce), df)
# also make a data.frame of QC data for easy plotting:
forQCplots <- as.data.frame(df)
forQCplots$batch <- factor(sce$batch)
# more intuitive to label cells that are kept rather than discarded:
forQCplots$keep <- !sce$discard
# This is figure 5.1B:
## features plot ####
featurePlot <-
  ggplot(forQCplots,
         aes(x = batch,
             y = detected,
             fill = keep)) +
  geom_violin(aes(group = batch), scale = 'width', colour = 'gray') +
  geom_point(
    position = position_jitter(width = 0.05),
    colour = 'black',
    shape = 21,
    alpha = 0.3
  ) +
  xlab('Batch') +
  ylab('Detected genes') +
  ggtitle(paste0("Detected features per cell")) +
  geom_segment(aes(
    y = attr(qc.nexprs, 'thresholds')['lower', 1],
    yend = attr(qc.nexprs, 'thresholds')['lower', 1],
    x = 0.75,
    xend = 1.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  geom_segment(aes(
    y = attr(qc.nexprs, 'thresholds')['lower', 2],
    yend = attr(qc.nexprs, 'thresholds')['lower', 2],
    x = 1.75,
    xend = 2.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  geom_segment(aes(
    y = attr(qc.nexprs, 'thresholds')['lower', 3],
    yend = attr(qc.nexprs, 'thresholds')['lower', 3],
    x = 2.75,
    xend = 3.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  ylim(0, 16000) +
  labs(fill = 'Pass\n(overall)') +
  guides(colour = F) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks = element_line(),
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.2))
  )
featurePlot
##end####
# This is figure 5.1A:
## libSize ####
libPlot <-
  ggplot(forQCplots,
         aes(x = batch,
             y = sum,
             fill = keep)) +
  geom_violin(aes(group = batch), scale = 'width', colour = 'gray') +
  geom_point(
    position = position_jitter(width = 0.05),
    colour = 'black',
    shape = 21,
    alpha = 0.3
  ) +
  xlab('Batch') +
  ylab('Library size') +
  ggtitle(paste0("Library size per cell")) +
  geom_segment(aes(
    y = attr(qc.lib, 'thresholds')['lower', 1],
    yend = attr(qc.lib, 'thresholds')['lower', 1],
    x = 0.75,
    xend = 1.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  geom_segment(aes(
    y = attr(qc.lib, 'thresholds')['lower', 2],
    yend = attr(qc.lib, 'thresholds')['lower', 2],
    x = 1.75,
    xend = 2.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  geom_segment(aes(
    y = attr(qc.lib, 'thresholds')['lower', 3],
    yend = attr(qc.lib, 'thresholds')['lower', 3],
    x = 2.75,
    xend = 3.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  # ylim(0,16000)+
  scale_y_log10() +
  labs(fill = 'Pass\n(overall)') +
  guides(colour = F) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks = element_line(),
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.2))
  )
libPlot

##end ####
# This is figure 5.1C:
## mito proportion ####
mitoPlot <-
  ggplot(forQCplots,
         aes(x = batch,
             y = subsets_Mito_percent,
             fill = keep)) +
  geom_violin(aes(group = batch), scale = 'width', colour = 'gray') +
  geom_point(
    position = position_jitter(width = 0.05),
    colour = 'black',
    shape = 21,
    alpha = 0.3
  ) +
  xlab('Batch') +
  ylab('Mitochondrial counts (%)') +
  ggtitle(paste0("Proportion of mitochondrial counts per cell")) +
  geom_segment(aes(
    y = attr(qc.mito.batch, 'thresholds')['higher', 1],
    yend = attr(qc.mito.batch, 'thresholds')['higher', 1],
    x = 0.75,
    xend = 1.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  geom_segment(aes(
    y = attr(qc.mito.batch, 'thresholds')['higher', 2],
    yend = attr(qc.mito.batch, 'thresholds')['higher', 2],
    x = 1.75,
    xend = 2.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  geom_segment(aes(
    y = attr(qc.mito.batch, 'thresholds')['higher', 3],
    yend = attr(qc.mito.batch, 'thresholds')['higher', 3],
    x = 2.75,
    xend = 3.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  # ylim(0,100)+
  # scale_y_log10()+
  labs(fill = 'Pass\n(overall)') +
  guides(colour = F) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks = element_line(),
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.2))
  )
mitoPlot

##end ####
# This is figure 5.1D:
## spikes proportion ####
spikesPlot <-
  ggplot(forQCplots,
         aes(x = batch,
             y = subsets_ERCCspike_percent,
             fill = keep)) +
  geom_violin(aes(group = batch), scale = 'width', colour = 'gray') +
  geom_point(
    position = position_jitter(width = 0.05),
    colour = 'black',
    shape = 21,
    alpha = 0.3
  ) +
  xlab('Batch') +
  ylab('Spike-in counts (%)') +
  ggtitle(paste0("Proportion of spike-in counts per cell")) +
  geom_segment(aes(
    y = attr(qc.spikes, 'thresholds')['higher', 1],
    yend = attr(qc.spikes, 'thresholds')['higher', 1],
    x = 0.75,
    xend = 1.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  geom_segment(aes(
    y = attr(qc.spikes, 'thresholds')['higher', 2],
    yend = attr(qc.spikes, 'thresholds')['higher', 2],
    x = 1.75,
    xend = 2.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  geom_segment(aes(
    y = attr(qc.spikes, 'thresholds')['higher', 3],
    yend = attr(qc.spikes, 'thresholds')['higher', 3],
    x = 2.75,
    xend = 3.25
  ),
  col = 'black',
  # linetype = 'dashed',
  size = rel(1.5)) +
  # ylim(0,100)+
  # scale_y_log10()+
  labs(fill = 'Pass\n(overall)') +
  guides(colour = F) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks = element_line(),
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.2))
  )
spikesPlot

##end ####

## Are we discarding a specific cell type that is erroneously being annotated as poor-quality?
# calculate average counts of every gene in discarded cells:
lost <- calculateAverage(counts(sce)[, sce$discard])
# Do the same for kept cells:
kept <- calculateAverage(counts(sce)[,!sce$discard])

# Calculated counts per million for a more easily comparible set of gene counts, and log2:
library(edgeR)
logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 1)
#Calculate fold-change of genes in lost cells vs kept cells (already logged cpm):
logFC <- logged[, 1] - logged[, 2]
#What is the overall abundance of each gene in the un-pruned dataset?
abundance <- rowMeans(logged) %>%
  data.frame('abundance' = .) %>%
  rownames_to_column('gene_id')
#Extract names of spike genes:
spikeGenes <- rownames(sce)[which(rownames(sce) %like% '^ERCC')]
# Label mitochondrial genes in this table:
abundance$mito <- abundance$gene_id %in% mitochondrialGenes
# and spike genes:
abundance$spike <- abundance$gene_id %in% spikeGenes
# Add logFC between lost and kept cells per gene:
abundance$logFC <- logFC
# annotate this data.frame
abundance <-
  left_join(abundance, annotations[, c('gene_id', 'gene_name')])

# which genes are over-discarded, above an arbitrary threshold of 2? (See graph)
over.discarded <- (logFC)[(logFC) > 2] %>%
  data.frame('logFC' = .) %>%
  rownames_to_column('gene_id')
over.discarded <- inner_join(over.discarded, annotations)
# Extract the over-discarded genes that aren't spike-ins or mitochondrial:
table.genes <- abundance %>%
  subset(mito == F & spike == F & logFC > 2) %>%
  dplyr::select(gene_name)
# Save for external gene enrichment analysis:
write.table(table.genes, 'over_discarded_genes.csv')

library(ggrepel)
library(gridExtra)
#split table for better visualisation:
table1 <- tableGrob(table.genes[1:17,], rows = NULL)
table2 <- tableGrob(table.genes[18:33,], rows = NULL)

# scatter plot of logFC vs abundance, colouring mitochondrial and spike-in genes:
plot <- ggplot(data = abundance, aes(x = abundance, y = logFC)) +
  geom_point(aes(alpha = logFC > 2)) +
  geom_point(data = subset(abundance, mito), aes(color = 'Mitochondrial genes')) +
  geom_point(data = subset(abundance, spike), aes(color = 'Spike-in controls')) +
  theme_classic() +
  xlab('Average log2(CPM + 1)') +
  ylab('Log-FC(lost/kept') +
  guides('alpha' = F) +
  labs('colour' = NULL) +
  scale_color_manual(values = list(
    'Mitochondrial genes' = 'red',
    'Spike-in controls' = 'blue'
  )) +
  scale_alpha_manual(values = c(0.3, 1)) +
  geom_hline(
    yintercept = 2,
    colour = 'red',
    size = 2,
    linetype = 'dashed'
  ) +
  # ylim(-3,6)+
  theme(axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)))


# arrange the scatter plot and tables in a single figure:
# This is figure 5.1E:
grid.arrange(plot,
             table1,
             table2,
             ncol = 3,
             widths = c(.8, 0.1, .1))

saveRDS(sce, 'sce_QC.rds')
