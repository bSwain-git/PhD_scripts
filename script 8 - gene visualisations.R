sce <- readRDS('sce_final.rds')
library(tidyverse)
library(scater)

corrected <- altExp(sce, 'MNN_corrected')

accessions <- list()
accessions[['ABCG2']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'ABCG2']
accessions[['SLC9A1']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'SLC9A1']
accessions[['SLC4A4']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'SLC4A4']
accessions[['SLC9A3']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'SLC9A3']
accessions[['SLC4A2']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'SLC4A2']
accessions[['GATA3']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'GATA3']

accessions[['TERT']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'TERT']
accessions[['KLF4']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'KLF4']
accessions[['BMP4']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'BMP4']
accessions[['PDPN']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'PDPN']
accessions[['THY1']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'THY1']
accessions[['CEACAM1']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'CEACAM1']

UMAPforPlot <- data.frame(
  'cluster' = sce$clustWalktrap,
  'UMAP1' = reducedDim(corrected, 'UMAP')[, 1],
  'Uncorr1' = reducedDim(sce, 'UMAP')[, 1],
  'UMAP2' = reducedDim(corrected, 'UMAP')[, 2],
  'Uncorr2' = reducedDim(sce, 'UMAP')[, 2],
  'batch' = sce$batch,
  'Cellstate' = sce$Characteristics.phenotype.,
  'ABCG2' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$ABCG2,],
  'ABCG2.logcounts' = assay(sce, 'logcounts')[rownames(corrected) %in% accessions$ABCG2,]
  ,
  'SLC9A1' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$SLC9A1,]
  ,
  'SLC4A4' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$SLC4A4,]
  ,
  'SLC9A3' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$SLC9A3,],
  'TERT.logcounts' = assay(sce, 'logcounts')[rownames(sce) %in% accessions$TERT,],
  'KLF4.logcounts' = assay(sce, 'logcounts')[rownames(sce) %in% accessions$KLF4,],
  'GATA3' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$GATA3,],
  'SLC4A2' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$SLC4A2,],
  'BMP4.logcounts' = assay(sce, 'logcounts')[rownames(sce) %in% accessions$BMP4,],
  'PDPN' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$PDPN,],
  'THY1' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$THY1,],
  'CEACAM1' = assay(corrected, 'reconstructed')[rownames(corrected) %in% accessions$CEACAM1,]
  
  
  
)

# UMAP plots for cell state/ gene expression ####
# This is figure 5.6A:

ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = Cellstate)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  annotate(
    geom = 'label',
    x = 3.5,
    y = 5,
    label = 'Naive',
    colour = '#F8766D',
    size = 8
  ) +
  annotate(
    geom = 'label',
    x = 0,
    y = -2.5,
    label = 'Primed',
    colour = '#00BFC4',
    size = 8
  ) +
  guides(colour = F)
#
#ABCG2
# This is figure 5.6H:

ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = ABCG2)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('ABCG2'), 'expression'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )
#
# ##NHE1
# This is figure 5.6F:
ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = SLC9A1)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('SLC9A1'), 'expression'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )


#SLC4A2
# This is figure 5.6G:

ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = SLC4A2)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('SLC4A2'), 'expression'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )

#GATA3
# This is figure 5.6I:
ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = GATA3)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('GATA3'), 'expression'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )



#PDPN - should be in naive and primed
# This is figure 5.6C:

ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = PDPN)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('PDPN'), 'expression'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )

#THY1 - should be in  primed only
# This is figure 5.6D:
ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = THY1)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('THY1'), 'expression'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )



#CEACAM1 - should be in  naive only
# This is figure 5.6E:
ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = CEACAM1)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('CEACAM1'), 'expression'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )

#BMP4 logcounts
# This is figure 5.8D:
ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = BMP4.logcounts)) +
  geom_point(alpha = 1) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('BMP4'), 'logcounts'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )

#TERT logcounts
# This is figure 5.8B:
ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = TERT.logcounts)) +
  geom_point(alpha = 1) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('TERT'), 'logcounts'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )
#ABCG2 logcounts
# This is figure 5.8A:
ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = ABCG2.logcounts)) +
  geom_point(alpha = 1) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('ABCG2'), 'logcounts'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )

#KLF4 logcounts
# This is figure 5.8C:

ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = KLF4.logcounts)) +
  geom_point(alpha = 1) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = expression(atop(italic('KLF4'), 'logcounts'))) +
  guides(colour = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  # guides(colour=F)+
  scale_colour_viridis_c(option = 'B') +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.8, 0.2)
  )


# end UMAPs ####


## for violin plots and cluster UMAP ####
library(ggpubr)
# by cluster
# This is figure 5.6B:
ggplot(UMAPforPlot, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
  geom_point(alpha = 0.5) +
  theme_void() +
  guides(colour = F)

#order the clusters more intuitively for violin plots:
UMAPforPlot$clusterForOrd <-
  factor(UMAPforPlot$cluster, levels = c('1', '2', '5', '3', '4', '6'))

# ABCG2 (5.6H)
ggplot(UMAPforPlot, aes(x = clusterForOrd, y = ABCG2, fill = cluster)) +
  geom_violin(width = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15),
    alpha = 0.1,
    shape = 21,
    aes(fill = cluster)
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('\nCluster') +
  ylab('ABCG2 expression\n') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), vjust = -1),
    axis.title.x = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.line = element_line()
  )

#PDPN (5.6C)
ggplot(UMAPforPlot, aes(x = clusterForOrd, y = PDPN, fill = cluster)) +
  geom_violin(width = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15),
    alpha = 0.1,
    shape = 21,
    aes(fill = cluster)
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('\nCluster') +
  ylab('PDPN expression\n') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), vjust = -1),
    axis.title.x = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.line = element_line()
  )

# THY1 (5.6D)

ggplot(UMAPforPlot, aes(x = clusterForOrd, y = THY1, fill = cluster)) +
  geom_violin(width = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15),
    alpha = 0.1,
    shape = 21,
    aes(fill = cluster)
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('\nCluster') +
  ylab('THY1 expression\n') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), vjust = -1),
    axis.title.x = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.line = element_line()
  )

# CEACAM1 (5.6E)
ggplot(UMAPforPlot, aes(x = clusterForOrd, y = CEACAM1, fill = cluster)) +
  geom_violin(width = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15),
    alpha = 0.1,
    shape = 21,
    aes(fill = cluster)
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('\nCluster') +
  ylab('CEACAM1 expression\n') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), vjust = -1),
    axis.title.x = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.line = element_line()
  )

# SLC9A1 (5.6F)
ggplot(UMAPforPlot, aes(x = clusterForOrd, y = SLC9A1, fill = cluster)) +
  geom_violin(width = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15),
    alpha = 0.1,
    shape = 21,
    aes(fill = cluster)
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('\nCluster') +
  ylab('SLC9A1 expression\n') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), vjust = -1),
    axis.title.x = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.line = element_line()
  )

# SLC4A2 (5.6G)
ggplot(UMAPforPlot, aes(x = clusterForOrd, y = SLC4A2, fill = cluster)) +
  geom_violin(width = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15),
    alpha = 0.1,
    shape = 21,
    aes(fill = cluster)
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('\nCluster') +
  ylab('SLC4A2 expression\n') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), vjust = -1),
    axis.title.x = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.line = element_line()
  )

# GATA3 (5.6I)
ggplot(UMAPforPlot, aes(x = clusterForOrd, y = GATA3, fill = cluster)) +
  geom_violin(width = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15),
    alpha = 0.1,
    shape = 21,
    aes(fill = cluster)
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('\nCluster') +
  ylab('GATA3 expression\n') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), vjust = -1),
    axis.title.x = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.line = element_line()
  )

# uncorrected abcg2 logcounts (5.8A)
ggplot(UMAPforPlot,
       aes(
         x = clusterForOrd,
         y = ABCG2.logcounts,
         colour = as.factor(batch)
       )) +
  geom_boxplot(width = 0.5) +
  geom_point(
    position = position_dodge(width = 0.5),
    alpha = 0.1,
    shape = 21,
    aes(fill = as.factor(batch))
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('Cluster') +
  ylab('ABCG2 logcounts (uncorrected)') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), margin = margin(t = 5)),
    axis.title.x = element_text(size = rel(1.5), margin = margin(t =
                                                                   15)),
    axis.title.y = element_text(
      size = rel(1.5),
      angle = 90,
      margin = margin(r = 15)
    ),
    axis.line = element_line(),
    axis.text.y = element_text(size = rel(1.5), margin = margin(r =
                                                                  5)),
    axis.ticks.y = element_line(size = rel(1.1))
  )

#corrected ABCG2
ggplot(UMAPforPlot, aes(
  x = clusterForOrd,
  y = ABCG2,
  colour = as.factor(batch)
)) +
  geom_boxplot(width = 0.5) +
  geom_point(
    position = position_dodge(width = 0.5),
    alpha = 0.1,
    shape = 21,
    aes(fill = as.factor(batch))
  ) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('Cluster') +
  ylab('ABCG2 batch-corrected expression') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), margin = margin(t = 5)),
    axis.title.x = element_text(size = rel(1.5), margin = margin(t =
                                                                   15)),
    axis.title.y = element_text(
      size = rel(1.5),
      angle = 90,
      margin = margin(r = 15)
    ),
    axis.line = element_line(),
    axis.text.y = element_text(size = rel(1.5), margin = margin(r =
                                                                  5)),
    axis.ticks.y = element_line(size = rel(1.1))
  )

# ABCG2 by cell state:
ggplot(UMAPforPlot, aes(x = Cellstate, y = ABCG2, color = Cellstate)) +
  geom_violin(width = 0.5) +
  geom_point(position = position_dodge(width = 0.5),
             alpha = 0.5,
             shape = 21) +
  theme_void() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  xlab('Cell state') +
  ylab('ABCG2 batch-corrected expression') +
  guides(colour = F, fill = F) +
  theme(
    axis.text.x = element_text(size = rel(1.5), margin = margin(t = 5)),
    axis.title.x = element_text(size = rel(1.5), margin = margin(t =
                                                                   15)),
    axis.title.y = element_text(
      size = rel(1.5),
      angle = 90,
      margin = margin(r = 15)
    ),
    axis.line = element_line(),
    axis.text.y = element_text(size = rel(1.5), margin = margin(r =
                                                                  5)),
    axis.ticks.y = element_line(size = rel(1.1))
  )



# end cluster violins ####
#

#
