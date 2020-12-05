sce <- readRDS('sce_final.rds')
set.seed(1811)
library(tidyverse)
library(scater)

corrected <- altExp(sce, 'MNN_corrected')
colData(corrected) <- colData(sce)
sce.naive <- sce[, sce$Characteristics.phenotype. == 'naive']
corrected.naive <-
  corrected[, corrected$Characteristics.phenotype. == 'naive']

corrected.naive <-
  runDiffusionMap(corrected.naive, dimred = 'corrected')

accessions <- list()
accessions[['ABCG2']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'ABCG2']
accessions[['SLC9A1']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'SLC9A1']
accessions[['SLC4A4']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'SLC4A4']
accessions[['SLC9A3']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'SLC9A3']
accessions[['SLC4A10']] <-
  rowData(sce)$accession[rowData(sce)$symbol %in% 'SLC4A10']
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

# plotDiffusionMap(corrected.naive,colour_by=accessions$ABCG2,by_exprs_values='reconstructed')

forDMplot <-
  data.frame(
    'DM1' = reducedDim(corrected.naive, 'DiffusionMap')[, 1],
    'DM2' = reducedDim(corrected.naive, 'DiffusionMap')[, 2],
    'ABCG2.recon' = assay(corrected.naive, 'reconstructed')[rownames(corrected.naive) == accessions$ABCG2, ],
    'ABCG2.lc' = assay(sce.naive, 'logcounts')[rownames(sce.naive) == accessions$ABCG2, ],
    'GATA3.lc' = assay(sce.naive, 'logcounts')[rownames(sce.naive) == accessions$GATA3, ],
    'GATA3.recon' = assay(corrected.naive, 'reconstructed')[rownames(corrected.naive) == accessions$GATA3, ],
    'SLC9A1.recon' = assay(corrected.naive, 'reconstructed')[rownames(corrected.naive) == accessions$SLC9A1, ],
    'cluster' = corrected.naive$clustWalktrap,
    'batch' = as.factor(corrected.naive$batch)
  )

#colour DM by cluster
ggplot(forDMplot, aes(x = DM1, y = DM2, fill = cluster)) +
  geom_point(
    alpha = 1,
    pch = 21,
    colour = 'gray',
    size = 2
  ) +
  theme_void() +
  guides() +
  labs(fill = 'Cluster') +
  theme(
    # axis.text.x = element_text(size=rel(1),margin=margin(t=5)),
    #     axis.text.y = element_text(size=rel(1),margin=margin(r=5)),
    # axis.title.x = element_markdown(size=rel(1.5),margin = margin(t=10,b=5)),
    # axis.title.y = element_markdown(size=rel(1.5),angle = 90,margin = margin(r=10)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.5, 0.3)
    # axis.ticks = element_line(size=rel(1.1)),
    # axis.line = element_line()
  ) +
  scale_fill_manual(values = list('1' = '#F8766D',
                                  '2' = '#C59900',
                                  '5' = '#E76BF3'))


#colour DM by ABCG2
ggplot(forDMplot, aes(x = DM1, y = DM2, fill = ABCG2.recon)) +
  geom_point(
    alpha = 1,
    pch = 21,
    colour = 'black',
    size = 2
  ) +
  theme_void() +
  guides() +
  labs(fill = expression(atop(italic('ABCG2'), 'expression'))) +
  guides(fill = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  scale_fill_viridis_c(option = 'B') +
  theme(
    # axis.text.x = element_text(size=rel(1),margin=margin(t=5)),
    #     axis.text.y = element_text(size=rel(1),margin=margin(r=5)),
    # axis.title.x = element_markdown(size=rel(1.5),margin = margin(t=10,b=5)),
    # axis.title.y = element_markdown(size=rel(1.5),angle = 90,margin = margin(r=10)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.5, 0.3)
    # axis.ticks = element_line(size=rel(1.1)),
    # axis.line = element_line()
  )


#colour DM by GATA3
ggplot(forDMplot, aes(x = DM1, y = -DM2, fill = GATA3.recon)) +
  geom_point(
    alpha = 1,
    pch = 21,
    colour = 'black',
    size = 2
  ) +
  theme_void() +
  guides() +
  labs(fill = expression(atop(italic('GATA3'), 'expression'))) +
  guides(fill = guide_colourbar(
    barwidth = 1,
    barheight = 5,
    title.position = 'left',
    size = rel(1.5)
  )) +
  scale_fill_viridis_c(option = 'B') +
  theme(
    # axis.text.x = element_text(size=rel(1),margin=margin(t=5)),
    #     axis.text.y = element_text(size=rel(1),margin=margin(r=5)),
    # axis.title.x = element_markdown(size=rel(1.5),margin = margin(t=10,b=5)),
    # axis.title.y = element_markdown(size=rel(1.5),angle = 90,margin = margin(r=10)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.5, 0.3)
    # axis.ticks = element_line(size=rel(1.1)),
    # axis.line = element_line()
  )
library(ggtext)

# ABCG2 vs GATA3:
ggplot(forDMplot, aes(x = GATA3.recon, y = ABCG2.recon)) +
  geom_point(alpha = 1) +
  geom_smooth(
    se = F,
    colour = 'red',
    method = 'loess',
    span = .9
  ) +
  theme_void() +
  # guides(fill=guide_colourbar(barwidth =1, barheight = 5,title.position = 'left',size=rel(1.5)))+
  theme(
    axis.text.x = element_text(size = rel(1), margin = margin(t = 5)),
    axis.text.y = element_text(size = rel(1), margin = margin(r = 5)),
    axis.title.x = element_markdown(size = rel(1.5), margin = margin(t =
                                                                       10, b = 5)),
    axis.title.y = element_markdown(
      size = rel(1.5),
      angle = 90,
      margin = margin(r = 10)
    ),
    axis.ticks = element_line(size = rel(1.1)),
    axis.line = element_line(),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(.3, 0.8)
  ) +
  labs(x = '*GATA3*',
       y = '*ABCG2*',
       fill = '- (Diffusion\ncomponent 1)') +
  scale_colour_manual(values = list('1' = '#F8766D',
                                    '2' = '#C59900',
                                    '5' = '#E76BF3'))
