sce <- readRDS('sce_final.rds')
library(tidyverse)
library(scater)

corrected <- altExp(sce, 'MNN_corrected')
uncorrected <- sce
rm(sce)

UMAPforPlot <- data.frame(
  'batch' = uncorrected$batch,
  'cUMAP1' = reducedDim(corrected, 'UMAP')[, 1],
  'cUMAP2' = reducedDim(corrected, 'UMAP')[, 2],
  'ucUMAP1' = reducedDim(uncorrected, 'UMAP')[, 1],
  'ucUMAP2' = reducedDim(uncorrected, 'UMAP')[, 2]
  
)

#uncorrected:
# This is figure 5.5A:
ggplot(UMAPforPlot, aes(
  x = ucUMAP1,
  y = ucUMAP2,
  colour = factor(batch)
)) +
  geom_point(aes(alpha = 0.5)) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = 'Batch') +
  guides(alpha = F) +
  # guides(colour=F)+
  # scale_colour_viridis_c(option='B')+
  theme(
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(0.6, 0.7)
  )


#corrected:
# This is figure 5.5B:
ggplot(UMAPforPlot, aes(
  x = cUMAP1,
  y = cUMAP2,
  colour = factor(batch)
)) +
  geom_point(aes(alpha = 0.5)) +
  theme_void() +
  # annotate(geom='label',x=2,y=1,label='Naive',colour='#F8766D',size=8)+
  # annotate(geom='label',x=15,y=-1,label='Primed',colour='#00BFC4',size=8)+
  labs(colour = 'Batch') +
  guides(alpha = F) +
  # guides(colour=F)+
  # scale_colour_viridis_c(option='B')+
  theme(
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5)),
    legend.position = c(0.6, 0.7)
  )
