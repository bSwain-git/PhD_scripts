source(
  'C:/Users/Brendan/OneDrive - University Of Cambridge/Data/R analysis/generalFunctions.R'
)
setwd.here()
library(flowCore)
library(tidyverse)
library(data.table)
set.seed(1811)

#read in files :
directory0 <- './Split raw data/'
files.list <- list.files(directory0, pattern = '.rds')
files.list
data_EtBr <- lapply(files.list[files.list %like% 'EtBr'],
                       function(x) {
                         readRDS(paste0(directory0, x))
                       })

#One dataset had no pH annotation (was done at normal buffer pH == 7.4)
pData(data_EtBr[[3]])$pH <- 7.4

##annotate each with replicate number:
for (i in 1:length(data_EtBr)) {
  sampleNames(data_EtBr[[i]]) <-
    sampleNames(data_EtBr[[i]]) %>%
    lapply(function(name0) {
      paste(name0, 'day', i)
    })

  pData(data_EtBr[[i]])$name <-
    pData(data_EtBr[[i]])$name %>%
    lapply(function(name0) {
      paste(name0, 'day', i)
    })
  

  pData(data_EtBr[[i]])$day <- as.factor(i)
  class(pData(data_EtBr[[1]])$name) <- 'character'
}

#combine the sets together
EtBr_combined <- rbind2(data_EtBr[[1]], data_EtBr[[2]])
EtBr_combined <- rbind2(EtBr_combined, data_EtBr[[3]])
EtBr_combined <- rbind2(EtBr_combined, data_EtBr[[4]])
EtBr_combined <- rbind2(EtBr_combined, data_EtBr[[5]])
pData(EtBr_combined)
FS <- EtBr_combined

#set the correct class for the sample metadata columns:
pData(FS)$name <- as.character(pData(FS)$name)
pData(FS)[!colnames(pData(FS)) == 'name'] <-
  lapply(pData(FS)[!colnames(pData(FS)) == 'name'], as.factor)


library(ggcyto)
library(flowViz)
library(flowStats)
library(scales)

#plotting to remove debris from analysis
AllEventsPars <-
  ggcyto_par_set(limits = list(x = c(0, 1000000), y = c(0, 1000000)))



##remove debris from analysis
rg <-
  rectangleGate(
    "FSC-A" = c(270000, 900000),
    'SSC-A' = c(120000, 600000),
    filterId = 'rectangle'
  )
ng <- norm2Filter("SSC-A",
                  "FSC-A",
                  scale.factor = 1.5,
                  filterId = 'norm')
lc <- ng %&% rg
liveCells <- filter(FS, lc)
summary(liveCells)
FSlive <- Subset(FS, liveCells)

FSsubset<-FS[(pData(FS)$day==2) & 
               !(pData(FS)$inhibitor=='unstained')]
FSliveSubset<-FSlive[pData(FSlive)$day==2 & 
                       !(pData(FSlive)$inhibitor=='unstained')]
gg_allEvents <- ggcyto(FSsubset, aes(`FSC-A`, `SSC-A`)) +
  geom_hex(bins = 200,alpha=0.5) +
  facet_grid(cellType*pH~inhibitor)+
  geom_hex(data=FSliveSubset,bins=200)+
  AllEventsPars +
  xlab('Forward-Scatter Area') +
  ylab('Side-Scatter Area')+
  geom_gate(rg)
# gg_allEvents
# ggsave('../Figs/EtBr/SSC-FSC plot EtBr.tiff',height = 6,width = 8,units = 'in')

##plotting to show single cells and doublets:
doubletsSingletsPars <-
  ggcyto_par_set(limits = list(x = c(300, 600), y = c(2e5, 8e5)))



# filter out doublets:
sg <- norm2Filter("FSC-W", "FSC-A", scale.factor = 2.5)


singletsOnly <- filter(FSlive, sg)
FSfiltered <- Subset(FSlive, singletsOnly)
FSfilteredSubset<-FSfiltered[(pData(FSfiltered)$day==2) & 
                       !(pData(FSfiltered)$inhibitor=='unstained')]


gg_singletsAndDoublets <-
  ggcyto(FSliveSubset, aes(x = `FSC-W`, y = `FSC-A`)) +
  geom_hex(bins = 200,alpha=0.5,colour='gray') +
  facet_grid(cellType*pH~inhibitor)+
    geom_hex(data=FSfilteredSubset,bins=200)+
  doubletsSingletsPars +
  xlab('Forward-Scatter Width') +
  ylab('Forward-Scatter Area')
# gg_singletsAndDoublets
# ggsave('../Figs/EtBr/FSCA-FSCW plot EtBr.tiff',height = 6,width = 8,units = 'in')


# log-transform EtBr fluorescence
FSfiltered <-
  flowCore::transform(FSfiltered, 'log_EtBr' = log10(`EtBr`))
library(ggpubr)

##Blank data by subtracting the fluorescence of 'unstained' sample
FS_BGsub <- FSfiltered

unstainedConditions<-sampleNames(FS_BGsub)[pData(FS_BGsub)$dye%like% 'unstained']

for (i in 1:length(unstainedConditions)){
  unstained0<-unstainedConditions[i]
  metaUnstained0<-pData(FS_BGsub[sampleNames(FS_BGsub)==unstained0])
  toBlank0<-sampleNames(FS_BGsub)[pData(FS_BGsub)$cellType== metaUnstained0$cellType&
                                    pData(FS_BGsub)$pH== metaUnstained0$pH&
                                    pData(FS_BGsub)$day== metaUnstained0$day]
  unstainedVal0<-
    FS_BGsub@frames[[unstained0]]@exprs[, c('EtBr')] %>%
    as.matrix() %>%
    median(na.rm = T)
  for (j in 1:length(toBlank0)){
    currentlyBlanking<-toBlank0[j]
    FS_BGsub@frames[[currentlyBlanking]]@exprs[, c('EtBr')]<-
      FS_BGsub@frames[[currentlyBlanking]]@exprs[, c('EtBr')]-unstainedVal0
  }
  
}

# saveRDS(FS_BGsub,'FS_final_EtBr.rds')
# FSfinal<-readRDS('FS_final_EtBr.rds')
FSfinal<-FS_BGsub

##plot unstained data:
ggplot(data = FSfiltered[(pData(FSfiltered)$day==2) & 
                           pData(FSfiltered)$dye=='unstained'], aes(x =EtBr, fill = dye)) +
  geom_density(alpha = 0.5) +
  facet_grid(~ cellType  , scales = 'fixed') +
  # scale_x_log10() +
  xlim(-2000,3000)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('EtBr fluorescence') +
  ylab('Population density ') +
  scale_fill_manual(values = c(
    # '#F8766D'
    # '#00BA38'
    '#619CFF'
  ))

##plot to show pH effect using data from day 1: 

ggplot(data = FSfinal[pData(FSfinal)$inhibitor%in% c('active','ko143')&
                        pData(FSfinal)$pH %in% c(6,7.4) &
                        pData(FSfinal)$day %in% c(1)], aes(x = log_EtBr, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid( pH~ cellType, scales = 'fixed') +
  # scale_x_log10(limits=c(1,1e6)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('EtBr fluorescence') +
  ylab('Population density ') +
  scale_fill_manual(values = c('#F8766D',
                               '#00BA38',
                               # 'purple' ,
                               '#619CFF'))

##plot to show nigericin effect using data from day 4 and 5:
ggplot(data = FSfinal[(pData(FSfinal)$inhibitor%in% c('active','nigericin')&
                        pData(FSfinal)$pH %in%  c(6)&
                        pData(FSfinal)$day %in% c(5))|
                        pData(FSfinal)$inhibitor%in% c('active','nigericin')&
                        pData(FSfinal)$pH %in%  c(7.4)&
                        pData(FSfinal)$day %in% c(4)], aes(x = log_EtBr, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid(pH~cellType, scales = 'fixed') +
  # scale_x_log10(limits=c(1,1e6)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('log10(EtBr fluorescence)') +
  ylab('Population density ') +
  scale_fill_manual(values = c('#F8766D', 
                                'purple','#619CFF'
                               ))
