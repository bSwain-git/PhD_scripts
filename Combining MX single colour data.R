source(
  'C:/Users/Brendan/OneDrive - University Of Cambridge/Data/R analysis/generalFunctions.R'
)
setwd.here()
library(flowCore)
library(tidyverse)
library(data.table)

#read in files :
directory0 <- './Split raw data/'
files.list <- list.files(directory0, pattern = '.rds')
files.list
data_MX <- lapply(files.list[files.list %like% 'MX'],
                       function(x) {
                         readRDS(paste0(directory0, x))
                       })


##annotate each with replicate number:
for (i in 1:length(data_MX)) {
  sampleNames(data_MX[[i]]) <-
    sampleNames(data_MX[[i]]) %>%
    lapply(function(name0) {
      paste(name0, 'rep', i)
    })
  #We don't need pH annotations for these data:
  
  pData(data_MX[[i]])$pH<-NULL
  
  
  pData(data_MX[[i]])$name <-
    pData(data_MX[[i]])$name %>%
    lapply(function(name0) {
      paste(name0, 'rep', i)
    })
  
  pData(data_MX[[i]])$replicate <- as.factor(i)
  class(pData(data_MX[[1]])$name) <- 'character'
}

#combine the sets together
MX_combined <- rbind2(data_MX[[1]], data_MX[[2]])
MX_combined <- rbind2(MX_combined, data_MX[[3]])

pData(MX_combined)
# sampleNames(MX_combined)<- stringr::str_replace(sampleNames(MX_combined),'un -k','unstained')
FS <- MX_combined


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

FSsubset<-FS[
  (pData(FS)$replicate==2) &
    !(pData(FS)$inhibitor=='unstained')]
FSliveSubset<-FSlive[
  pData(FSlive)$replicate==2 &
    !(pData(FSlive)$inhibitor=='unstained')]

gg_allEvents <- ggcyto(FSsubset, aes(`FSC-A`, `SSC-A`)) +
  geom_hex(bins = 200,alpha=0.5) +
  geom_hex(data=FSliveSubset,bins=200)+
  facet_grid(cellType~inhibitor)+
    AllEventsPars +
  xlab('Forward-Scatter Area') +
  ylab('Side-Scatter Area')+
  geom_gate(rg)

# gg_allEvents

##plotting to show single cells and doublets:
doubletsSingletsPars <-
  ggcyto_par_set(limits = list(x = c(300, 600), y = c(2e5, 8e5)))



# filter out doublets:
sg <- norm2Filter("FSC-W", "FSC-A", scale.factor = 2.5)


singletsOnly <- filter(FSlive, sg)
FSfiltered <- Subset(FSlive, singletsOnly)
FSfilteredSubset<-FSfiltered[(pData(FSfiltered)$replicate==2) & 
                               !(pData(FSfiltered)$inhibitor=='unstained')]



gg_singletsAndDoublets <-
  ggcyto(FSliveSubset, aes(x = `FSC-W`, y = `FSC-A`)) +
  geom_hex(bins = 200,alpha=0.5,colour='gray') +
  geom_hex(data=FSfilteredSubset,bins=200)+
  facet_grid(cellType~inhibitor)+
    doubletsSingletsPars +
  xlab('Forward-Scatter Width') +
  ylab('Forward-Scatter Area')
# gg_singletsAndDoublets

# log-transform MX fluorescence,
FSfiltered <-
  flowCore::transform(FSfiltered, 'log_MX' = log(`MX`))
library(ggpubr)


#plot unstained data:
ggplot(data = FSfiltered[(pData(FSfiltered)$replicate==2) & 
                           pData(FSfiltered)$dye=='unstained'], aes(x =MX, fill = dye)) +
  geom_density(alpha = 0.5) +
  facet_grid(~ cellType  , scales = 'fixed') +
  # scale_x_log10() +
  xlim(-500,1000)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('MX fluorescence') +
  ylab('Population density ') +
  scale_fill_manual(values = c(
    # '#F8766D'
    # '#00BA38'
    '#619CFF'
  ))

##Blank data by subtracting the fluorescence of 'unstained' sample
FS_BGsub <- FSfiltered

unstainedConditions<-sampleNames(FS_BGsub)[pData(FS_BGsub)$dye=='unstained']

for (i in 1:length(unstainedConditions)){
  unstained0<-unstainedConditions[i]
  metaUnstained0<-pData(FS_BGsub[sampleNames(FS_BGsub)==unstained0])
  toBlank0<-sampleNames(FS_BGsub)[pData(FS_BGsub)$cellType== metaUnstained0$cellType&
                       pData(FS_BGsub)$replicate== metaUnstained0$replicate]
  unstainedVal0<-
    FS_BGsub@frames[[unstained0]]@exprs[, c('MX')] %>%
    as.matrix() %>%
    median(na.rm = T)
  for (j in 1:length(toBlank0)){
    currentlyBlanking<-toBlank0[j]
    FS_BGsub@frames[[currentlyBlanking]]@exprs[, c('MX')]<-
      FS_BGsub@frames[[currentlyBlanking]]@exprs[, c('MX')]-unstainedVal0
  }
  
  }

#plot blanked data:
ggplot(data = FS_BGsub[!pData(FS_BGsub)$dye=='unstained'&
                         pData(FS_BGsub)$replicate == 2], aes(x = log_MX, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid(. ~ cellType  , scales = 'fixed') +
  # scale_x_log10() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('log10(MX fluorescence)') +
  ylab('Population density ')+
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF'))

saveRDS(FS_BGsub,'FS_final_MX.rds')


