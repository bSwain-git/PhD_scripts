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
CFSE.files<-files.list[files.list %like% 'CFSE']
data_CFSE <- lapply(CFSE.files[c(-2)],
                       function(x) {
                         readRDS(paste0(directory0, x))
                       })

##annotate each with replicate number:
for (i in 1:length(data_CFSE)) {
  sampleNames(data_CFSE[[i]]) <-
    sampleNames(data_CFSE[[i]]) %>%
    lapply(function(name0) {
      paste0(name0, '.', i)
    })

  pData(data_CFSE[[i]])$name <-
    pData(data_CFSE[[i]])$name %>%
    lapply(function(name0) {
      paste0(name0, '.', i)
    })
  

  pData(data_CFSE[[i]])$replicate <- as.factor(i)
  class(pData(data_CFSE[[1]])$name) <- 'character'
}

for (j in 1:length(data_CFSE)){
  frame0<-data_CFSE[[j]]
  pData0<-pData(frame0)
  pData0$inhibitor[pData0$inhibitor == 'inhibited']<-'ko143'
  pData(data_CFSE[[j]])<-pData0
}


lapply(data_CFSE,pData)

# pData(data_CFSE[[2]])$cal<-F
# pData(data_CFSE[[2]])$pH<-NA
colOrder<-colnames(pData(data_CFSE[[1]]))
pData(data_CFSE[[2]])<-(pData(data_CFSE[[2]]))[colOrder]
# pData(data_CFSE[[3]])<-(pData(data_CFSE[[3]]))[colOrder]
# # 

data_CFSE_combined <- rbind2(data_CFSE[[1]], data_CFSE[[2]])
# data_CFSE_combined <- rbind2(data_CFSE_combined, data_CFSE[[3]])

pData(data_CFSE_combined)
FS <- data_CFSE_combined

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

 FSsubset<-FS[!(pData(FS)$inhibitor=='unstained') &
                !pData(FS)$cal == T &
                pData(FS)$replicate == 2]
 FSlivesubset<-FSlive[!(pData(FSlive)$inhibitor=='unstained') &
                !pData(FSlive)$cal == T &
                  pData(FSlive)$replicate == 2]
gg_allEvents <- ggcyto(FSsubset, aes(`FSC-A`, `SSC-A`)) +
  geom_hex(bins = 200,alpha=0.5) +
  facet_grid(cellType~inhibitor)+
  geom_hex(data=FSlivesubset,bins=200)+
  AllEventsPars +
  xlab('Forward-Scatter Area') +
  ylab('Side-Scatter Area')+
  geom_gate(rg)
# gg_allEvents
# ggsave('../Figs/CFSE/SSC-FSC plot CFSE.tiff',height = 6,width = 8,units = 'in')

##plotting to show single cells and doublets:
doubletsSingletsPars <-
  ggcyto_par_set(limits = list(x = c(300, 600), y = c(2e5, 8e5)))

# filter out doublets:
sg <- norm2Filter("FSC-W", "FSC-A", scale.factor = 2.5)


singletsOnly <- filter(FSlive, sg)
FSfiltered <- Subset(FSlive, singletsOnly)
FSfilteredsubset<-FSfiltered[!(pData(FSfiltered)$inhibitor=='unstained') &
               !pData(FSfiltered)$cal == T &
               pData(FSfiltered)$replicate == 2]


gg_singletsAndDoublets <-
  ggcyto(FSlivesubset, aes(x = `FSC-W`, y = `FSC-A`)) +
  geom_hex(bins = 200,alpha=0.5,colour='gray') +
  facet_grid(cellType~inhibitor)+
    geom_hex(data=FSfilteredsubset,bins=200)+
  doubletsSingletsPars +
  xlab('Forward-Scatter Width') +
  ylab('Forward-Scatter Area')
# gg_singletsAndDoublets
# ggsave('./Figs/CFSE/pH/FSCA-FSCW plot CFSE.tiff',height = 6,width = 8,units = 'in')

##plot unstained data:
ggplot(data = FSfiltered[(pData(FSfiltered)$replicate==2) & 
                           pData(FSfiltered)$dye=='unstained'], aes(x =CFSE_iso, fill = dye)) +
  geom_density(alpha = 0.5) +
  facet_grid(~ cellType  , scales = 'fixed') +
  # scale_x_log10() +
  xlim(-200,200)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('CFSE fluorescence (isosbestic emission)') +
  ylab('Population density ') +
  scale_fill_manual(values = c(
    # '#F8766D'
    # '#00BA38'
    '#619CFF'
  ))
ggplot(data = FSfiltered[(pData(FSfiltered)$replicate==2) & 
                           pData(FSfiltered)$dye=='unstained'], aes(x =CFSE_pH, fill = dye)) +
  geom_density(alpha = 0.5) +
  facet_grid(~ cellType  , scales = 'fixed') +
  # scale_x_log10() +
  xlim(-200,200)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('CFSE fluorescence (pH-sensitive emission)') +
  ylab('Population density ') +
  scale_fill_manual(values = c(
    # '#F8766D'
    # '#00BA38'
    '#619CFF'
  ))


library(ggpubr)

##Blank data by subtracting the LINEAR fluorescence of 'unstained' sample for iso comparisons
FS_BGsub <- FSfiltered

unstainedConditions<-sampleNames(FS_BGsub)[pData(FS_BGsub)$dye%like% 'unstained']

for (i in 1:length(unstainedConditions)){
  unstained0<-unstainedConditions[i]
  metaUnstained0<-pData(FS_BGsub[sampleNames(FS_BGsub)==unstained0])
  toBlank0<-sampleNames(FS_BGsub)[pData(FS_BGsub)$cellType== metaUnstained0$cellType&
                                    pData(FS_BGsub)$replicate== metaUnstained0$replicate]
  for (k in c('CFSE_pH','CFSE_iso')){
    
  unstainedVal0<-
    FS_BGsub@frames[[unstained0]]@exprs[, c(k)] %>%
    as.matrix() %>%
    median(na.rm = T)
  for (j in 1:length(toBlank0)){
    currentlyBlanking<-toBlank0[j]
    FS_BGsub@frames[[currentlyBlanking]]@exprs[, c(k)]<-
      FS_BGsub@frames[[currentlyBlanking]]@exprs[, c(k)]-unstainedVal0
  }
}
}

# saveRDS(FS_BGsub,'FS_final_data_CFSE.rds')
# FSfinal<-readRDS('FS_final_CFSE.rds')
FSfinal<-FS_BGsub

FSfinal <-
  flowCore::transform(FSfinal, 'log_CFSE_pH' = log10(`CFSE_pH`))
FSfinal <-
  flowCore::transform(FSfinal, 'log_CFSE_iso' = log10(`CFSE_iso`))
#isosbestic plot:
ggplot(data = FSfinal[pData(FSfinal)$cal == F &
                        !pData(FSfiltered)$dye == 'unstained' &
                        pData(FSfinal)$replicate==1],aes(x = log_CFSE_iso, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid(~cellType, scales = 'fixed') +
  xlim(3,5.5)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('log10(CFSE fluorescence), isosbestic emission (617 nm)') +
  ylab('Population density ') +
  scale_fill_manual(values = c('#F8766D',
                               '#00BA38',
                               # 'purple' ,
                               '#619CFF'))

#pH only plot:
ggplot(data = FSfinal[pData(FSfinal)$cal == F &
                        !pData(FSfiltered)$dye == 'unstained' &
                        pData(FSfinal)$replicate==1],aes(x =log_CFSE_pH, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid(~cellType, scales = 'fixed') +
  xlim(3,5.5)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('log10(CFSE fluorescence), pH-sensitive emission (525 nm)') +
  ylab('Population density ') +
  scale_fill_manual(values = c('#F8766D',
                               '#00BA38',
                               # 'purple' ,
                               '#619CFF'))


##ratio
FSfinal <-
  flowCore::transform(FSfinal, 'CFSE_ratio' = `CFSE_pH`/`CFSE_iso`)

#ratio plot:
ggplot(data = FSfinal[pData(FSfinal)$cal == F &
                        !pData(FSfiltered)$dye == 'unstained' &
                        pData(FSfinal)$replicate==1],aes(x = CFSE_ratio, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid(~cellType, scales = 'fixed') +
  xlim(1,3)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('CFSE ratio') +
  ylab('Population density ') +
  scale_fill_manual(values = c('#F8766D',
                               '#00BA38',
                               # 'purple' ,
                               '#619CFF'))


ggplot(data = FSfinal[pData(FSfinal)$cal == T],aes(x = CFSE_ratio, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid(pH~cellType, scales = 'fixed') +
  xlim(0,5)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('CFSE ratio') +
  ylab('Population density ') +
  scale_fill_manual(values = c('#F8766D',
                               '#00BA38',
                               # 'purple' ,
                               '#619CFF'))



calCons<-sampleNames(FSfinal)[pData(FSfinal)$cal==T]
calVals<-matrix(ncol = 2,nrow=length(calCons))
colnames(calVals)<-c('CFSE_pH','CFSE_iso')
rownames(calVals)<-(calCons)
for (i in 1:length(calCons)){
  calCon0<-calCons[i]
  medians0<-c(median(FSfinal[[calCon0]]@exprs[,'CFSE_pH'],na.rm = T),
              median(FSfinal[[calCon0]]@exprs[,'CFSE_iso'],na.rm = T)
  )
  calVals[i,]<-medians0
}
pHcals<-as.data.frame(calVals)%>%
  rownames_to_column(var='sample')
pHcals$pH<-numextract(pHcals$sample)
pHcals$pH[pHcals$pH == 7.4]<-7.33
pHcals$cellType<-str_sub(pHcals$sample,end=2)
pHcals$cellType[pHcals$cellType == 'EQ']<-'E211Q'
pHcals$inhibitor[pHcals$sample %like% '-k'] <- 'active'
pHcals$inhibitor[pHcals$sample %like% '\\+k'] <- 'ko143'
pHcals$ratio<-pHcals$CFSE_pH/pHcals$CFSE_iso

ggplot(pHcals,aes(x=pH,y=CFSE_pH/CFSE_iso))+
  geom_point()+
  facet_grid(inhibitor~cellType)+
  geom_smooth(method='lm')+
  ylab('CFSE ratio')

FSinterp<-FSfinal
#create a dummy column to hold interpolated data:
FSinterp <-
  flowCore::transform(FSinterp, 'CFSE_pH_interp' = 0*`CFSE_pH`)

## using a linear model combinig all cals:
pHmodel<-pHcals%>%
  lm(pH ~ ratio,.)
library(Hmisc)

for (i in 1:length(sampleNames(FSinterp))){
  sample0<-sampleNames(FSinterp)[i]

  ratio0<-FSinterp[[sample0]]@exprs[, c('CFSE_ratio')]%>%
    as.numeric()
  
  meta0<-pData(FSinterp)[pData(FSinterp)$name == sample0,]
  
  if (meta0$cal ==F & !meta0$inhibitor == 'unstained'){
    cellType0<-meta0$cellType%>%as.character()
    inhibitor0<-(meta0$inhibitor)%>%as.character()
    FF0<-exprs(FSinterp[[sample0]])%>%as.matrix()
    interpRatio0<-
      pHmodel%>%
      # subset(pHcals,(cellType == cellType0 & inhibitor == (inhibitor0)))%>%
      # lm(pH ~ ratio,.)%>%
      predict(.,newdata=data.frame(ratio=ratio0))%>%
      data.frame('CFSE_pH'=.)%>%
      asNumericMatrix()
    FF0[,'CFSE_pH_interp']<-interpRatio0
    exprs(FSinterp[[sample0]])<-FF0
  }
}


#interpolated pH plot:
ggplot(data = FSinterp[pData(FSinterp)$cal == F&
                        !pData(FSinterp)$dye == 'unstained'],aes(x = CFSE_pH_interp, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid(~cellType, scales = 'fixed') +
  xlim(4,8)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('Interpolated pH') +
  ylab('Population density ') +
  scale_fill_manual(values = c('#F8766D',
                               '#00BA38',
                               # 'purple' ,
                               '#619CFF'))

interppHList<-list()
for (f0 in 1:length(sampleNames(FSinterp))){
  frameName0<-sampleNames(FSinterp)[f0]
  interppHList[[f0]]<-FSinterp@frames[[frameName0]]@exprs[,'CFSE_pH_interp']%>%median()
  
}
interppH<-data.frame('sampleName'=sampleNames(FSinterp),'pH'=unlist(interppHList))
interppH<-interppH[!interppH$pH == 0,]

interppH$replicate<-str_sub(interppH$sampleName,start=-1)
interppH$cellType<-str_sub(interppH$sample,end=2)
interppH$cellType[interppH$cellType == 'EQ']<-'E211Q'
interppH$inhibitor[interppH$sample %like% '-k'] <- 'active'
interppH$inhibitor[interppH$sample %like% '\\+k'] <- 'ko143'

group_by(interppH,replicate,cellType,inhibitor)%>%
  summarise('meanpH' = mean(pH))%>%
ggplot(.,aes(x=inhibitor,y=meanpH,colour=cellType))+
  geom_point()+
  geom_boxplot()+
  stat_compare_means()+
  # facet_grid(~replicate)+
  scale_colour_manual(values = c('#00BA38','#F8766D'
                               
                               # 'purple' ,
                               # '#619CFF'
                               ))


#create a dummy column to hold interpolated data:
FSinterp <-
  flowCore::transform(FSinterp, 'CFSE_pH_interp_individual' = 0*`CFSE_pH`)

## using a linear model for each sample type:


for (i in 1:length(sampleNames(FSinterp))){
  sample0<-sampleNames(FSinterp)[i]
  
  ratio0<-FSinterp[[sample0]]@exprs[, c('CFSE_ratio')]%>%
    as.numeric()
  
  meta0<-pData(FSinterp)[pData(FSinterp)$name == sample0,]
  
  if (meta0$cal ==F & !meta0$inhibitor == 'unstained'){
    
    calValsRelevant<-subset(pHcals,(cellType == meta0$cellType & inhibitor == meta0$inhibitor))
    pHmodelIndividual<-calValsRelevant%>%
      lm(pH ~ ratio,.)
    
    cellType0<-meta0$cellType%>%as.character()
    inhibitor0<-(meta0$inhibitor)%>%as.character()
    FF0<-exprs(FSinterp[[sample0]])%>%as.matrix()
    interpRatio0<-
      pHmodelIndividual%>%
      predict(.,newdata=data.frame(ratio=ratio0))%>%
      data.frame('CFSE_pH'=.)%>%
      asNumericMatrix()
    FF0[,'CFSE_pH_interp_individual']<-interpRatio0
    exprs(FSinterp[[sample0]])<-FF0
  }
}


#interpolated pH plot:
ggplot(data = FSinterp[pData(FSinterp)$cal == F&
                         !pData(FSinterp)$dye == 'unstained'],aes(x = CFSE_pH_interp_individual, fill = inhibitor)) +
  geom_density(alpha = 0.5) +
  facet_grid(~cellType, scales = 'fixed') +
  xlim(4,9)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab('Interpolated pH') +
  ylab('Population density ') +
  scale_fill_manual(values = c('#F8766D',
                               '#00BA38',
                               # 'purple' ,
                               '#619CFF'))

interppHListIndividual<-list()
for (f0 in 1:length(sampleNames(FSinterp))){
  frameName0<-sampleNames(FSinterp)[f0]
  interppHListIndividual[[f0]]<-FSinterp@frames[[frameName0]]@exprs[,'CFSE_pH_interp_individual']%>%median()
  
}
interppHIndividual<-data.frame('sampleName'=sampleNames(FSinterp),'pH.int'=unlist(interppHListIndividual))
interppHIndividual<-interppHIndividual[!interppHIndividual$pH == 0,]

interppHfinal<-inner_join(interppH,interppHIndividual)

interppHfinal$cellType
interppHfinal%>%mutate(cellType=fct_relevel(cellType,
                              'WT','E211Q'))%>%
group_by(replicate,cellType,inhibitor)%>%

    summarise('meanpH' = mean(pH))%>%
  ggplot(.,aes(x=cellType,y=meanpH,fill=inhibitor,color=inhibitor,alpha=cellType))+
  # geom_point()+
  geom_boxplot(size=1)+
  ylab((expression(Calculated~pH[IN])))+
  xlab('ABCG2 variant')+
  # stat_compare_means()+
  # facet_grid(~replicate)+
  scale_fill_manual(values = c('#F8766D','#00BA38'),guide=F)+
  scale_color_manual(values = c('#F8766D','#00BA38'))+
  scale_alpha_discrete(range=c(0,0.5),guide=F)+
  ylim(6.5,8)+
  theme(panel.background = element_blank(),
        axis.title = element_text(face='bold'),
        axis.text  = element_text(face='bold'),
        panel.grid.major.x = element_blank(),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                                 colour = "lightgray"),
    panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                    colour = "gray"))
saveRDS(interppHfinal,'interpolated_pH_flow.rds')
