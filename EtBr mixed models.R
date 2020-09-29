source(
  'C:/Users/Brendan/OneDrive - University Of Cambridge/Data/R analysis/generalFunctions.R'
)
setwd.here()
library(flowCore)
library(tidyverse)
library(data.table)

FSfinal<-readRDS('FS_final_EtBr.rds')

#the only data sets that exhibit a mixed population are the stained Wt samples without ko143:
samplesForMM<-sampleNames(FSfinal)[pData(FSfinal)$cellType =='WT' &
                                     pData(FSfinal)$inhibitor == 'active']
##include the inhibited equivalent samples for comparison, using simple normal distribution: 
samplesForNorm<-sampleNames(FSfinal)[pData(FSfinal)$cellType =='WT' &
                                       pData(FSfinal)$inhibitor%in%c('nigericin','ko143')]

library(ggplot2)
library(MASS)
library(mixtools)

set.seed(1811)
MMDistList<-list()
MMDataList<-list()
for (i in 1:length(samplesForMM)){
  
  sample0<-samplesForMM[i]
  MMdata0 <- FSfinal@frames[[sample0]]@exprs[, c('log_EtBr')] %>%
    as.matrix()%>%
    .[!is.na(.)]
  # pData0<-pData(FSfinal[pData(FSfinal)$name == sample0])
  # 
  # pDataInhibNames<-pData(FSfinal)$name[(pData(FSfinal)$cellType == 'WT'&
  #                  pData(FSfinal)$inhibitor == 'ko143'&
  #                  pData(FSfinal)$dye == 'EtBr'&
  #                  pData(FSfinal)$pH == pData0$pH)] 
  # pDataInhib<-FSfinal[pData(FSfinal)$name %in% pDataInhibNames]
  # medians0<-fsApply(pDataInhib, function(x){
  #   x@exprs[,c('EtBr')]%>%mean()
  # })%>%median()
    
    MM0 <- normalmixEM(x = MMdata0,
                     fast=FALSE,
                     maxit=10000,
                     epsilon = 1e-16,
                     # mean.constr = c(NA,medians0),
                     maxrestarts=1000,
                     k = 2,
                     verb=F)
  MMDistList[[i]]<-MM0
  MMDataList[[i]]<-MMdata0
  names(MMDistList)[i]<-sample0
  names(MMDataList)[i]<-sample0
}

highLam<-list()
lowLam<-list()
nameList<-list()
for (it in 1:length(MMDistList)){
  whichLam<-which.max(MMDistList[[it]]$mu)
  
  
  highLam[[it]]<-MMDistList[[it]]$lambda[whichLam]
  lowLam[[it]]<-MMDistList[[it]]$lambda[!c(1,2) ==whichLam]
  
  nameList[[it]]<-names(MMDistList)[it]
}

EtMMlambdas<-data.frame('high'=unlist(highLam),'low'=unlist(lowLam),name=unlist(nameList))
write.csv(EtMMlambdas,('Lambdas/Et.csv'))

NormDistList<-list()
NormDataList<-list()
for (i in 1:length(samplesForNorm)){
  
  sample0<-samplesForNorm[i]
  Normdata0 <- FSfinal@frames[[sample0]]@exprs[, c('log_EtBr')] %>%
    as.matrix()%>%
    .[!is.na(.)]
  Norm0 <- fitdistr(Normdata0,'normal')
  NormDistList[[i]]<-Norm0
  NormDataList[[i]]<-Normdata0
  names(NormDistList)[i]<-sample0
  names(NormDataList)[i]<-sample0
}



plot_mix_comps <- function(x, mu, sigma, lam) {
  (lam * dnorm((x), mu, sigma))
}

####Plot ko143 example plots####
##example ko143 inhibition plot samples come from 'day 2':
normExamples<-(NormDataList)[names(NormDataList) %like% 'day 2']
normDists<-(NormDistList)[names(NormDistList) %like% 'day 2']
MMexamples<-(MMDataList)[names(MMDataList) %like% 'day 2']
MMDists<-(MMDistList)[names(MMDistList) %like% 'day 2']

Ko143PlotList<-list()

for(j in 1:length(MMexamples)){
  MMexample0<-MMexamples[j]
  MMDist0<-MMDists[[j]]
  if (names(MMexample0) %like% '6'){
    pH0<-'6'
    matchedNormExamples<-normExamples[names(normExamples) %like% '6']
    matchedNormDists<-normDists[names(normDists) %like% '6']
  }else{
    pH0<-'7.4'
    matchedNormExamples<-normExamples[!names(normExamples) %like% '6']
    matchedNormDists<-normDists[!names(normDists) %like% '6']
  }
  
  ko143Example<-matchedNormExamples[names(matchedNormExamples) %like% '\\+k']
  ko143Dist<-matchedNormDists[names(matchedNormDists) %like% '\\+k']
  
  Ko143PlotList[[pH0]]<-
    
    ##histograms of raw data:####
  ggplot() +
    geom_histogram(
      data = data.frame('x'=unlist(ko143Example)),
      aes(x,..density.., 
          fill = '+ 1 uM ko143'),
      binwidth = 0.005,
      alpha = 0.5
    ) +
    geom_histogram(
      data = data.frame('x'=unlist(MMexample0)),
      aes(x,..density.., 
          fill = '+ DMSO'),
      binwidth = 0.005,
      alpha = 0.5
    )+
    ##end####
  ##lines showing fitted normal distribution for inhibited sample:####
  stat_function(
    geom = "line",
    fun = plot_mix_comps,
    args = list(
      mu = ko143Dist[[1]]$estimate['mean'],
      sigma = ko143Dist[[1]]$estimate['sd'],
      lam = 1
    ),
    lwd = 1,
    aes(colour = 'inhibited'),
    show.legend = F
  ) +
    ##end####
  ##lines showing mixed models for active samples:####
  stat_function(
    geom = "line",
    fun = plot_mix_comps,
    args = list(MMDist0$mu[which.max(c(MMDist0$mu[1],MMDist0$mu[2]))], 
                MMDist0$sigma[which.max(c(MMDist0$mu[1],MMDist0$mu[2]))], 
                lam = MMDist0$lambda[which.max(c(MMDist0$mu[1],MMDist0$mu[2]))]),
    lwd = 1,
    aes(colour = 'untransfected')
  ) +
    stat_function(
      geom = "line",
      fun = plot_mix_comps,
      args = list(MMDist0$mu[which.min(c(MMDist0$mu[1],MMDist0$mu[2]))],
                  MMDist0$sigma[which.min(c(MMDist0$mu[1],MMDist0$mu[2]))],
                  lam = MMDist0$lambda[which.min(c(MMDist0$mu[1],MMDist0$mu[2]))]),
      lwd = 1,
      aes(colour = 'transfected')
    )+
    ##end####
  ##Vertical lines indicating various means:####
  geom_vline(xintercept = ko143Dist[[1]]$estimate['mean'], color = 'forestgreen',size=1) +
    geom_vline(xintercept = MMDist0$mu[which.max(c(MMDist0$mu[1],MMDist0$mu[2]))], color = 'violetred4',size=1) +
    geom_vline(xintercept = MMDist0$mu[which.min(c(MMDist0$mu[1],MMDist0$mu[2]))], color = 'red',size=1) +
    geom_vline(xintercept = mean(as.numeric(unlist(MMexample0))), color = 'black',size=1)+
    ##end####
  ##Label the means:####
  annotate('label',
           x = MMDist0$mu[which.max(c(MMDist0$mu[1],MMDist0$mu[2]))]+0.5,
           y = 1,
           label = paste('\'untransfected\' mean =\n', signif(MMDist0$mu[which.max(c(MMDist0$mu[1],MMDist0$mu[2]))], 4))
           , colour = 'violetred4')+
    annotate('label',
             x = mean(as.numeric(unlist(MMexample0)))-0.25,
             y = 1.5,
             label = paste('combined mean\n(+ DMSO) =', signif(mean(as.numeric(unlist(MMexample0))), 4))
             , colour = 'black')+
    annotate('label',
             x = MMDist0$mu[which.min(c(MMDist0$mu[1],MMDist0$mu[2]))]-0.25,
             y = .75,
             label = paste('\'transfected\' mean =\n', signif(MMDist0$mu[which.min(c(MMDist0$mu[1],MMDist0$mu[2]))], 4))
             , colour = 'red')+
    annotate('label',
             x = ko143Dist[[1]]$estimate['mean']+0.5,
             y = 2,
             label = paste('inhibited mean =\n', signif(ko143Dist[[1]]$estimate['mean'], 4))
             , colour = 'darkgreen')+
    ##end####
  ##aesthetic changes:####
  xlab('log10(EtBr fluorescence)') +
    ylab('Population density')+
    ggtitle(paste0('Ko143 inhibition at pH ',pH0))+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
    # xlim(c(0,2.5))+
    # ylim(c(0,1.5))+
    scale_fill_manual(values = c('+ DMSO'='#F8766D','+ 1 uM ko143'='#00BA38' ))+
    scale_color_manual(values=c('inhibited'='forestgreen','untransfected'='violetred4','transfected'='red'))
  ##end####
}


library(gridExtra)
ko143PlotGrid<-grid.arrange(Ko143PlotList[[2]],Ko143PlotList[[1]],nrow=2)
# ggsave(filename=paste0('./Figs/EtBr/linear blanking/ko143 example plots.png'),plot=ko143PlotGrid,width=300,height=380,units='mm')
##end####


####Plot nigericin example plots####
##example Nigericin inhibition plot samples come from 'day 4':
normExamples<-(NormDataList)[names(NormDataList) %like% 'day 4']
normDists<-(NormDistList)[names(NormDistList) %like% 'day 4']
MMexamples<-(MMDataList)[names(MMDataList) %like% 'day 4']
MMDists<-(MMDistList)[names(MMDistList) %like% 'day 4']

NPlotList<-list()

for(j in 1:length(MMexamples)){
  MMexample0<-MMexamples[j]
  MMDist0<-MMDists[[j]]
  if (names(MMexample0) %like% '6'){
    pH0<-'6'
    matchedNormExamples<-normExamples[names(normExamples) %like% '6']
    matchedNormDists<-normDists[names(normDists) %like% '6']
  }else{
    pH0<-'7.4'
    matchedNormExamples<-normExamples[!names(normExamples) %like% '6']
    matchedNormDists<-normDists[!names(normDists) %like% '6']
  }
  
  NigericinExample<-matchedNormExamples[names(matchedNormExamples) %like% '\\+N']
  NigericinDist<-matchedNormDists[names(matchedNormDists) %like% '\\+N']
  
  NPlotList[[pH0]]<-
    
    ##histograms of raw data:####
  ggplot() +
    geom_histogram(
      data = data.frame('x'=unlist(NigericinExample)),
      aes(x,..density.., 
          fill = '+ 10 uM nigericin'),
      binwidth = 0.01,
      alpha = 0.5
    ) +
    geom_histogram(
      data = data.frame('x'=unlist(MMexample0)),
      aes(x,..density.., 
          fill = '+ EtOH'),
      binwidth = 0.01,
      alpha = 0.5
    )+
    ##end####
  ##lines showing fitted normal distribution for inhibited sample:####
  stat_function(
    geom = "line",
    fun = plot_mix_comps,
    args = list(
      mu = NigericinDist[[1]]$estimate['mean'],
      sigma = NigericinDist[[1]]$estimate['sd'],
      lam = 1
    ),
    lwd = 1,
    aes(colour = 'inhibited'),
    show.legend = F
  ) +
    ##end####
  ##lines showing mixed models for active samples:####
  stat_function(
    geom = "line",
    fun = plot_mix_comps,
    args = list(MMDist0$mu[2], MMDist0$sigma[2], lam = MMDist0$lambda[2]),
    lwd = 1,
    aes(colour = 'untransfected')
  ) +
    stat_function(
      geom = "line",
      fun = plot_mix_comps,
      args = list(MMDist0$mu[1], MMDist0$sigma[1], lam = MMDist0$lambda[1]),
      lwd = 1,
      aes(colour = 'transfected')
    )+
    ##end####
  ##Vertical lines indicating various means:####
  geom_vline(xintercept = NigericinDist[[1]]$estimate['mean'], color = 'forestgreen',size=1) +
    geom_vline(xintercept = MMDist0$mu[2], color = 'violetred4',size=1) +
    geom_vline(xintercept = MMDist0$mu[1], color = 'red',size=1) +
    geom_vline(xintercept = mean(as.numeric(unlist(MMexample0))), color = 'black',size=1)+
    ##end####
  ##Label the means:####
  annotate('label',
    x = MMDist0$mu[2]+0.5,
    y = 1,
    label = paste('\'untransfected\' mean =\n', signif(MMDist0$mu[2], 4))
  , colour = 'violetred4')+
    annotate('label',
      x = mean(as.numeric(unlist(MMexample0)))-0.25,
      y = 1.5,
      label = paste('combined mean\n(+ DMSO) =', signif(mean(as.numeric(unlist(MMexample0))), 4))
    , colour = 'black')+
    annotate('label',
      x = MMDist0$mu[1]-0.25,
      y = .75,
      label = paste('\'transfected\' mean =\n', signif(MMDist0$mu[1], 4))
    , colour = 'red')+
    annotate('label',
      x = NigericinDist[[1]]$estimate['mean']+0.5,
      y = 2,
      label = paste('inhibited mean =\n', signif(NigericinDist[[1]]$estimate['mean'], 4))
    , colour = 'darkgreen')+
    ##end####
  ##aesthetic changes:####
  xlab('log10(EtBr fluorescence)') +
    ylab('Population density')+
    ggtitle(paste0('Nigericin inhibition at pH ',pH0))+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
    # xlim(c(0,2.5))+
    # ylim(c(0,1.5))+
    scale_fill_manual(values = c('+ EtOH'='#F8766D','+ 10 uM nigericin'='#00BA38' ))+
    scale_color_manual(values=c('inhibited'='forestgreen','untransfected'='violetred4','transfected'='red'))
  ##end####
}


NigericinPlotGrid<-grid.arrange(NPlotList[[2]],NPlotList[[1]],nrow=2)
# ggsave(filename=paste0('./Figs/EtBr/Nigericin example plots.tiff'),plot=NigericinPlotGrid,width=300,height=380,units='mm')
##end####


###exporting mean values for mixed models:####
names(MMDistList)<- gsub("-k|-N", "no inhibitor", names(MMDistList))

MMparamsMatrix<-matrix(nrow=length(MMDistList),ncol=2)
colnames(MMparamsMatrix)<-c('transfected','untransfected')
rownames(MMparamsMatrix)<-names(MMDistList)

for(k in 1:length(MMDistList)){
  model0<-MMDistList[[k]]
  rownames(MMparamsMatrix)[k]<-names(MMDistList[k])
  means0<-model0$mu
  MMparamsMatrix[k,'untransfected']<-means0[which.max(means0)]
  MMparamsMatrix[k,'transfected']<-means0[which.min(means0)]
}

write.csv(MMparamsMatrix,'EtBr parameters/MMparams_EtBr.csv')

##end####

##exporting combined means: ####
MMCombinedMean<-lapply(MMDataList,mean)
names(MMCombinedMean)<- gsub("-k|-N", "no inhibitor", names(MMCombinedMean))
MMCombinedMean<-as.data.frame(MMCombinedMean)%>%
  t()
colnames(MMCombinedMean)<-'combined mean'
write.csv(MMCombinedMean,'EtBr parameters/MM_combinedMean_EtBr.csv')


samplesNormAll<-sampleNames(FSfinal)[!((pData(FSfinal)$cellType =='WT' &
                                                     pData(FSfinal)$inhibitor == 'active') | pData(FSfinal)$dye =='unstained')]
NormListAll<-list()
for (i in 1:length(samplesNormAll)){

  sample0<-samplesNormAll[i]
  Normdata0 <- FSfinal@frames[[sample0]]@exprs[, c('log_EtBr')] %>%
    as.matrix()%>%
    .[!is.na(.)]
  Norm0 <- fitdistr(Normdata0,'normal')
  NormListAll[[i]]<-Norm0
  names(NormListAll)[i]<-sample0
}
names(NormListAll)<- gsub("\\+k", "plus ko143", names(NormListAll))
names(NormListAll)<- gsub("\\+N", "plus Nigericin", names(NormListAll))

names(NormListAll)<- gsub("-k|-N", "no inhibitor", names(NormListAll))


NormParams<-lapply(NormListAll,function(x)
  x$estimate['mean'])
NormParams<- as.data.frame(NormParams)%>%
  t()
NormParams<-cbind(NormParams,pData(FSfinal)[!((pData(FSfinal)$cellType =='WT' &
                                                 pData(FSfinal)$inhibitor == 'active') | pData(FSfinal)$dye =='unstained'),])
write.csv(NormParams,'EtBr parameters/NormParams_EtBr.csv')

##end####




