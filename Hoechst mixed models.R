source(
  'C:/Users/Brendan/OneDrive - University Of Cambridge/Data/R analysis/generalFunctions.R'
)
setwd.here()
library(flowCore)
library(tidyverse)
library(data.table)

FSfinal<-readRDS('FS_final_Hoechst.rds')

#the only data sets that exhibit a mixed population are the stained Wt samples without ko143:
samplesForMM<-sampleNames(FSfinal)[sampleNames(FSfinal)%like% 'WT Ho -k']
##include the inhibited equivalent samples for comparison, using simple normal distribution: 
samplesForNorm<-sampleNames(FSfinal)[sampleNames(FSfinal)%like% 'WT Ho \\+k']

library(ggplot2)
library(MASS)
library(mixtools)

set.seed(1811)
MMList<-list()
MMDataList<-list()
for (i in 1:length(samplesForMM)){
  
  sample0<-samplesForMM[i]
  MMdata0 <- FSfinal@frames[[sample0]]@exprs[, c('log_Hoechst')] %>%
    as.matrix()%>%
    .[!is.na(.)]
  MM0 <- normalmixEM(x = MMdata0,
                     fast=FALSE,
                     maxit=10000,
                     epsilon = 1e-16,
                     maxrestarts=1000,
                     k = 2)
  MMList[[i]]<-MM0
  MMDataList[[i]]<-MMdata0
  names(MMList)[i]<-sample0
  names(MMDataList)[i]<-sample0
}

highLam<-list()
lowLam<-list()
nameList<-list()
for (it in 1:length(MMList)){
  whichLam<-which.max(MMList[[it]]$mu)
  
  
  highLam[[it]]<-MMList[[it]]$lambda[whichLam]
  lowLam[[it]]<-MMList[[it]]$lambda[!c(1,2) ==whichLam]
  
  nameList[[it]]<-names(MMList)[it]
}

HoMMlambdas<-data.frame('high'=unlist(highLam),'low'=unlist(lowLam),name=unlist(nameList))
write.csv(HoMMlambdas,('Lambdas/Ho.csv'))


NormList<-list()
NormDataList<-list()
for (i in 1:length(samplesForNorm)){
  
  sample0<-samplesForNorm[i]
  Normdata0 <- FSfinal@frames[[sample0]]@exprs[, c('log_Hoechst')] %>%
    as.matrix()%>%
    .[!is.na(.)]
  Norm0 <- fitdistr(Normdata0,'normal')
  NormList[[i]]<-Norm0
  NormDataList[[i]]<-Normdata0
  names(NormList)[i]<-sample0
  names(NormDataList)[i]<-sample0
}




plot_mix_comps <- function(x, mu, sigma, lam) {
  (lam * dnorm((x), mu, sigma))
}


####Plot####
for(i in 1:3){
  plot0<-ggplot() +
    geom_histogram(
      data = data.frame(x=NormDataList[[i]]),
      aes(x,..density.., 
          fill = '+ 2 uM ko143'),
      binwidth = 0.01,
      # colour = "gray",
      alpha = 0.5
    ) +
    geom_histogram(
      data = data.frame(x=MMDataList[[i]]),
      aes(x, ..density.., fill = '+ DMSO'),
      binwidth = 0.01,
      # colour = "gray",
      alpha = 0.5
    ) +
    xlab('log10(Hoechst fluorescence)') +
    ylab('Population density')+
    stat_function(
      geom = "line",
      fun = plot_mix_comps,
      args = list(
        mu = NormList[[i]]$estimate['mean'],
        sigma = NormList[[i]]$estimate['sd'],
        lam = 1
      ),
      lwd = 1,
      aes(colour = 'inhibited'),
      show.legend = F
    ) +
    stat_function(
      geom = "line",
      fun = plot_mix_comps,
      args = list(MMList[[i]]$mu[2], MMList[[i]]$sigma[2], lam = MMList[[i]]$lambda[2]),
      lwd = 1,
      aes(colour = 'untransfected')
    ) +
    stat_function(
      geom = "line",
      fun = plot_mix_comps,
      args = list(MMList[[i]]$mu[1], MMList[[i]]$sigma[1], lam = MMList[[i]]$lambda[1]),
      lwd = 1,
      aes(colour = 'transfected')
    ) +
    geom_vline(xintercept = NormList[[i]]$estimate['mean'], color = 'darkgreen')+
    geom_vline(xintercept = MMList[[i]]$mu[2], color = 'violetred4') +
    geom_vline(xintercept = MMList[[i]]$mu[1], color = 'red') +
    geom_vline(xintercept = mean(as.numeric(unlist(MMDataList[[i]]))), color = 'black') +
    
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    geom_label(aes(
      x = MMList[[i]]$mu[2] + 0.75,
      y = .75,
      label = paste('\'untransfected\' mean =\n', signif(MMList[[i]]$mu[2], 4))
    ), colour = 'violetred4') +
    geom_label(aes(
      x = mean(as.numeric(unlist(MMDataList[[i]])))-0.5,
      y = 2.5,
      label = paste('combined mean\n(+ DMSO) =', signif(mean(as.numeric(unlist(MMDataList[[i]]))
      ), 4))
    ), colour = 'black') +
    geom_label(aes(
      x = MMList[[i]]$mu[1]-0.5,
      y = 1.5,
      label = paste('\'transfected\' mean =\n', signif(MMList[[i]]$mu[1], 4))
    ), colour = 'red') +
    geom_label(aes(
      x = NormList[[i]]$estimate['mean']+0.75,
      y = 1.5,
      label = paste('\'inhibited\' mean =\n', signif(NormList[[i]]$estimate['mean'], 4))
    ), colour = 'darkgreen') +
    # ylim(c(0,0.0001))+
    xlim(c(1.5,6))+
    scale_fill_manual(values = c('+ DMSO'='#F8766D','+ 2 uM ko143'='#00BA38'))+
    scale_color_manual(values=c('inhibited'='forestgreen','untransfected'='violetred4','transfected'='red'))
  
  ggsave(filename=paste0('./Figs/Hoechst/linear blanking/Ho mixed model graph replicate',i,'.tiff'),plot=plot0,width=300,height=190,units='mm')
}

names(MMList)<- gsub("-k", "no ko143", names(MMList))

MMparamsWT<-lapply(MMList,function(x)
  x$mu)%>%
  as.data.frame

rownames(MMparamsWT)<-c('transfected','untransfected')
write.csv(MMparamsWT,'./Hoechst parameters/MMparams_Hoechst.csv')
MMCombinedMean<-lapply(MMDataList,mean)
names(MMCombinedMean)<- gsub("-k", "no ko143", names(MMCombinedMean))
write.csv(MMCombinedMean,'./Hoechst parameters/MM_combinedMean_Hoechst.csv')


samplesNormAll<-sampleNames(FSfinal)[!(sampleNames(FSfinal)%like% 'WT Ho -k' | sampleNames(FSfinal)%like% 'unstained' )]
NormListAll<-list()
for (i in 1:length(samplesNormAll)){
  
  sample0<-samplesNormAll[i]
  Normdata0 <- FSfinal@frames[[sample0]]@exprs[, c('log_Hoechst')] %>%
    as.matrix()%>%
    .[!is.na(.)]
  Norm0 <- fitdistr(Normdata0,'normal')
  NormListAll[[i]]<-Norm0
  names(NormListAll)[i]<-sample0
}
names(NormListAll)<- gsub("\\+k", "plus ko143", names(NormListAll))
names(NormListAll)<- gsub("-k", "no ko143", names(NormListAll))

NormParams<-lapply(NormListAll,function(x)
  x$estimate['mean'])
NormParams<- as.data.frame(NormParams)

write.csv(NormParams,'./Hoechst parameters/NormParams_Hoechst.csv')



