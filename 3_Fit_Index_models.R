###########################################################################
# Load Packages -----------------------------------------------------------
###########################################################################

rm(list = ls())
library(tidyverse)
library(sf)
library(gamm4)
library(hydroGOF)
library(gridExtra)
library(boot)
library(performance)
library(caret)
library(patchwork)
library(DHARMa)
library(mgcViz)
mem.maxVSize(vsize=Inf)

###########################################################################
# Load Data ---------------------------------------------------------------
###########################################################################
#nrwqn and soe data prepared in script 1
load("./Data/invNiwaWq.RData")
load("./Data/niwaMCIWq.RData")
load("./Data/soeMCI.RData")
load("./Data/dfSoe.RData")

#NZ climate zones used for random effects in models
climateZones<-st_read("./Data/ClimateZones/NZ_6regions_EPSG4326.shp")

###########################################################################
# Combine SOE and NRWQN into single frame ---------------------------------
###########################################################################
invNiwaWq%>%
  filter(DRP>0.001,NNN>0.001,NH4N>0.001)%>%
  mutate(rcid="NIWA",
         Source="NIWA",
         sampleID=paste(sid,as.character(date),sep="_"))%>%
  dplyr::select(-NNN,-NH4N)->invNiwaWq

dfSoe%>%
  filter(DRP>0.001,DIN>0.002)%>%
  dplyr::select(-nzreach,-rcsid)%>%
  mutate(Source="SOE")->dfSoe

bind_rows(invNiwaWq,dfSoe)->dfTaxa

#prep mci dfs
niwaMCIWq%>%
  filter(DRP>0.001,NNN>0.001,NH4N>0.001)%>%
  mutate(rcid="NIWA",
         Source="NIWA",
         sampleID=paste(sid,as.character(date),sep="_"))%>%
  dplyr::select(-NNN,-NH4N)->niwaMCIWq

soeMCI%>%
  filter(DRP>0.001,DIN>0.002)%>%
  left_join(unique(dfSoe[,c("sid",
                            "easting",
                            "northing",
                            "Source",
                            "rcid")]))->soeMCI

bind_rows(niwaMCIWq,soeMCI)->dfMetrics

###########################################################################
# Initialise climate zone and VC and other vectors-------------------------
###########################################################################
#convert df to sf object 
dfMetrics%>%
  dplyr::select(nzsegment,easting,northing)%>%
  distinct()%>%
  st_as_sf(coords=c("easting","northing"),crs=27200)%>%
  st_transform(crs=4326)->siteSegments

#join climate zones to segments spatially 
st_join(siteSegments,climateZones%>%
          dplyr::select(Location,geometry),join=st_nearest_feature)%>%
  st_drop_geometry()->siteSegments

#join with dfMetrics
dfMetrics%>%
  left_join(unique(siteSegments))->dfMetrics

dfMetrics%>%
  rename(ClimateZone="Location")%>%
  mutate(ClimateZone=case_when(ClimateZone=="WSI"~"SI West",
                               ClimateZone=="ESI"~"SI East",
                               ClimateZone=="ENI"~"NI East",
                               ClimateZone=="NNI"~"NI North",
                               ClimateZone=="NSI"~"SI North",
                               ClimateZone=="WNI"~"NI West"))->dfMetrics

#VC data exist in two different columns - combining here into one
dfMetrics%>%
  mutate(ClarDisc=case_when(is.na(CLAR)~BDISC,TRUE~CLAR))->dfMetrics

dfMetrics%>%
  mutate(year=lubridate::year(date))%>%
  #exclude >=2020 due to increase in taxa in nrwqn data
  filter(year<2020)%>%
  filter(!is.na(DIN),!is.na(ClarDisc),!is.na(DRP))%>%
  mutate(lnDIN=log10(DIN),lnClar=log10(ClarDisc),lnDRP=log10(DRP))->dfMetrics

#calculate ASPM following NEMs scaling
dfMetrics%>%
  mutate(ASPM=(MCI/200+per_EPT_abun/1+Num_EPT/29)/3)->dfMetrics
  
###########################################################################
# Fit MIMs for MCI and QMCI -----------------------------------------------
###########################################################################
#use gamm4 assuming, gaussian for qmci, mci
source("./Functions/fitGammsUMIMSScaled.R")

dfMetrics%>%
  mutate(sid=as.factor(sid),
         year=as.factor(year),
         ClimateZone=as.factor(ClimateZone))->dfMetricsFold

#scale stressor covariates
dfMetricsFold%>%
  ungroup()%>%
  mutate(rsDIN=(DIN-min(DIN))/(max(DIN)-min(DIN)),
         rsDRP=(DRP-min(DRP))/(max(DRP)-min(DRP)),
         rsClar=(ClarDisc-min(ClarDisc))/
           (max(ClarDisc)-min(ClarDisc)))->dfMetricsScaled

#get same folds as the SSDM from Stoffels et al.
load(paste0("./Data/","Acanthophlebia","_Scaled_StackedModelCV.RData"))

dfMetricsScaled%>%
  left_join(unique(modelCV$Data[,c("sid",".folds")]))->dfMetricsScaled

for (UMIM in c("QMCI","MCI")){
  DF<-dfMetricsScaled
  
  #fit models
  modelCV<-fitGammsUMIMScaled(DF,K=6,umim=UMIM)
  save(modelCV,file=paste0("./Models/UMIM/ModelCV/",UMIM,
                           "_Scaled_ModelCV.RData"))
}

###########################################################################
# Fit ASPM MIM ------------------------------------------------------------
###########################################################################
#use a beta mgcv gam model for aspm
source("./Functions/fitGammsUMIMSScaledBeta.R")
#create vector of year nested within cz to allow nested random smoothers in mgcv
dfMetricsScaled%>%
  mutate(ClimateZone_year=as.factor(paste(ClimateZone,
                                          year,sep=":")))->dfMetricsScaled

for (UMIM in c("ASPM")){
  DF<-dfMetricsScaled
  
  #fit models
  modelCV<-fitGammsUMIMScaledBeta(DF,K=6,umim=UMIM)
  save(modelCV,file=paste0("./Models/UMIM/ModelCV/",UMIM,
                           "_Scaled_Beta_ModelCV.RData"))
}
###########################################################################
# Compile MIM cross-validation statistics ---------------------------------
###########################################################################
#compile model selection information
modelcvByUMIM<-list()
for (UMIM in c("MCI","QMCI","ASPM")){
  
  if (UMIM%in%c("ASPM")){
    load(paste0("./Models/UMIM/ModelCV/",UMIM,"_Scaled_Beta_ModelCV.RData"))
    modelcvByUMIM[[UMIM]]<-modelCV$ModPerformance
  }
  else{
    load(paste0("./Models/UMIM/ModelCV/",UMIM,"_Scaled_ModelCV.RData"))
    modelcvByUMIM[[UMIM]]<-modelCV$ModPerformance
  }
}
save(modelcvByUMIM,file="./Models/UMIM/ModelCV/modelCVStats.RData")

#generate table of umim statistics
modelcvByUMIM%>%
  bind_rows(.id="Metric")%>%
  filter(Metric%in%c("MCI","QMCI","ASPM"))%>%
  mutate(Model=case_when(modelName_fixed=="threeVar"~"DIN+DRP+VC",
                         modelName_fixed=="threeVarInt"~"DIN:DRP:VC",
                         modelName_fixed=="twoVarInt_DIN_DRP"~"DIN:DRP",
                         modelName_fixed=="twoVarInt_clar_DIN"~"DIN:VC",
                         modelName_fixed=="twoVarInt_clar_DRP"~"DRP:VC",
                         modelName_fixed=="twoVar_DIN_DRP"~"DIN+DRP",
                         modelName_fixed=="twoVar_clar_DIN"~"DIN+VC",
                         modelName_fixed=="twoVar_clar_DRP"~"DRP+VC",
                         modelName_fixed=="null"~"NULL",
                         modelName_fixed=="clar"~"VC",
                         TRUE~modelName_fixed),
         Complexity=case_when(Model=="NULL"~1,
                              Model%in%c("DIN","DRP","VC")~2,
                              Model%in%c("DIN+DRP","DIN+VC","DRP+VC")~3,
                              Model%in%c("DIN:DRP","DIN:VC","DRP:VC")~4,
                              Model%in%c("DIN+DRP+VC")~5,
                              Model%in%c("DIN:DRP:VC")~6))%>%
  group_by(Metric)%>%
  mutate(rankAIC=rank(AIC))->modelRanks

modelRanks%>%
  group_by(Metric)%>%
  arrange(rankAIC,.by_group=T)%>%
  dplyr::select(Metric,Model,AIC,wAIC,rmse_CV,nse)%>%
  mutate(AIC=round(AIC,2),
         wAIC=round(wAIC,2),
         rmse_CV=round(rmse_CV,2))->fullUMIMModelSelect

modelRanks%>%
  filter(rankAIC==1)%>%
  mutate(AIC=round(AIC,2),
         wAIC=round(wAIC,2),
         rmse_CV=round(rmse_CV,2),
         nse=round(nse,2))%>%
  dplyr::select(Metric,
                Model,
                AIC,
                wAIC,
                rmse_CV)->modelSelectCVTable

write.csv(modelSelectCVTable,
          file="./Models/UMIM/ModelCV/SmallmodelSelectCVTable.csv",row.names=F)
write.csv(fullUMIMModelSelect,
          file="./Models/UMIM/ModelCV/modelSelectCVTable.csv",row.names=F)

###########################################################################
# Fit best model for each using full data ---------------------------------
###########################################################################
#note, aspm model not fitted due to best model being the null
source("./Functions/fitBestUMIMSScaled.R")

for (umim in c("MCI","QMCI")){
  load(paste0("./Models/UMIM/ModelCV/",umim,"_Scaled_ModelCV.RData"))
  
  bestModel<-fitBestUMIMScaled(cvData=modelCV,metric="AIC",umim=umim)
  
  save(bestModel,file=paste0("./Models/UMIM/BestModel/",
                             umim,
                             "_Scaled_UMIMBestModel.RData"))
}

#aspm model
load(paste0("./Models/UMIM/ModelCV/","ASPM",
            "_Scaled_Beta_ModelCV.RData"))

bestModel<-list()
bestModel[["cvData"]]<-modelCV$Data
bestModel[["BestModel"]]<-gam(ASPM~s(rsDRP,bs="cr",k=6)+
                                s(sid,bs='re')+
                                s(ClimateZone,bs='re')+
                                s(ClimateZone_year,bs='re'),
                              data=modelCV$Data,
                              family=betar(link="logit"),
                              method="REML")

save(bestModel,file=paste0("./Models/UMIM/BestModel/",
                           "ASPM",
                           "_Scaled_UMIMBestModel.RData"))

###########################################################################
# MIM Residual diagnostics ------------------------------------------------
###########################################################################
source("./Functions/plotDharmaResiduals.R")
#######mci
load(paste0("./Models/UMIM/BestModel/","MCI","_Scaled_UMIMBestModel.RData"))
mci<-tibble(res=resid(bestModel$BestModel$mer,type="pearson"),
            drp=(rank(bestModel$Data$DRP)/length(bestModel$Data$DRP)),
            din=(rank(bestModel$Data$DIN)/length(bestModel$Data$DIN)),
            vc=(rank(bestModel$Data$ClarDisc)/length(bestModel$Data$ClarDisc)),
            fitted=rank(fitted(bestModel$BestModel$mer))/
              length(fitted(bestModel$BestModel$mer)))

mci%>%
  ggplot(aes(sample=res))+
  stat_qq()+
  stat_qq_line()+
  ggtitle("MCI Residual QQ plot")+
  ylab("Observed")+
  xlab("Expected")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->mciQQ

mci%>%
  ggplot(aes(x=fitted,y=res))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("MCI Residual versus predicted")+
  ylab("Residual")+
  xlab("Predicted (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->resVpred

mci%>%
  ggplot(aes(x=drp,y=res))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("MCI Residual versus DRP")+
  ylab("Residual")+
  xlab("DRP (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->resVdrp

mci%>%
  ggplot(aes(x=din,y=res))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("MCI Residual versus DIN")+
  ylab("Residual")+
  xlab("DIN (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->resVdin

mci%>%
  ggplot(aes(x=vc,y=res))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("MCI Residual versus VC")+
  ylab("Residual")+
  xlab("VC (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->resVvc

ggsave(filename = './Figures/UMIMAnalysis/ResidualPlots/mciResiduals.tiff',
       plot=mciQQ+resVpred+plot_spacer()+
         resVdin+resVdrp+resVvc+plot_layout(ncol=3), 
       dpi = 300, width = 16, height = 0.7*16, units = 'cm',compression='lzw')

####qmci
load(paste0("./Models/UMIM/BestModel/","QMCI","_Scaled_UMIMBestModel.RData"))

qmci<-tibble(res=resid(bestModel$BestModel$mer,type="pearson"),
             drp=(rank(bestModel$Data$DRP)/length(bestModel$Data$DRP)),
             din=(rank(bestModel$Data$DIN)/length(bestModel$Data$DIN)),
             vc=(rank(bestModel$Data$ClarDisc)/length(bestModel$Data$ClarDisc)),
             fitted=rank(fitted(bestModel$BestModel$mer))/length(fitted(bestModel$BestModel$mer)))
qmci%>%
  ggplot(aes(sample=res))+
  stat_qq()+
  stat_qq_line()+
  ggtitle("QMCI Residual QQ plot")+
  ylab("Observed")+
  xlab("Expected")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->qmciQQ

qmci%>%
  ggplot(aes(x=fitted,y=res))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("QMCI Residual versus predicted")+
  ylab("Residual")+
  xlab("Predicted (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->qresVpred

qmci%>%
  ggplot(aes(x=drp,y=res))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("QMCI Residual versus DRP")+
  ylab("Residual")+
  xlab("DRP (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->qresVdrp

qmci%>%
  ggplot(aes(x=din,y=res))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("QMCI Residual versus DIN")+
  ylab("Residual")+
  xlab("DIN (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->qresVdin

qmci%>%
  ggplot(aes(x=vc,y=res))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("QMCI Residual versus VC")+
  ylab("Residual")+
  xlab("VC (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->qresVvc

ggsave(filename = './Figures/UMIMAnalysis/ResidualPlots/qmciResiduals.tiff',
       plot =qmciQQ+qresVpred+plot_spacer()+
         qresVdin+qresVdrp+qresVvc+plot_layout(ncol=3), 
       dpi = 300, width = 16, height = 0.7*16, units = 'cm',compression='lzw')

####ASPM
load(paste0("./Models/UMIM/BestModel/","ASPM","_Scaled_UMIMBestModel.RData"))
model<-bestModel$BestModel
Data<-bestModel$cvData
Res<-simulateResiduals(fittedModel = model, plot = F,re.form=NULL)

drp<-(rank(Data$DRP)/length(Res$fittedPredictedResponse))
din<-(rank(Data$DIN)/length(Res$fittedPredictedResponse))
vc<-(rank(Data$ClarDisc)/length(Res$fittedPredictedResponse))

#to get qqplot manually
u<-Res$scaledResiduals
n<-length(u)
m<-(1:n)/(n+1)
sx<-sort(m)
sy<-sort(u)
qqDF<-tibble(x=sx,y=sy)

ggplot(data=qqDF,aes(x=sx,y=sy))+
  geom_point()+
  geom_abline()+
  ggtitle("ASPM Residual QQ plot")+
  ylab("Observed")+
  xlab("Expected")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->qq

#to get residuals versus fitted manually
rX<-(rank(Res$fittedPredictedResponse)/length(Res$fittedPredictedResponse))
rrDF<-tibble(x=rX,y=Res$scaledResiduals,DRP=drp,DIN=din,VC=vc)
ggplot(data=rrDF,aes(x=x,y=y))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("ASPM Residual versus predicted")+
  ylab("Residual")+
  xlab("Predicted (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->rr

ggplot(data=rrDF,aes(x=DIN,y=y))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("ASPM Residual versus DIN")+
  ylab("Residual")+
  xlab("DIN (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->rDIN

ggplot(data=rrDF,aes(x=DRP,y=y))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("ASPM Residual versus DRP")+
  ylab("Residual")+
  xlab("DRP (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->rDRP

ggplot(data=rrDF,aes(x=VC,y=y))+
  geom_point(shape=21,alpha=0.5)+
  ggtitle("ASPM Residual versus VC")+
  ylab("Residual")+
  xlab("VC (rank transformed)")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8))->rVC

ggsave(filename = './Figures/UMIMAnalysis/ResidualPlots/aspmResiduals.tiff',
       plot =qq+rr+plot_spacer()+rDIN+rDRP+rVC+plot_layout(ncol=3), 
       dpi = 300, width = 16, height = 0.7*16, units = 'cm',compression='lzw')
