###########################################################################
# Load packages/data ------------------------------------------------------
###########################################################################
rm(list = ls())
library(tidyverse)
library(gamm4)
library(patchwork)
library(boot)
library(trend)

#load nems mci tolerance score table---------------------------------------
mciTolerance<-read_csv(file="./Data/nems_mci_2022_liz.csv")%>%
  mutate(HB=as.numeric(str_remove(HB,"-")),
         SB=as.numeric(str_remove(SB,"-")))%>%
  #note this column is named species in keeping with a previous analysis,
  #but does not specify taxa at species level
  rename(species="MCILevel")%>%
  mutate(species=case_when(species=="Physa = Physella"~"Physa",
                           species=="Hyridella = Echyridella"~"Echyridella",
                           species=="Melanopsis = Zemelanopsis"~"Zemelanopsis",
                           species=="Gundlachia = Ferrissia"~"Ferrissia",
                           species=="Glyptophysa = Physastra"~"Glyptophysa",
                           species=="Hydropsyche - Aoteapsyche"~"Aoteapsyche",
                           species=="Hydropsyche - Orthopsyche"~"Orthopsyche",
                           TRUE~species))

mem.maxVSize(vsize=Inf)

options(tibble.width=Inf)                                     
options(tibble.print_max=250)

###########################################################################
# Test sensitivity of MCIs to stressors -----------------------------------
###########################################################################
source("./Functions/plotGammsUMIMSScaledStates.R")

#load best mci model and plot partial effect of each stressor
load(paste0("./Models/UMIM/BestModel/","MCI","_Scaled_UMIMBestModel.RData"))
mciMod<-bestModel
mci<-plotUMIMS(model.gamm=mciMod$BestModel,Data=mciMod$Data,umim="MCI",
               save=FALSE)

#load best qmci model and plot partial effect of each stressor#qmci
load(paste0("./Models/UMIM/BestModel/","QMCI","_Scaled_UMIMBestModel.RData"))
qmciMod<-bestModel

qmci<-plotUMIMS(model.gamm=qmciMod$BestModel,Data=qmciMod$Data,umim="QMCI",
                save=FALSE)

#aspm
load(paste0("./Models/UMIM/BestModel/","ASPM","_Scaled_UMIMBestModel.RData"))
aspmMod<-bestModel

# Covariate medians:
aspmMod$cvData %>%
  ungroup()%>%
  dplyr::summarise(DIN_med = median(rsDIN,na.rm=T),
                   DRP_med = median(rsDRP,na.rm=T),
                   CLAR_med = median(rsClar,na.rm=T)) -> Medians
# Covariate domains:
aspmMod$cvData %>%
  ungroup()%>%
  dplyr::summarise(DIN_min = min(rsDIN,na.rm=T),
                   DRP_min = min(rsDRP,na.rm=T),
                   CLAR_min = min(rsClar,na.rm=T),
                   DIN_Q05 = quantile(rsDIN,probs=0.025,na.rm=T),
                   DRP_Q05 = quantile(rsDRP,probs=0.025,na.rm=T),
                   CLAR_Q05 = quantile(rsClar,probs=0.025,na.rm=T),
                   DIN_Q95 = quantile(rsDIN,probs=0.975,na.rm=T),
                   DRP_Q95 = quantile(rsDRP,probs=0.975,na.rm=T),
                   CLAR_Q95 = quantile(rsClar,probs=0.975,na.rm=T),
                   DIN_max = max(rsDIN,na.rm=T),
                   DRP_max = max(rsDRP,na.rm=T),
                   CLAR_max = max(rsClar,na.rm=T)) -> Domains

rsDRP <- data.frame(rsDIN = rep(Medians$DIN_med,50),
                           rsDRP = seq(Domains$DRP_Q05,Domains$DRP_Q95,length.out=50),
                           rsClar = rep(Medians$CLAR_med,50),
                    ClimateZone=NA,
                    sid=NA,
                    ClimateZone_year=NA)

rsDRP$y<-predict.gam(aspmMod$BestModel, newdata = rsDRP,
                       exclude=c("s(ClimateZone","s(sid)","s(ClimateZone_year)"),
                       newdata.guaranteed=TRUE)

rsDRP$se<-predict.gam(aspmMod$BestModel, newdata = rsDRP,
                     exclude=c("s(ClimateZone","s(sid)","s(ClimateZone_year)"),
                     newdata.guaranteed=TRUE,se.fit=T)$se.fit

rsDRP%>%
  mutate(upper=inv.logit(y+1.96*se),
         lower=inv.logit(y-1.96*se),
         ASPM=inv.logit(y))->rsDRP

xname<-"DRP (mg/L)"
rsDRP%>%
  mutate(rsDRP=rsDRP*(max(aspmMod$cvData$DRP,na.rm=T)-min(aspmMod$cvData$DRP,na.rm=T))+min(aspmMod$cvData$DRP,na.rm=T))->rsDRP

aspmMod$cvData%>%
  filter(rsDRP>=Domains$DRP_Q05,
         rsDRP<=Domains$DRP_Q95)->rugDF

rugDF%>%
  ungroup()%>%
  mutate(rsDIN=rsDIN*(max(aspmMod$cvData$DIN,na.rm=T)-min(aspmMod$cvData$DIN,na.rm=T))+min(aspmMod$cvData$DIN,na.rm=T),
         rsDRP=rsDRP*(max(aspmMod$cvData$DRP,na.rm=T)-min(aspmMod$cvData$DRP,na.rm=T))+min(aspmMod$cvData$DRP,na.rm=T),
         rsClar=rsClar*(max(aspmMod$cvData$ClarDisc,na.rm=T)-min(aspmMod$cvData$ClarDisc,na.rm=T))+min(aspmMod$cvData$ClarDisc,na.rm=T))->rugDF

states<-tibble(State=c("A","B","C","D"),ymn=c(0.6,0.4,0.3,-Inf),ymx=c(Inf,0.599999,0.399999,0.299999),xmn=rep(-Inf,4),xmx=rep(Inf,4))

ggplot() +
  geom_rect(data=states,aes(fill = State,ymin=ymn,ymax=ymx,xmin=xmn,xmax=xmx),alpha=0.7)+
  geom_ribbon(data=rsDRP,aes(x = rsDRP,ymin = lower,ymax=upper),fill="grey80") +
  geom_line(data = rsDRP, mapping = aes(x = rsDRP, y = ASPM)) +
  geom_rug(data = rugDF, aes(x = rsDRP, y = ASPM), sides = "b")+
  labs(
    y="ASPM",
    x=xname
  )+
  scale_fill_brewer(name="NPSFM state",palette="YlOrRd")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7, face = 'bold'),
        legend.key.width = unit(10,'mm'),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.direction = 'horizontal',
        legend.title.align = 0,
        legend.box = 'vertical',
        legend.box.just = 'right',
        panel.grid.minor=element_blank(),
        legend.margin = margin(0, 0, 0, 0),legend.position = 'right',legend.justification.right="right")->aspmFig

#compile panels
b<-mci$plot_rsDRP+ggtitle("b")+theme(legend.position="none",plot.tag=element_text(face="bold",size=7),axis.text.x=element_blank())+xlab(NULL)
c<-mci$plot_rsDIN+ggtitle("c")+theme(legend.position="none",plot.tag=element_text(face="bold",size=7),axis.text.x=element_blank(),axis.text.y=element_blank())+ylab(NULL)+xlab(NULL)
d<-mci$plot_rsClar+ggtitle("d")+theme(legend.position="none",plot.tag=element_text(face="bold",size=7),axis.text.x=element_blank(),axis.text.y=element_blank())+ylab(NULL)+xlab(NULL)+scale_x_reverse()
e<-qmci$plot_rsDRP+ggtitle("e")+theme(legend.position="none",plot.tag=element_text(face="bold",size=7))
f<-mci$plot_rsDIN+ggtitle("f")+theme(legend.position="none",plot.tag=element_text(face="bold",size=7),axis.text.y=element_blank())+ylab(NULL)
g<-mci$plot_rsClar+ggtitle("g")+theme(legend.position="none",plot.tag=element_text(face="bold",size=7),axis.text.y=element_blank())+ylab(NULL)+scale_x_reverse()+xlab("SFS (VC, m)")
a<-aspmFig+ggtitle("a")+theme(plot.tag=element_text(face="bold",size=7),axis.text.x=element_blank())+xlab(NULL)

layout<-"
ABB
CDE
FGH
"

ggsave(filename = './Figures/umimPartialWBands1.tiff',plot =
         a+guide_area()+b+c+d+e+f+g+plot_layout(guides="collect",design=layout)&guides(fill=guide_legend(title.position="top"))&theme(legend.title = element_text(size = 5),
                                                                                                                                      legend.text = element_text(size = 5)), 
       dpi = 600, width = 11, height = 11, units = 'cm',compression='lzw')

###########################################################################
# Do taxon stressor tolerances correlate with MCI tolerance scores? -------
###########################################################################
#load taxon ttCrit estimates from Stoffels and White 2024
ttCrit<-read_csv("./Data/ttCritTableStoffelsWhite.csv")
Data<-mciMod$Data
ttCrit%>%
  #select species and median ttcrit columns
  dplyr::select(species,MedianDIN,MedianDRP,MedianVC)%>%
  #correct excel conversion of -Inf to error code
  mutate(MedianVC=as.numeric(case_when(MedianVC=="#NAME?"~"-Inf",
                            TRUE~MedianVC)))%>%
  #join mci tolerance values
  left_join(mciTolerance[,c("species","HB")])%>%
  mutate(MedianDIN=case_when(MedianDIN==Inf~quantile(Data$DIN,probs=0.975),
                             TRUE~MedianDIN),
         MedianDRP=case_when(MedianDRP==Inf~quantile(Data$DRP,probs=0.975),
                             TRUE~MedianDRP),
         MedianVC=case_when(MedianVC==-Inf~quantile(Data$ClarDisc,probs=0.025),
                            TRUE~MedianVC))%>%
  #take inverse of vc ttcrit such that its values scale positively with tolerance
  mutate(MedianVC=-MedianVC)%>%
  mutate(MeanTTCrit=(MedianDIN+MedianDRP+MedianVC)/3)%>%
  #take the inverse MCI tolerance score such that high values indicate insensitive taxa
  mutate(HB=-HB)->ttCrit

ttCrit%>%
  left_join(mciTolerance[,c("species","Order")])%>%
  mutate(Order2=case_when(is.na(Order)~species,TRUE~Order))%>%
  mutate(Group=case_when(Order%in%c("Ephemeroptera","Plecoptera","Trichoptera")~"EPT",TRUE~"Non-EPT"))->ttCrit

ttCrit%>%
  ggplot(aes(y=-HB,x=MedianDIN,colour=Group))+
  geom_point(alpha=0.3)+
  theme_bw()+
  scale_colour_manual(values=c("red","orange"))+
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7, face = 'bold'),
        legend.position = 'top',
        legend.key.width = unit(10,'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.direction = 'horizontal',
        legend.title.align = 0,
        legend.box = 'vertical',
        legend.box.just = 'left')+
  xlab(expression(paste("DIN ",italic(TT[crit])," (mg/L)")))+
  ggtitle("a")+
  scale_y_reverse()+
  ylab("NZMCI tolerance value")->dinSensHB

#drp
ttCrit%>%
  ggplot(aes(y=-HB,x=MedianDRP,colour=Group))+
  geom_point(alpha=0.3)+
  theme_bw()+
  scale_colour_manual(values=c("red","orange"))+
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7, face = 'bold'),
        legend.position = 'top',
        legend.key.width = unit(10,'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.direction = 'horizontal',
        legend.title.align = 0,
        legend.box = 'vertical',
        legend.box.just = 'left')+
  xlab(expression(paste("DRP ",italic(TT[crit])," (mg/L)")))+
  ggtitle("b")+
  scale_y_reverse()+
  ylab(NULL)->drpSensHB

#vc
ttCrit%>%
ggplot(aes(y=-HB,x=-MedianVC,colour=Group))+
  geom_point(alpha=0.3)+
  theme_bw()+
  scale_colour_manual(values=c("red","orange"))+
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7, face = 'bold'),
        legend.position = 'top',
        legend.key.width = unit(10,'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.direction = 'horizontal',
        legend.title.align = 0,
        legend.box = 'vertical',
        legend.box.just = 'left')+
  xlab(expression(paste("SFS ",italic(TT[crit])," (VC, m)")))+
  ggtitle("c")+
  scale_x_reverse()+
  scale_y_reverse()+
  ylab(NULL)->vcSensHB

#plot mean xc95~tolerance scores
ttCrit%>%
  ggplot(aes(x=MeanTTCrit,y=-HB,colour=Group))+
  geom_point(alpha=0.3)+
  theme_bw()+
  scale_colour_manual(values=c("red","orange"))+
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7, face = 'bold'),
        legend.position = 'top',
        legend.key.width = unit(10,'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.direction = 'horizontal',
        legend.title.align = 0,
        legend.box = 'vertical',
        legend.box.just = 'left')+
  xlab(expression(paste("Mean stressor ",italic(TT[crit]))))+
  ggtitle("d")+
  scale_y_reverse()+
  ylab(NULL)->meanXC95vTV


ggsave(filename = 
         paste0('./Figures/ssdVsMCITolerance.tiff'),
       plot =dinSensHB+
         drpSensHB+
         vcSensHB+
         meanXC95vTV+plot_layout(ncol=4,guides="collect")&theme(legend.position="bottom",legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10)),
       dpi = 600, width = 16, height = 5, units = 'cm',compression='lzw')

#3 quantify correlations between tolerance scores and xc95
dinCorrelation<-cor.test(x=ttCrit$MedianDIN,y=ttCrit$HB,method="kendall")
drpCorrelation<-cor.test(x=ttCrit$MedianDRP,y=ttCrit$HB,method="kendall")
vcCorrelation<-cor.test(x=ttCrit$MedianVC,y=ttCrit$HB,method="kendall")
meanTTCritCorrelation<-cor.test(x=ttCrit$MeanTTCrit,y=ttCrit$HB,method="kendall")

###########################################################################
# Are species stressor tolerances positively correlated? ------------------
###########################################################################
ttCrit%>%
  ggplot(aes(x=MedianDRP,y=MedianDIN,colour=Group))+
  geom_point(alpha=0.3)+
  theme_bw()+
  scale_colour_manual(values=c("red","orange"))+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7, face = 'bold'),
        legend.position = 'top',
        legend.key.width = unit(10,'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.direction = 'horizontal',
        legend.title.align = 0,
        legend.box = 'vertical',
        legend.box.just = 'left')+
  ggtitle("a")+
  xlab(expression(paste("DRP ",italic(TT[crit])," (mg/L)")))+
  ylab(expression(paste("DIN ",italic(TT[crit])," (mg/L)")))->DinDrp

#Din vs VC
ttCrit%>%
  ggplot(aes(x=-MedianVC,y=MedianDIN,colour=Group))+
  geom_point(alpha=0.3)+
  theme_bw()+
  scale_colour_manual(values=c("red","orange"))+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7, face = 'bold'),
        legend.position = 'top',
        legend.key.width = unit(10,'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.direction = 'horizontal',
        legend.title.align = 0,
        legend.box = 'vertical',
        legend.box.just = 'left')+
  ggtitle("b")+
  scale_x_reverse()+
  xlab(expression(paste("SFS ",italic(TT[crit])," (VC, m)")))+
  ylab(expression(paste("DIN ",italic(TT[crit])," (mg/L)")))->DinVC

#DRP vs VC
ttCrit%>%
  ggplot(aes(y=MedianDRP,x=-MedianVC,colour=Group))+
  geom_point(alpha=0.3)+
  theme_bw()+
  scale_colour_manual(values=c("red","orange"))+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7, face = 'bold'),
        legend.position = 'top',
        legend.key.width = unit(10,'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.direction = 'horizontal',
        legend.title.align = 0,
        legend.box = 'vertical',
        legend.box.just = 'left')+
  ggtitle("c")+
  scale_x_reverse()+
  xlab(expression(paste("SFS ",italic(TT[crit])," (VC, m)")))+
  ylab(expression(paste("DRP ",italic(TT[crit])," (mg/L)")))->DrpVC

ggsave(filename=paste0('./Figures/multiStressSensitivityCovariance.tiff'),
       plot=DinDrp+DinVC+DrpVC+plot_layout(guides="collect")&theme(legend.position="bottom",legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10)), dpi = 300, width = 14, 
       height = 5, units = 'cm',compression='lzw')

#correlation coefficients
dinDrpCorr<-cor.test(x=ttCrit$MedianDIN,y=ttCrit$MedianDRP,method="kendall")
drpVCCorr<-cor.test(x=ttCrit$MedianDRP,y=ttCrit$MedianVC,method="kendall")
dinVCCorr<-cor.test(x=ttCrit$MedianVC,y=ttCrit$MedianDIN,method="kendall")

