###########################################################################
# Load packages -----------------------------------------------------------
###########################################################################
rm(list = ls())
library(tidyverse)
library(httr)
library(jsonlite)
library(lubridate)

options(tibble.width=Inf)                                     
options(tibble.print_max=250)

###########################################################################
# Load data ---------------------------------------------------------------
###########################################################################
#wq data from Dupree et al 2017
DF_nrwqn<-read_csv(file="./Data/Merged datasets_TM_CWH.csv",
                   na = c('','NA','na'),
                   guess_max = 10^5,
                   name_repair = "minimal")%>%
  dplyr::select(1:328)

DF_nrwqn%>%
  dplyr::select(sid:rcid,
                easting:nzreach,
                DIN_median_12,
                DRP_median_12,
                CLAR_median_12,
                TURB_median_12,
                sdate.inverts,
                Acarina:Polychaeta)%>%
  gather(species,count,Acarina:Polychaeta)%>%
  mutate(date=dmy(sdate.inverts),
         sampleID=paste0(sid,"_",date),
         present=as.integer(count>0))%>%
  mutate(year=as.numeric(format(date,'%Y')))%>%
  rename(DIN="DIN_median_12",
         DRP="DRP_median_12",
         CLAR="CLAR_median_12",
         TURB="TURB_median_12")%>%
  dplyr::select(-sdate.inverts,-count,-sdate.peri)%>%
  mutate(DIN=DIN/1000,
         DRP=DRP/1000)->oldNiwa

#load tidied nrwqn data from script 1
load("./Data/invNiwa.RData")
load("./Data/niwaMCI.RData")

#WQ data downloaded from hydrowebportal https://niwa.co.nz/niwa-hydro-web-portal
load(file="./Data/nrwqnDataRaw.RData")

###########################################################################
# Initial wrangling -------------------------------------------------------
###########################################################################
bind_rows(wqData)->wqData

#tidy into long format
wqData%>%
  dplyr::select(timestamp,value,sid,parameter)%>%
  pivot_wider(names_from="parameter",values_from="value")%>%
  #calculate DIN
  rename(DRP="DRP-P",
         TURB=`Turbidity (Nephelom)`,
         CLAR=`Visual Water Clarity.Black Disc`,
         NH4N=`NH4-N (Dis)`,
         NNN=`NO3+NO2 as N`)%>%
  #apply detection limit corrections
  mutate(DRP=case_when(DRP<=1~1,TRUE~DRP),
         NH4N=case_when(NH4N<=1~1,TRUE~NH4N),
         TURB=case_when(TURB<=0.05~0.05,TRUE~TURB),
         CLAR=case_when(CLAR<=0.01~0.01,TRUE~CLAR),
         NNN=case_when(NNN<=1~1,TRUE~NNN))%>%
  mutate(DIN=NNN+NH4N)%>%
  dplyr::select(timestamp,sid,DRP,TURB,CLAR,DIN,NNN,NH4N)%>%
  #correct units
  mutate(DRP=DRP/1000,
         DIN=DIN/1000)%>%
  pivot_longer(cols=c("DRP","TURB","CLAR","DIN","NNN","NH4N"))->wqData

###########################################################################
# Calculate 12-monthly medians for samples --------------------------------
###########################################################################
invNiwa%>%
  mutate(date=as.Date(date))%>%
  distinct(sid,date)->wqMedians

sites<-unique(wqMedians$sid)

for (site in sites){
  print(paste(which(sites==site)/length(sites)*100,"% sites completed"))
  samples<-wqMedians[wqMedians$sid==site,"date"]
  for (sampleDate in as.character(samples$date)){
    wqData%>%
      filter(sid==site,
             !timestamp>as.POSIXct(sampleDate,tz="UTC"),
             !timestamp<as.POSIXct(sampleDate,tz="UTC")-31536000)->wq
    wq%>%
      group_by(name)%>%
      dplyr::summarise(md=median(value,na.rm=T),
                       months=n_distinct(timestamp))->wqMd
    
    wqMedians[wqMedians$sid==site & wqMedians$date==as.Date(sampleDate),
              "DIN"]<-wqMd[wqMd$name=="DIN",]$md
    wqMedians[wqMedians$sid==site & wqMedians$date==as.Date(sampleDate),
              "CLAR"]<-wqMd[wqMd$name=="CLAR",]$md
    wqMedians[wqMedians$sid==site & wqMedians$date==as.Date(sampleDate),
              "DRP"]<-wqMd[wqMd$name=="DRP",]$md
    wqMedians[wqMedians$sid==site & wqMedians$date==as.Date(sampleDate),
              "TURB"]<-wqMd[wqMd$name=="TURB",]$md
    wqMedians[wqMedians$sid==site & wqMedians$date==as.Date(sampleDate),
              "NNN"]<-wqMd[wqMd$name=="NNN",]$md
    wqMedians[wqMedians$sid==site & wqMedians$date==as.Date(sampleDate),
              "NH4N"]<-wqMd[wqMd$name=="NH4N",]$md
  }
}

###########################################################################
# Join wq to invert data and save -----------------------------------------
###########################################################################
invNiwa%>%
  left_join(wqMedians)->invNiwaWq

niwaMCI%>%
  left_join(wqMedians)->niwaMCIWq

save(niwaMCIWq,file="./Data/niwaMCIWq.RData")
save(invNiwaWq,file="./Data/invNiwaWq.RData")

