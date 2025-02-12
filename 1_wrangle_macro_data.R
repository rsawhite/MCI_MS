###########################################################################
# Load Packages -----------------------------------------------------------
###########################################################################
rm(list = ls())
library(tidyverse)
library(readxl)
library(lubridate)

options(tibble.width=Inf)                                     
options(tibble.print_max=200)
###########################################################################
# Load data ---------------------------------------------------------------
###########################################################################
#load NRWQN data and wrangle into dataframe--------------------------------
inverts<-read_excel(path="./Data/NRWQN Inv 1990 to 2022.xlsx",skip=3)

inverts%>%
  pivot_longer(cols=c(3:length(inverts)),
               names_to = "name",
               values_to = "count")->invLong

ntaxa<-length(unique(invLong$Taxa))

sites<-read_excel(path="./Data/NRWQN Inv 1990 to 2022.xlsx",n_max=3)

sites%>%
  dplyr::select(-"Site code",-"...2")%>%
  slice(2)%>%
  unlist(.,use.names=FALSE)->siteNames

dates<-read_excel(path="./Data/NRWQN Inv 1990 to 2022.xlsx",
                  n_max=1,
                  col_types="date")

dates[,1774][[1]]<-as_datetime(x="2018-08-13",format=NULL)

dates%>%
  dplyr::select(-"Site code",-"...2")%>%
  pivot_longer(cols=c(1:1947))->datesLong

invLong%>%
  mutate(sid=rep(siteNames,times=ntaxa),
         date=rep(datesLong$value,times=ntaxa))%>%
  #filter out this site which was never sampled
  filter(!sid=="DN4")->invLong

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

#load SoE data-------------------------------------------------------------
DF_soe<-read_csv(
  file="./Data/d_final1_metrics_and_inverts_MCIlevel_tosend_withClasses.csv",
  na = c('','NA','na'),
  guess_max = 10^5)

DF_soe%>%
  dplyr::select(sampleID:NZREACH,
                rcid:NZMGN,
                Turbidity,
                MCI_hb:pEPTabundEXCLHydropt,
                Clarity,
                BDISC,
                DRP,
                DIN,
                Acanthophlebia:Uropetala)%>%
  gather(species,
         count,
         Acanthophlebia:Uropetala)%>%
  rename(nzreach="NZREACH",
         easting="NZMGE",
         northing="NZMGN",
         sid="rcsid",
         TURB="Turbidity",
         CLAR="Clarity")%>%
  mutate(date=as.Date(date),
         present=as.integer(count>0),
         year=year(date))->soeData

#load an earlier NRQWN dataset - used for wrangling only, not analysis-----
DF_nrwqn<-read_csv(file="./Data/Merged datasets_TM_CWH.csv",
                   na = c('','NA','na'),
                   guess_max = 10^5,
                   name_repair = "minimal")%>%
  #the whole data set is duplicated col wise. Selecting one set
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
  gather(species,
         count,
         Acarina:Polychaeta)%>%
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
         DRP=DRP/1000)->niwaData

###########################################################################
# Correct an issue taxonomic identities in NRQWN data ---------------------
###########################################################################
invLong%>%
  filter(count>0)%>%
  mutate(year=format(date,format="%Y"))%>%
  group_by(Taxa,Group)%>%
  dplyr::summarise(firstYr=min(year))->newTaxaByYr

#there is a sudden increase in the number of new taxa identified after 2020
#this was likely caused by a changing in lab identification protocols rather
#than true immigration of new taxa, especially given the gap in new taxa between
#2013 to 2020. 
newTaxaByYr%>%
  ggplot(aes(x=firstYr))+
  geom_histogram(stat="count")+
  xlab("Year")+
  ylab("New taxa")

invLong%>%
  filter(count>0)%>%
  mutate(year=format(date,format="%Y"))%>%
  group_by(Taxa)%>%
  dplyr::mutate(firstYr=min(year))%>%
  ungroup()%>%
  filter(firstYr>=2020)%>%
  group_by(Taxa)%>%
  dplyr::summarise(sites=paste0(sid,collapse="_"),
                   TaxaFirstYr=min(firstYr))->newTaxaList

invLong%>%
  filter(count>0)%>%
  mutate(year=format(date,format="%Y"))%>%
  mutate(genus_family=word(Taxa,1))%>%
  group_by(genus_family)%>%
  dplyr::summarise(genusFamilyFirstFound=min(year))->genusFirstFound

#most new taxa are those identified to species which were previously identified
#to genus level or lower resolution following MCI standards. Rather than 
#resolving this issue, we removed data for years>=2020 in later steps of 
#of the analysis.
newTaxaList%>%
  mutate(genus_family=word(Taxa,1))%>%
  left_join(genusFirstFound)%>%
  dplyr::select(Taxa,TaxaFirstYr,genusFamilyFirstFound,sites)->newTaxaList


###########################################################################
# Ensure taxonomy is consistent across NRQWN and SOE data -----------------
###########################################################################
#note that the soe data and older nrwqn data were wrangled previously for 
#Depree et al. 2017 for analyses of MCI. Hence the focus here is to ensure all
#newer nrqwn data use the same MCI level taxonomy.

#MCI identifies species to genus or lower resolution, therefore start by 
#identifying matching genera across data sets
invLong%>%
  mutate(genus=word(Taxa,1))%>%
  mutate(species=genus)->invLong

#columns indicating if genus was present in previous nrwqn dataset or soe data
invLong[invLong$genus%in%unique(niwaData$species),"inPreviousNIWAData"]<-"Yes"
invLong[invLong$genus%in%unique(soeData$species),"inPreviousSOEData"]<-"Yes"

#genera that werent in the old datasets
invLong%>%
  filter(is.na(inPreviousNIWAData))%>%
  filter(is.na(inPreviousSOEData))%>%
  distinct(Taxa,species)%>%
  arrange()

#change species names in current nrqwn data to match level used in SoE and 
#previous nrqwn data
invLong[invLong$Taxa=="Sinelobus stanfordi","species"]<-"Tanaidacea"
invLong[invLong$Taxa%in%c("Phreatogammarus sp.",
                          "Phreatogammarus fragilis"),"species"]<-"Amphipoda"
invLong[invLong$Taxa%in%c("Namalycastis tiriteae"),"species"]<-"Polychaeta"
invLong[invLong$Taxa%in%c("Daphniidae"),"species"]<-"Cladocera"
invLong[invLong$Taxa%in%c("Naididae"),"species"]<-"Oligochaeta"
invLong[invLong$Taxa%in%c("Tubificidae"),"species"]<-"Oligochaeta"
invLong[invLong$Taxa%in%c("Phreodrilidae"),"species"]<-"Oligochaeta"
invLong[invLong$Taxa%in%c("Lumbriculidae"),"species"]<-"Oligochaeta"
invLong[invLong$Taxa%in%c("Enchytraeidae"),"species"]<-"Oligochaeta"
invLong[invLong$Taxa%in%c("Eiseniella spp."),"species"]<-"Oligochaeta"
invLong[invLong$Taxa%in%c("Austropeplea tomentosa"),"species"]<-"Lymnaeidae"
invLong[invLong$Taxa%in%c("Saldula sp."),"species"]<-"Saldidae"
invLong[invLong$Taxa%in%c("Paratrichocladius pluriserialis"),
        "species"]<-"Orthocladiinae"
invLong[invLong$Taxa%in%c("Parachironomus cylindricus"),
        "species"]<-"Chironomidae"
invLong[invLong$Taxa%in%c("Orthoclad sp. B 'Tongue'"),
        "species"]<-"Orthocladiinae"
invLong[invLong$Taxa%in%c("Naonella forsythi"),"species"]<-"Chironomidae"
invLong[invLong$Taxa%in%c("Nr. Kaniwhaniwhanus I"),"species"]<-"Kaniwhaniwhanus"
invLong[invLong$Taxa%in%c("Eukiefferiella spp."),"species"]<-"Orthocladiinae"
invLong[invLong$Taxa%in%c("Chironomini"),"species"]<-"Chironomidae"
invLong[invLong$Taxa%in%c("Chironominae"),"species"]<-"Chironomidae"
invLong[invLong$Taxa%in%c("Oeconesus spp."),"species"]<-"Oeconesidae"
invLong[invLong$genus%in%c("Cricotopus"),"species"]<-"Orthocladiinae"
invLong[invLong$genus%in%c("Diamesinae"),"species"]<-"Chironomidae"

#check for remaining differences between new niwa data and old soe data
unique(c(setdiff(invLong$species,
                 soeData$species),
         setdiff(invLong$species,niwaData$species)))

#remove taxa from the nrwqn data that are not an mci level taxon, or are not 
#identified to sufficient mci level
invLong%>%
  filter(!species%in%c("Gripopterygidae",
                       "Leptophlebiidae",
                       "Hydrobiosidae",
                       "Hydroptilidae",
                       "Polycentropodidae",
                       "Blephariceridae",
                       "Plecoptera",
                       "Peleocorhynchidae",
                       "Osmylidae",
                       "Kaniwhaniwhanus",
                       "Pirara"))->invLong

#All taxa in the nrqwn data are assigned to correct mci taxon
setdiff(invLong$species,mciTolerance$species)

#three taxa in soedata do not have an mci taxonomic entry. Ampullariidae do not
#have a score at any level, while Ecnomidae was not identified to sufficient 
#level (genus) to recieve a score.
setdiff(soeData$species,mciTolerance$species)

#resummarise nrwqn data to get observations by species (MCI level) column 
invLong%>%
  group_by(species,sid,date)%>%
  dplyr::summarise(count=sum(count))%>%
  ungroup()->invNiwa

###########################################################################
# Expand absent taxa across datasets --------------------------------------
###########################################################################
#species present in soe but not in nrwqn data
soeSpp<-setdiff(soeData$species,invNiwa$species)

#add 0s for soeSpp taxa in all site and sample combinations of nrwqn data
for (site in unique(invNiwa$sid)){
  samples<-unique(invNiwa[invNiwa$sid==site,"date"])
  expand_grid(species=soeSpp,sid=site,date=samples$date,count=0)
  invNiwa%>%
    add_row(expand_grid(species=soeSpp,
                        sid=site,
                        date=samples$date,
                        count=0))->invNiwa
}

#species present in niwa but not in soe - there are none
niwaSpp<-setdiff(invNiwa$species,soeData$species)

###########################################################################
# Tidy SOE data -----------------------------------------------------------
###########################################################################
DF_soe%>%
  #remove nrwqn site
  filter(!rcid=="NRWQN1")%>%
  #remove samples with unknown dates
  filter(!is.na(date))%>%
  #create unique site name - this is so we can calculate metrics per site. 
  #Cant do this by reach since there are multiple sites per reach 
  mutate(sid=case_when(is.na(rcsid)~paste(NZREACH,rcid,sep="_"),
                       TRUE~paste(rcid,rcsid,sep="_")))%>%
  dplyr::select(sampleID:NZREACH,
                rcid:NZMGN,
                Turbidity,
                MCI_hb,
                EPTrichEXCLHydropt:pEPTabundEXCLHydropt,
                Clarity,
                BDISC,
                DRP,
                DIN,
                Acanthophlebia:Uropetala,
                sid)%>%
  gather(species,count,Acanthophlebia:Uropetala)%>%
  rename(nzreach="NZREACH",
         easting="NZMGE",
         northing="NZMGN",
         TURB="Turbidity",
         CLAR="Clarity")%>%
  mutate(date=as.Date(date),
         present=as.integer(count>0),
         year=year(date))->dfSoe

#remove samples with wonky dates (i.e., year = 0004) 
dfSoe%>%
  mutate(timeDiff=2023-year)%>%
  filter(timeDiff<2000)%>%
  #remove samples with missing coords
  filter(!is.na(easting))%>%
  dplyr::select(-timeDiff)->dfSoe

###########################################################################
# Calculate mci metrics for nrwqn and soe datasets ------------------------
###########################################################################
source("./Functions/myMetricFunctions.R")

invNiwa%>%
  left_join(mciTolerance%>%dplyr::select(species,HB,Order,Family))%>%
  rename(Abundance="count",
         Taxa="species")->niwaMCI

#calculate metrics for niwa data
niwaMCI%>%
  #following nems, hydroptilids are excluded from %ept calculations
  mutate(Order=case_when(Family=="Hydroptilidae"~"notTrichoptera",
                         TRUE~Order))%>%
  group_by(sid,date)%>%
  filter(Abundance!=0)%>%
  do(summarise(.,
               MCI = mci_score_HB(.),
               QMCI = qmci_score_HB(.),
               Num_species = num_species(.),
               Num_EPT = num_EPT(.),
               per_EPT_taxa = per_EPT_taxa(.),
               per_EPT_abun = per_EPT_abun(.)
  ))%>%
  ungroup()%>%
  #calculate ASPM
  do(mutate(.,
            ASPM=getASPM(.)))->niwaMCI

#recalculate mci for soe data to check its good
dfSoe%>%
  left_join(mciTolerance%>%dplyr::select(species,HB,Order,Family))%>%
  rename(Abundance="count",
         Taxa="species")->soeMCI

soeMCI%>%
  mutate(Order=case_when(Family=="Hydroptilidae"~"notTrichoptera",
                         TRUE~Order))%>%
  group_by(sid,date,sampleID)%>%
  filter(Abundance!=0)%>%
  do(summarise(.,
               MCI = mci_score_HB(.),
               QMCI = qmci_score_HB(.),
               Num_species = num_species(.),
               Num_EPT = num_EPT(.),
               per_EPT_taxa = per_EPT_taxa(.),
               per_EPT_abun = per_EPT_abun(.)
  ))%>%
  #calculate ASPM
  ungroup()%>%
  do(mutate(.,
            ASPM=getASPM(.)))%>%
  #left Join to get wq values
  left_join(unique(dfSoe[,c("sid",
                            "date",
                            "sampleID",
                            "TURB","DIN",
                            "DRP",
                            "CLAR",
                            "BDISC")]))->soeMCI

###########################################################################
# Identify DRN segments for each site -------------------------------------
###########################################################################
#segment assignments made using https://shiny.niwa.local/nzsegments/ using 
#coordinates of data
segAssign<-read_csv("./Data/nzsegment_assignments_2022-06-14.csv")

niwaData%>%
  mutate(sid=str_remove(sid, "NAT-"))%>%
  mutate(sid=case_when(sid=="DN10"~sid,TRUE~str_remove(sid,"0")))%>%
  dplyr::select(sid,easting,northing)%>%
  distinct()%>%
  left_join(segAssign%>%dplyr::select(easting,northing,nzsegment))%>%
  add_row(sid="TU1",nzsegment=7154460,easting=2699800,northing=6249001)->recKey

invNiwa%>%
  left_join(recKey)->invNiwa

niwaMCI%>%
  left_join(recKey)->niwaMCI

#get nzsegments for soe Data
dfSoe%>%
  left_join(segAssign%>%dplyr::select(easting,northing,nzsegment))%>%
  dplyr::select(-MCI_hb,
                -EPTrichEXCLHydropt,
                -pEPTrichEXCLHydropt,
                -pEPTabundEXCLHydropt,
                -present,
                -year)->dfSoe

soeMCI%>%
  left_join(unique(dfSoe[,c("sampleID","sid","date","nzsegment")]))->soeMCI

###########################################################################
# Partition sites by region -----------------------------------------------
###########################################################################
#load REC data
load("./Data/REC24_REC2Utility_RiverNames_HydroPreds_CorrectedJune2018.RData")
REC2Utility%>%
  dplyr::select(nzsegment,BestNZReach,rcID)%>%
  mutate(Region=case_when(rcID==1~"Northland",
                          rcID==2~"Auckland",
                          rcID==3~"Waikato",
                          rcID==4~"Bay of Plenty",
                          rcID==5~"Gisborne",
                          rcID==6~"Taranaki",
                          rcID==7~"Manawatu",
                          rcID==8~"Hawkes Bay",
                          rcID==9~"Wellington",
                          rcID==10~"Tasman",
                          rcID==11~"Marlborough",
                          rcID==12~"West Coast",
                          rcID==13~"Canterbury",
                          rcID==14~"Otago",
                          rcID==15~"Southland",
                          TRUE~"DontKnow"))%>%
  rename(nzreach=BestNZReach)->REC2Utility

invNiwa%>%
  left_join(REC2Utility[,c("nzsegment","Region")])->invNiwa

niwaMCI%>%
  left_join(REC2Utility[,c("nzsegment","Region")])->niwaMCI

dfSoe%>%
  left_join(REC2Utility[,c("nzsegment","Region")])->dfSoe

soeMCI%>%
  left_join(REC2Utility[,c("nzsegment","Region")])->soeMCI

###########################################################################
# Save dataframes ---------------------------------------------------------
###########################################################################
save(niwaMCI,file="./Data/niwaMCI.RData")
save(invNiwa,file="./Data/invNiwa.RData")
save(dfSoe,file="./Data/dfSoe.RData")
save(soeMCI,file="./Data/soeMCI.RData")