rm(list=ls())
library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()

#Load in REEF functions
source(here('code','reef_functions.R'))

#Load in REEF data
fish<- read.csv(here('data','fish.csv'),header=T)
colnames(fish)=c('id','formid','fish_memberid','region','site','speciesid','familyid','abundance','timecode','site4','site1','comment')
R<- subset(fish,region=='TWA')
Turtles<- subset(R,familyid=='111')
R<- subset(fish,region=='HAW')
Turtles<- subset(R,familyid=='066')
R<- subset(fish,region=='CIP')
Turtles<- subset(R,familyid=='110')

R<- subset(fish,region=='EAM')  # 0 obs
Turtles<- subset(R,familyid=='047')  # 0 obs
R<- subset(fish,region=='IORS') # 102,822 obs
Turtles<- subset(R,familyid=='47') # 0 obs

R<- subset(fish,region=='NE') # 36,755 obs
Turtles<- subset(R,familyid=='111') # 0 obs
R<- subset(fish,region=='PAC')
Turtles<- subset(R,familyid=='079')
R<- subset(fish,region=='SAS')
Turtles<- subset(R,familyid=='111')
R<- subset(fish,region=='SOP')
Turtles<- subset(R,familyid=='110')
R<- subset(fish,region=='TEP')
Turtles<- subset(R,familyid=='069')

turtlesTWA<- reef_filter_sp(R=R,GZ=R$site1,sp=c(0655:0661),invert=0)



surveys<- read.csv(here('data','TWAsurveys.csv'),header=F)
colnames(surveys)=c('formid','batch','fish_memberid','exp','type','geogr','geogr4','geogr1','latitude','longitude','date','stemp','btemp','btime','start','visibility','averagedepth','maxdepth','current','habitat')

Turtles[,13:31]<- surveys[match(Turtles$formid,surveys$formid),2:20] # WHAT SHOULD THESE NUMBERS REFER TO EXACTLY
TurtleSurveys<- match(R$formid,surveys$formid) # WHAT SHOULD THESE NUMBERS REFER TO EXACTLY
R$formid %in% surveys$formid

# let's get the year & month from the Date field (requires lubridate)
Turtles$date<-ymd(Turtles$date)#put into proper date format                       # Warning: 3116 failed to parse #
Turtles<-cbind(Turtles,year=year(Turtles$date))
Turtles<-cbind(Turtles,month=month(Turtles$date))
Turtles<-cbind(Turtles,day=day(Turtles$date))

# Thin raw data to exclude dives shorter than 20min, longer than 120min
Turtles<-Turtles[as.numeric(Turtles$btime)>20,]               # Warning: In `[.data.frame`(R, as.numeric(R$btime) > 20, ) : NAs introduced by coercion #
Turtles<-Turtles[as.numeric(Turtles$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
Turtles<-Turtles[as.numeric(Turtles$start)>5,]
Turtles<-Turtles[as.numeric(Turtles$start)<20,]

TWATurtle<- reef_filter_sp(Turtles=Turtles,GZ=c(0655:0660),sp=7,invert=0)





#Sites for Salish Sea
SS_sites<- read.csv(here('data','Dive Site Map Coordinates.csv'),header=T)
CB<- subset(SS_sites,Basin=='CB') #Central Basin sites
HC<- subset(SS_sites,Basin=='HC') #Hood Canal
SOG<- subset(SS_sites,Basin=='SOG') #Straight of Georgia
SPS<- subset(SS_sites,Basin=='SPS') #South Puget Sound
WB<- subset(SS_sites,Basin=='WB') #Whidbey Basin

R<- subset(R_PAC, site %in% SS_sites$Zone.Code) #540k observations
R$basin<- SS_sites$Basin[match(R$site,SS_sites$Zone.Code)] #Codes for the basins of the Salish Sea

R_surv<- R %>% group_by(site) %>% summarize(n.surv=n_distinct(formid))
SS_sites$n.surveys=R_surv$n.surv[match(SS_sites$Zone.Code,R_surv$site)]
SS_sites<- subset(SS_sites,n.surveys>=10)
SS_sites<- subset(SS_sites,is.na(Lat)==F)
write.csv(SS_sites,'SS_survey_sites.csv')
summary(as.factor(R$basin)) #Surveys by Basin

pac_surveys<- read.csv(here('data','PACsurveys.csv')) #Survey-level data
R[,15:34]<- pac_surveys[match(R$formid,pac_surveys$formid),2:21]

# let's get the year & month from the Date field (requires lubridate)
R$date<-ymd(R$date)#put into proper date format
R<-cbind(R,year=year(R$date))
R<-cbind(R,month=month(R$date))
R<-cbind(R,day=day(R$date))

# Thin  raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]

species_freq<- R %>% group_by(speciesid) %>% summarize(n=n())
species_freq$freq<- species_freq$n/length(unique(R$formid))

species_sub<- subset(species_freq,freq>0.01)
species_list<- read.csv(here('data','TWAspecies.csv'))
species_sub[,4:8]<- species_list[match(species_sub$speciesid,species_list$speciesid),2:6]
species_sub<- species_sub[order(species_sub$familyid),]

#Use REEF function for these species
for(i in 1:nrow(species_sub)){
  sp_ss<- reef_filter_sp(R=R,sp=species_sub$speciesid[i],GZ=SS_sites$Zone.Code,invert=species_sub$invert[i])
  sp_cb<- subset(sp_ss,basin=='CB')
  sp_hc<- subset(sp_ss,basin=='HC')
  sp_sog<- subset(sp_ss,basin=='SOG')
  sp_sps<- subset(sp_ss,basin=='SPS')
  sp_wb<- subset(sp_ss,basin=='WB')
  write.csv(sp_ss,here('data','species data','Salish Sea',paste(species_sub$commonname[i],'_ss','.csv',sep='')))
  write.csv(sp_cb,here('data','species data','CB',paste(species_sub$commonname[i],'_cb','.csv',sep='')))
  write.csv(sp_hc,here('data','species data','HC',paste(species_sub$commonname[i],'_hc','.csv',sep='')))
  write.csv(sp_sog,here('data','species data','SOG',paste(species_sub$commonname[i],'_sog','.csv',sep='')))
  write.csv(sp_sps,here('data','species data','SPS',paste(species_sub$commonname[i],'_sps','.csv',sep='')))
  write.csv(sp_wb,here('data','species data','WB',paste(species_sub$commonname[i],'_wb','.csv',sep='')))
}

for(i in 1:nrow(species_sub)){
  sp_ss<- reef_filter_sp(R=R,sp=species_sub$speciesid[i],GZ=SS_sites$Zone.Code,invert=species_sub$invert[i])
  sp_cb<- subset(sp_ss,basin=='CB')
  sp_hc<- subset(sp_ss,basin=='HC')
  sp_sog<- subset(sp_ss,basin=='SOG')
  sp_sps<- subset(sp_ss,basin=='SPS')
  sp_wb<- subset(sp_ss,basin=='WB')
  write.csv(sp_ss,here('data','species data','Salish Sea',paste(species_sub$commonname[i],'_ss','.csv',sep='')))
  write.csv(sp_cb,here('data','species data','CB',paste(species_sub$commonname[i],'_cb','.csv',sep='')))
  write.csv(sp_hc,here('data','species data','HC',paste(species_sub$commonname[i],'_hc','.csv',sep='')))
  write.csv(sp_sog,here('data','species data','SOG',paste(species_sub$commonname[i],'_sog','.csv',sep='')))
  write.csv(sp_sps,here('data','species data','SPS',paste(species_sub$commonname[i],'_sps','.csv',sep='')))
  write.csv(sp_wb,here('data','species data','WB',paste(species_sub$commonname[i],'_wb','.csv',sep='')))
}

#Multi-state models
surfperch<-  reef_filter_sp(R=R,sp=c(116,117,122,123),GZ=SS_sites$Zone.Code,invert=0)
write.csv(surfperch,here('data','species data','multispecies','surfperch.csv'))
rockfish<-  reef_filter_sp(R=R,sp=c(45,48,52,57,58,63,65),GZ=SS_sites$Zone.Code,invert=0)
write.csv(rockfish,here('data','species data','multispecies','rockfish.csv'))
greenling<-  reef_filter_sp(R=R,sp=c(29,30,32,33),GZ=SS_sites$Zone.Code,invert=0)
write.csv(greenling,here('data','species data','multispecies','greenling.csv'))
flatfish<-  reef_filter_sp(R=R,sp=c(24,22,25,26,28,239),GZ=SS_sites$Zone.Code,invert=0)
write.csv(flatfish,here('data','species data','multispecies','flatfish.csv'))
sculpin<-  reef_filter_sp(R=R,sp=c(68,69,72,75,76,78,246,258,262,278,298,304,352),GZ=SS_sites$Zone.Code,invert=0)
write.csv(sculpin,here('data','species data','multispecies','sculpin.csv'))
gunnels<- reef_filter_sp(R=R,sp=c(species_sub$speciesid[species_sub$familyid==51]),GZ=SS_sites$Zone.Code,invert=0)
write.csv(gunnels,here('data','species data','multispecies','gunnels.csv'))
pricklebacks<- reef_filter_sp(R=R,sp=c(species_sub$speciesid[species_sub$familyid==67]),GZ=SS_sites$Zone.Code,invert=0)
write.csv(pricklebacks,here('data','species data','multispecies','pricklebacks.csv'))

