rm(list=ls())
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here('reef_functions.R'))

###
fish<- read.csv(here('data','fish.csv'),header=F)
colnames(fish)=c('id','formid','fish_memberid','region','site','speciesid','familyid','abundance','timecode','site4','site1','comment')
R<- subset(fish,region=='TWA')
surveys<- read.csv(here('data','TWAsurveys.csv'),header=F)
colnames(surveys)=c('formid','batch','fish_memberid','exp','type','geogr','geogr4','geogr1','latitude','longitude','date','stemp','btemp','btime','start','visibility','averagedepth','maxdepth','current','habitat')

R[,13:31]<- surveys[match(R$formid,surveys$formid),2:20]

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

goliath<- reef_filter_sp(R=R,GZ=c(2101,2201,2301,3101,3201,3301,3403:3411),sp=96,invert=0)
write.csv(goliath,'Goliath_grouper_all_surveys.csv')
goliath<- subset(goliath,year<2022)
goliath<- subset(goliath,abundance<4)
gg_occs<- goliath

gg_occs$habitat2<- NA
gg_occs<- gg_occs%>%
  mutate(
    habitat2=ifelse(habitat%in%c(0,1,11,12,7,8,9),'mixed',habitat2),
    habitat2=ifelse(habitat%in%c(4,5,6),'dropoff',habitat2),
    habitat2=ifelse(habitat%in%c(10),'artificial',habitat2),
    habitat2=ifelse(habitat%in%c(2),'highreef',habitat2),
    habitat2=ifelse(habitat%in%c(3),'lowreef',habitat2)
  )

write.csv(gg_occs,'Goliath_grouper_2021.csv')