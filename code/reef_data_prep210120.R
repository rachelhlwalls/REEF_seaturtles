rm(list=ls())
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here('Dropbox/#REEF/ExploringData/reef_functions.R'))
# create new project in rstudio and create proj and folder gonna work on that proj with then in that if you put subfolder called data and throw all REEF data in, then the 'data' code works, have data input and output files here(data, TWA, species)
###
fish<- read.csv(here('Dropbox/#REEF/ExploringData/fish.csv'),header=F) # removed 'data', what was that for?
colnames(fish)=c('id','formid','fish_memberid','region','site','speciesid','familyid','abundance','timecode','site4','site1','comment')
R<- subset(fish,region=='TWA')
Turtles<- subset(R,familyid=='111')
surveys<- read.csv(here('Dropbox/#REEF/ExploringData/TropicalWestAtlantic/TWAsurveys.csv'),header=F) # removed 'data', from here('data', 'TWAsurveys.csv') and added Dropbox/#REEF/ExploringData/TropicalWestAtlantic/
colnames(surveys)=c('formid','batch','fish_memberid','exp','type','geogr','geogr4','geogr1','latitude','longitude','date','stemp','btemp','btime','start','visibility','averagedepth','maxdepth','current','habitat')

R[,13:31]<- surveys[match(R$formid,surveys$formid),2:20]

# let's get the year & month from the Date field (requires lubridate)
R$date<-ymd(R$date)#put into proper date format                       # Warning: 3116 failed to parse #
R<-cbind(R,year=year(R$date))
R<-cbind(R,month=month(R$date))
R<-cbind(R,day=day(R$date))

# Thin raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]               # Warning: In `[.data.frame`(R, as.numeric(R$btime) > 20, ) : NAs introduced by coercion #
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]

TWAturtle<- reef_filter_sp(R=R,GZ=c(3101,4102,4201,5705,3408),sp=c(0212,0213,0456:0460),invert=0) # TWA regions, not including Natador depressus, sp. not present in Hawaii or TWA datasets # Error during wrapup: non-numeric argument to binary operator # Error: no more error handlers available (recursive errors?); invoking 'abort' restart

### edited to here ###

write.csv(TWAturtle,'TWAsurveys.csv')
TWAturtle<- subset(TWAturtle,year<2022)
TWAturtle<- subset(TWAturtle,abundance<4)
gg_occs<- TWAturtle

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