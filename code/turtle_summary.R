rm(list=ls())
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here('code', 'reef_functions.R'))
#source(here('code', 'reef_functions2.R')) # >1 species in reef_filter_sp at a time

### Tropical West Atlantic ###
fish<- read.csv(here('data','fish.csv'),header=T)
#colnames(fish)=c('id','formid','fish_memberid','region','site','speciesid','familyid','abundance','timecode','site4','site1','comment')
R<- subset(fish,region=='TWA')
#R<- subset(R,speciesid=='655') # Unidentified spp only
#R<- subset(R,speciesid=='662') # FP greens only
surveys<- read.csv(here('data','TWAsurveys.csv'),header=T)
#colnames(surveys)=c('formid','batch','fish_memberid','exp','type','geogr','geogr4','geogr1','latitude','longitude','date','stemp','btemp','btime','start','visibility','averagedepth','maxdepth','current','habitat')

R[,13:31]<- surveys[match(R$formid,surveys$formid),2:20] # this adds new columns pairing them to formid

# let's get the year & month from the Date field (requires lubridate)
R$date<-lubridate::ymd(R$date)#put into proper date format
R<-cbind(R,year=lubridate::year(R$date))
R<-cbind(R,month=lubridate::month(R$date))
R<-cbind(R,day=lubridate::day(R$date))
R<- R[complete.cases(R$date),]

# Thin  raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]

# Thin data to remove any sightings from before 2002 (turtles only included from summer 2001 on)
R<-R[as.numeric(R$year)>2001,]

# Geog zones
reef_geog<- read.csv(here('data','TWAgeog.csv'),na.strings='NULL')
reef_geog$region.id<- substr(reef_geog$geogid, 1,4) #Get the region id (first four digits)
reef_geog<- reef_geog[complete.cases(reef_geog),] #drop out sites with missing geographic coordinates
GZ=unique(reef_geog$region.id) # there are 227
#GZDS=unique(reef_geog$geogid) # there are 8970


# Unidentified turtle spp
spp=list()
spp_ts=list()
for(i in 1:227){
  spp[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=655,invert=0)
  if(nrow(spp[[i]])>0){
    spp_ts[[i]]=ts_reef(spp[[i]])
  }
}


# Loggerhead
logger=list()
logger_ts=list()
for(i in 1:227){
  logger[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=657,invert=0)
  if(nrow(logger[[i]])>0){
    logger_ts[[i]]=ts_reef(logger[[i]])
  }
}

# Green
green_3=list()
green_ts_3=list()
for(i in 1:227){
  green_3[[i]]=reef_filter_sp_3(R=R,GZ=GZ[i],sp=c(661,662),invert=0)
  if(nrow(green_3[[i]])>0){
    green_ts_3[[i]]=ts_reef(green_3[[i]])
  }
}


#Green sea turtles - all TWA
greens_twa=reef_filter_sp_nogz(R,sp=c(661,662))

length(unique(greens_twa$site))

greens_site_n<- greens_twa %>% group_by(site) %>% summarize(n=n())
greens_site4_n<- greens_twa %>% group_by(site4) %>% summarize(n.surveys=n_distinct(formid),n.sites=n_distinct(site))

greens_site4_ny<- greens_twa %>% group_by(site4,year) %>% summarize(n.surveys=n_distinct(formid),n.sites=n_distinct(site))

survs_8503<- subset(greens_site4_ny,site4==8503)
survs_3403<- subset(greens_site4_ny,site4==3403)
plot(n.surveys~year,data=survs_8503,type='n',bty='l',ylim=c(0,max(survs_8503$n.surveys)))
lines(n.surveys~year,data=survs_8503,lwd=2,col='navy')
lines(n.surveys~year,data=survs_3403,lwd=2,col='goldenrod')

hist(greens_site_n$n)


#green_ds[[i]]=length(unique(green[[i]]$site)) # total # of dive sites
#green[[124]]<-green[[124]][as.numeric(green[[124]]$abundance)>0,] # removing 0s/non-occurrences

green[[124]] %>%
  mutate(occurrence = ifelse(abundance==0,0,1)) %>%
  summarise(propOcc = sum(occurrence)/n())

#green[[124]] <- green[[124]] %>%
#  filter(abundance>0)

ggplot(green[[124]], aes(x = year, fill = green[[124]]$occurrence)) +
  geom_bar()

green[[124]] %>%
  summarise(propOcc = sum(occurrence)/n())

# Add column that is presence absence 1,0, for occurrence or not in for loop
for(i in 1:227){
  green[[i]] <- green[[i]] %>%  
    filter(abundance>0)
}

length(green[[1]]$geogr)
dvsites=list()
for(i in 1:45){
  if(nrow(green[[i]])>0){
    dvsites[[i]]=length(unique(green[[i]]$geogr))
  }
}

# Sum the contents of the list to give total # of dive sites across region:
sum(unlist(dvsites)) # in this case, 158

# Plot # of occurrences per year as barplot
counts<-table(R$year)
barplot(counts)


# Hawksbill
hawks=list()
hawks_ts=list()
for(i in 1:227){
  hawks[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=656,invert=0)
  if(nrow(hawks[[i]])>0){
    hawks_ts[[i]]=ts_reef(hawks[[i]])
  }
}


logger_sites=list()
for(i in 1:227){
  logger_sites[[i]]=unique(logger[[i]]$site4)
}

green_sites=list()
for(i in 1:227){
  green_sites[[i]]=unique(green[[i]]$site4)
}

hawks_sites=list()
for(i in 1:227){
  hawks_sites[[i]]=unique(hawks[[i]]$site4)
}

#logger_t.occ=list()
#for(i in 1:227){
#  if(nrow(logger_t.occ[[i]])>0){
#    logger_t.occ[[i]]=sum(logger_ts[[i]]$n.occ)
#  }
#}


hawks8503<- reef_filter_sp(R=R,GZ=c(8503),sp=656,invert=0)
hawks_ts8503=ts_reef(hawks8503)

green8503<- reef_filter_sp(R=R,GZ=c(8503),sp=661,invert=0) 
green_ts8503=ts_reef(green8503)

logger3301<- reef_filter_sp(R=R,GZ=c(3301),sp=657,invert=0)
logger_ts3301=ts_reef(logger3301)

#table(R$site4, R$abundance>0) # this gives a table of all site #s then false = column of total # of 0s and true = column of total # above 0 (1-4)









### Hawaii ###
R<- subset(fish,region=='HAW')
#R<- subset(R,speciesid=='212') # Green w FP
R<- subset(R,speciesid=='455') # Unidentified spp
surveys<- read.csv(here('data','HAWsurveys.csv'),header=T)
#colnames(surveys)=c('formid','batch','fish_memberid','exp','type','geogr','geogr4','geogr1','latitude','longitude','date','stemp','btemp','btime','start','visibility','averagedepth','maxdepth','current','habitat')

R[,13:31]<- surveys[match(R$formid,surveys$formid),2:20] # this adds new columns pairing them to formid

# let's get the year & month from the Date field (requires lubridate)
R$date<-lubridate::ymd(R$date)#put into proper date format
R<-cbind(R,year=lubridate::year(R$date))
R<-cbind(R,month=lubridate::month(R$date))
R<-cbind(R,day=lubridate::day(R$date))
R<- R[complete.cases(R$date),]

# Thin  raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]

# Thin data to remove any sightings from before 2002 (turtles only included from summer 2001 on)
R<-R[as.numeric(R$year)>2001,]


reef_geog<- read.csv(here('data','HAWgeog.csv'),na.strings='NULL')
reef_geog$region.id<- substr(reef_geog$geogid, 1,4) #Get the region id (first four digits)
reef_geog<- reef_geog[complete.cases(reef_geog),] #drop out sites with missing geographic coordinates
GZ=unique(reef_geog$region.id) # there are 45
#GZDS=unique(reef_geog$geogid) # there are 508

# Green
green=list()
green_ts=list()
for(i in 1:45){
  green[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=c(212,213),invert=0)
  if(nrow(green[[i]])>0){
    green_ts[[i]]=ts_reef(green[[i]])
  }
}

length(green[[1]]$geogr)
dvsites=list()
for(i in 1:45){
  if(nrow(green[[i]])>0){
    dvsites[[i]]=length(unique(green[[i]]$geogr))
  }
}

# Sum the contents of the list to give total # of dive sites across region:
sum(unlist(dvsites)) # in this case, 158

# Plot # of occurrences per year as barplot
counts<-table(R$year)
barplot(counts)


# Hawksbill
hawks=list()
hawks_ts=list()
for(i in 1:45){
  hawks[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=456,invert=0)
  if(nrow(hawks[[i]])>0){
    hawks_ts[[i]]=ts_reef(hawks[[i]])
  }
}

length(hawks[[1]]$geogr)
dvsitesHH=list()
for(i in 1:45){
  if(nrow(hawks[[i]])>0){
    dvsitesHH[[i]]=length(unique(hawks[[i]]$geogr))
  }
}

# Sum the contents of the list to give total # of dive sites across region:
sum(unlist(dvsitesHH)) # in this case, 10


#green_sites=list()
#for(i in 1:45){
#  green_sites[[i]]=unique(green[[i]]$site4)
#}
plots=list()
for(i in 1:45){
  if(nrow(green[[i]])>0){
  plots[[i]]<-plot(green_ts[[i]]$year,green_ts[[i]]$n.survs, type="l", xlab="year", ylab="# surveys")
  }
}



### Tropical Eastern Pacific ###
R<- subset(fish,region=='TEP')
#R<- subset(R,speciesid=='463') # Black turtle only
#R<- subset(R,speciesid=='462') # Green w FP only
R<- subset(R,speciesid=='455') # Unidentified spp only
surveys<- read.csv(here('data','TEPsurveys.csv'),header=T)
#colnames(surveys)=c('formid','batch','fish_memberid','exp','type','geogr','geogr4','geogr1','latitude','longitude','date','stemp','btemp','btime','start','visibility','averagedepth','maxdepth','current','habitat')

R[,13:31]<- surveys[match(R$formid,surveys$formid),2:20] # this adds new columns pairing them to formid

# let's get the year & month from the Date field (requires lubridate)
R$date<-lubridate::ymd(R$date)#put into proper date format
R<-cbind(R,year=lubridate::year(R$date))
R<-cbind(R,month=lubridate::month(R$date))
R<-cbind(R,day=lubridate::day(R$date))
R<- R[complete.cases(R$date),]

# Thin  raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]

# Thin data to remove any sightings from before 2002 (turtles only included from summer 2001 on)
R<-R[as.numeric(R$year)>2001,]

reef_geog<- read.csv(here('data','TEPgeog.csv'),na.strings='NULL')
reef_geog$region.id<- substr(reef_geog$geogid, 1,4) #Get the region id (first four digits)
reef_geog<- reef_geog[complete.cases(reef_geog),] #drop out sites with missing geographic coordinates
GZ=unique(reef_geog$region.id) # there are 124


# Green
green=list()
green_ts=list()
for(i in 1:124){
  green[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=c(461,462,463),invert=0)
  if(nrow(green[[i]])>0){
    green_ts[[i]]=ts_reef(green[[i]])
  }
}

black=list()
black_ts=list()
for(i in 1:124){
  black[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=c(463),invert=0)
  if(nrow(black[[i]])>0){
    black_ts[[i]]=ts_reef(black[[i]])
  }
}

green_sites=list()
for(i in 1:124){
  green_sites[[i]]=unique(green[[i]]$site4)
}



# Plot # of occurrences per year as barplot
counts<-table(R$year)
barplot(counts)






### Central Indo Pacific ###
R<- subset(fish,region=='CIP')
#R<- subset(R,familyid=='111')
R<- subset(R,speciesid=='287')
surveys<- read.csv(here('data','CIPsurveys.csv'),header=T)
#colnames(surveys)=c('formid','batch','fish_memberid','exp','type','geogr','geogr4','geogr1','latitude','longitude','date','stemp','btemp','btime','start','visibility','averagedepth','maxdepth','current','habitat')

R[,13:31]<- surveys[match(R$formid,surveys$formid),2:20] # this adds new columns pairing them to formid

# let's get the year & month from the Date field (requires lubridate)
R$date<-lubridate::ymd(R$date)#put into proper date format
R<-cbind(R,year=lubridate::year(R$date))
R<-cbind(R,month=lubridate::month(R$date))
R<-cbind(R,day=lubridate::day(R$date))
R<- R[complete.cases(R$date),]

# Thin  raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]

# Thin data to remove any sightings from before 2002 (turtles only included from summer 2001 on)
R<-R[as.numeric(R$year)>2001,]

reef_geog<- read.csv(here('data','CIPgeog.csv'),na.strings='NULL')
reef_geog$region.id<- substr(reef_geog$geogid, 1,4) #Get the region id (first four digits)
reef_geog<- reef_geog[complete.cases(reef_geog),] #drop out sites with missing geographic coordinates
GZ=unique(reef_geog$region.id) # there are 184


# Green
green=list()
green_ts=list()
for(i in 1:184){
  green[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=285,invert=0)
  if(nrow(green[[i]])>0){
    green_ts[[i]]=ts_reef(green[[i]])
  }
}

# Hawksbill
hawks=list()
hawks_ts=list()
for(i in 1:184){
  hawks[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=286,invert=0)
  if(nrow(hawks[[i]])>0){
    hawks_ts[[i]]=ts_reef(hawks[[i]])
  }
}

# Loggerhead
logger=list()
logger_ts=list()
for(i in 1:184){
  logger[[i]]=reef_filter_sp(R=R,GZ=GZ[i],sp=3069,invert=0)
  if(nrow(hawks[[i]])>0){
    logger_ts[[i]]=ts_reef(logger[[i]])
  }
}


green_sites=list()
for(i in 1:184){
  green_sites[[i]]=unique(green[[i]]$site4)
}

hawks_sites=list()
for(i in 1:184){
  hawks_sites[[i]]=unique(hawks[[i]]$site4)
}


# Plot # of occurrences per year as barplot
counts<-table(R$year)
barplot(counts)




logger<- subset(logger,abundance<4)
logger<- subset(logger,year<2022)
#write.csv(goliath,'Goliath_grouper_all_surveys.csv')
log_occs<- logger

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
###Functions####
reef_filter_sp = function(R,GZ,sp,agg){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% subset(site4 %in% GZ) %>% select('formid','speciesid','abundance',everything())
  TempDat_sp<- TempDat %>% subset(speciesid==sp$speciesid) %>% subset(abundance<4)
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n=n(),n.year=length(unique(year))) %>% subset(n.year>2)
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  Zeros<- complete(TempDat2,formid,nesting(speciesid),
                   fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat2$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat2[m,4:ncol(TempDat2)] #Replicate the survey-level data (except abundance)
  TempDat3<- rbind(TempDat2,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat3$occ=ifelse(TempDat3$abundance>0,1,0) #Code presence/absence based on abundance
  TempDat3<- TempDat3 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat3$abundance2<- TempDat3$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat3$abund_trans <-NA #Transform abundance categories into the minimum counts
  TempDat3<- TempDat3%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
 surveyors<- TempDat3 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
 surveyors_trim<- subset(surveyors,n>=3) #only keep surveys by members with 3 or more dives
  
 # TempDat3<- subset(TempDat2, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveys with less than 5
  site_subsets<-TempDat3 %>% group_by(geogr) %>% summarize(n=n_distinct(formid)) %>% subset(n>=5) #Calculate surveys per site
  
  TempDat4<- TempDat3 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  survs<- TempDat4 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat4<- TempDat4 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  if(agg==T){
    TempDat4<- TempDat4 %>% subset(abundance>2) %>% subset(year>1993)
  }else{
    TempDat4<- TempDat4 %>% subset(abundance<3) %>% subset(year>1993)
  }
  
  
  occ_dat<- subset(TempDat4,speciesid==sp$speciesid) #Subset out each species in the provided dataframe
  return(occ_dat)
}

#ts_reef = function(X,sp){
#    sp_x<-X #Subset out each species in the provided dataframe
#    occ_by_year<- sp_x %>% group_by(year) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n#.surv,sp=sp$commonname[match(unique(sp_x$speciesid),sp$speciesid)])
#    abun_by_year<- sp_x %>% group_by(year,geogr) %>% summarize(site_abun=mean(abund_trans),sd_abund=sd(abund_trans),n.surv=n()) %>% group_by(year) %>% summarize(mean_abund=mean(site_abun),n.survs=sum(n.surv),n.sites=n())
#    total_sd<- sp_x %>% group_by(year) %>% summarize(sd=sd(abund_trans))
    
#    comb<- left_join(occ_by_year,abun_by_year)
#    comb2<- left_join(comb,total_sd)
#    ts_dat<- comb2
#  return(ts_dat)
#}

abund_tranfs<- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*2+probs[4]*11+probs[5]*101
  return(sum)
}

ord_to_n<- function(x,c){
  abund_x<- numeric(length(x))
  p= matrix(ncol=5,nrow=length(x))
  if(ncol(c)==2){
      p[,1]=plogis(c[,1]-x)
      p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
      p[,3]=1-plogis(c[,2]-x)
      p[,4]=0
      p[,5]=0
      for(i in 1:length(x)){
        abund_x[i]=abund_tranfs(p[i,])  
      }
  }
  if(ncol(c)==3){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=1-plogis(c[,3]-x)
    p[,5]=0
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  if(ncol(c)==4){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=plogis(c[,4]-x[,i])-plogis(c[,3]-x[,i])
    p[,5]=1-plogis(c[,4]-x[,i])
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  return(abund_x)
}

TS_stan_state_only_plot_MARSS<- function(sp,GZ,params1,TT,ts){
  pdf(paste('./figures/',paste(sp,GZ,'state_space',sep='_'),'.pdf',sep=''),width=8,height=6)
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,1200))
    
    if(ncol(params1$c)==2){
    
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
        reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
        reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0

      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==3){

        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
        reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i]-params1$a)-plogis(params1$c[,1]-params1$x[,i]-params1$a)
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i]-params1$a)-plogis(params1$c[,2]-params1$x[,i]-params1$a)
        reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,i]-params1$a)
        reef_coef[,10]<- 0
    
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==4){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr2[,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr2[,i])-plogis(params1$c[,1]-params1$a_yr2[,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr2[,i])-plogis(params1$c[,2]-params1$a_yr2[,i])
        reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr2[,i])-plogis(params1$c[,3]-params1$a_yr2[,i])
        reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr2[,i])
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i]-params1$a)
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i]-params1$a)-plogis(params1$c[,1]-params1$x[,i]-params1$a)
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i]-params1$a)-plogis(params1$c[,2]-params1$x[,i]-params1$a)
        reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,i]-params1$a)-plogis(params1$c[,3]-params1$x[,i]-params1$a)
        reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,i]-params1$a)
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
  }  
  
  
  x_mat<- data.frame(median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    x_mat[i,1]=median(lambda_mat[[i]]$lambda.x)
    x_mat[i,2]=quantile(lambda_mat[[i]]$lambda.x,0.05)
    x_mat[i,3]=quantile(lambda_mat[[i]]$lambda.x,0.95)
  }
  
  y_mat<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  
  for(i in 1:TT){
    y_mat[i,2]=median(lambda_mat[[i]]$lambda.y)
    y_mat[i,3]=quantile(lambda_mat[[i]]$lambda.y,0.05)
    y_mat[i,4]=quantile(lambda_mat[[i]]$lambda.y,0.95)
  }
  
  par(xpd=T)
  plot(y_mat$median.reef~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(y_mat[,2]),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col='darkcyan')
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  polygon(x, y1, col = adjustcolor('darkcyan', alpha = 0.1), border=NA) # Add uncertainty polygon
  
  #for(z in 1:26){
  #   lines(c(y_mat[z,6],y_mat[z,7])~rep(y_mat[z,1],2),lwd=1,col=adjustcolor('dodgerblue4',alpha.f=0.6))
  #   lines(c(y_mat[z,3],y_mat[z,4])~rep(y_mat[z,1],2),lwd=1,col=adjustcolor('firebrick4',alpha.f=0.6))
  # }
  lines(y_mat$median.reef~y_mat$year,col='dodgerblue4',lwd=2)
  points(y_mat$median.reef~y_mat$year,col='white',pch=21,bg='dodgerblue4',cex=1.5)
  text(y=rep(0,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  #text(y=1.02*max(c(max(y_mat[,2]),max(x_mat))),x=max(y_mat$year)+1.5,'No. surveys',cex=0.5)
 # text(y=0.95*max(c(max(y_mat[,2]),max(x_mat))),x=ts$year,ts$n.sites,cex=0.7)
  
#  legend(2013,c(max(c(max(y_mat[,2]),max(x_mat)))*1.15),c('Est. RVC surveys','Est. REEF surveys'),text.col=c(adjustcolor('navy',alpha.f=0.5),adjustcolor('darkred',alpha.f=0.5),'dodgerblue4','firebrick4'),bty='n')
  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
  
}


###REEF Data Filtering####

##
gg_occs<- read.csv(here('data','Goliath_Grouper_surveys_fl.csv'))

####Stan model - ordinal abundance####
SS_trend_ord<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
}
parameters {
  ordered[K-1] c; //cutpoints
  real x0; //initial popn size

  //deviations from intercept
  vector[Z] beta; //effort coefficients 
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r;
  real<lower = 0> sd_q;
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr] obs_dev; //observation deviations 
}

transformed parameters{
  vector[TT] x;
  vector[N_yr] a_yr;

  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t];
  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + obs_dev[i]; 
  }
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  x0 ~ normal(0,5); //initial state

  //variance terms
  sd_q ~ inv_gamma(2,0.25);
  sd_r ~ inv_gamma(2,0.25);
  sd_site ~ inv_gamma(2, 1);
  sd_dv ~ inv_gamma(2, 1);
  sd_dmy ~ inv_gamma(2, 1);
  
  //varying intercepts
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev[t] ~ normal(0,sd_r);
  }

  y ~ ordered_logistic(a_yr[year_id]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X*beta,c);
  
}
"

SS_trend_ord1<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=1> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=1> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int Q; //columns in year matrix
  matrix[N,Q] X_yr; // design matrix of years X
  int K; //ordinal levels
  int<lower=1> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
}
parameters {
  ordered[K-1] c; //cutpoints
  
  //deviations from intercept
  vector[Z] beta; //effort coefficients 
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
  vector[N_yr] a_yr; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_yr;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
}
model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  
  //variance terms
  sd_yr ~ inv_gamma(2,0.25);
  sd_dv ~ inv_gamma(2, 1);
  sd_dmy ~ inv_gamma(2, 1);
  
  //varying intercepts
  a_yr ~ normal(0, sd_yr);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);


  y ~ ordered_logistic(a_yr[year_id]+a_dv[diver]+a_dmy[dmy]+X*beta,c);
  
}
"

SS_trend_ord1<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=1> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=1> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=1> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int Q; //columns in year matrix
  matrix[N,Q] X_yr; // design matrix of years X
  int K; //ordinal levels
  int<lower=1> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
}
parameters {
  ordered[K-1] c; //cutpoints
  real b_yr;

  //deviations from intercept
  vector[Z] beta; //effort coefficients 
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
  vector[N_yr] a_yr; //deviation between site day clusters
 
  //varying slope/intercepts for site
  vector[2] b_site[N_site]; //varying intercepts and slopes of year between sites
  vector<lower = 0>[2] sd_site; // standard dev. of intercept and slope for site
  vector[2] mu_site; //intercept and slope hyper-prior
  corr_matrix[2] Omega; //correlation matrix
 
  //variance on the deviance components
  real<lower = 0> sd_yr;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
}
model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  b_yr ~ normal(0,0.5); //global slope
  
  //varying intercept/slope for site
  mu_site ~ normal(0, 3);
  sd_site ~ inv_gamma(2, 1);
  Omega ~ lkj_corr(2);
  b_site ~ multi_normal(mu_site, quad_form_diag(Omega, sd_site));
  
  //variance terms
  sd_yr ~ inv_gamma(2,0.25);
  sd_dv ~ inv_gamma(2, 1);
  sd_dmy ~ inv_gamma(2, 1);
  
  //varying intercepts
  a_yr ~ normal(0, sd_yr);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);

  y ~ ordered_logistic(a_yr[year_id]+X_yr*b_site[site]+a_dv[diver]+a_dmy[dmy]+X*beta,c);

  
}
"

SS_trend_ord1<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=1> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=1> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=1> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int Q; //columns in year matrix
  matrix[N,Q] X_yr; // design matrix of years X
  int K; //ordinal levels
  int<lower=1> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
}
parameters {
  ordered[K-1] c; //cutpoints
  real b_yr;

  //deviations from intercept
  vector[Z] beta; //effort coefficients 
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
  vector[N_yr] a_yr; //deviation between site day clusters
 
  //varying slope/intercepts for site
  vector[2] b_site[N_site]; //varying intercepts and slopes of year between sites
  vector<lower = 0>[2] sd_site; // standard dev. of intercept and slope for site
  vector[2] mu_site; //intercept and slope hyper-prior
  corr_matrix[2] Omega; //correlation matrix
 
  //variance on the deviance components
  real<lower = 0> sd_yr;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
}
model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  b_yr ~ normal(0,0.5); //global slope
  
  //varying intercept/slope for site
  mu_site ~ normal(0, 3);
  sd_site ~ inv_gamma(2, 1);
  Omega ~ lkj_corr(2);
  b_site ~ multi_normal(mu_site, quad_form_diag(Omega, sd_site));
  
  //variance terms
  sd_yr ~ inv_gamma(2,0.25);
  sd_dv ~ inv_gamma(2, 1);
  sd_dmy ~ inv_gamma(2, 1);
  
  //varying intercepts
  a_yr ~ normal(0, sd_yr);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);

  y ~ ordered_logistic(a_yr[year_id]+X_yr*b_site[site]+a_dv[diver]+a_dmy[dmy]+X*beta,c);

  
}
"


site_trend_model1<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int K; //ordinal levels
}
parameters {
  ordered[K-1] c; //cutpoints

  //varying intercept
  vector[Z] beta; //effort coefficients - RVC
 
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
 
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real alpha;
  
}
model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  
  //variance terms
  sd_dv ~ inv_gamma(1, 0.5);
  sd_dmy ~ inv_gamma(2, 1);
  
  //varying intercepts
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(i in 1:N){
    y ~ ordered_logistic(a_dv[diver]+a_dmy[dmy]+X*beta,c);
  }

  
}
"


var_slope_test<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=1> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=1> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  matrix[N,2] X_yr; // design matrix of years X
  int K; //ordinal levels
    // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
}
parameters{
  ordered[K-1] c; //cutpoints

  //deviations from intercept
  vector[Z] beta; //effort coefficients 
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //varying slope/intercepts for site
  vector[2] b;
  vector<lower=0>[M_1] sd_1; 
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
 
  //variance on the deviance components
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
}
transformed parameters{
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  // compute actual group-level effects
  r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)';
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  
  vector[N] mu = X_yr*b;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
  }
}
model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  
  //varying intercept/slope for site
  sd_1 ~ inv_gamma(4, 2);
  L_1 ~ lkj_corr_cholesky(1);
  to_vector(z_1) ~ std_normal();
  
  //variance terms
  sd_dv ~ inv_gamma(4, 2);
  sd_dmy ~ inv_gamma(4, 2);
  
  //varying intercepts
  a_dv ~ std_normal();
  a_dmy ~ std_normal();
  
   y ~ ordered_logistic(mu+a_dv[diver]*sd_dv+a_dmy[dmy]*sd_dmy+X*beta,c);
}
generated quantities {
  vector<lower=-1,upper=1>[NC_1] cor_1;
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
   
   for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
"

##Testing site trend model
gg_occs$year.num<- gg_occs$year-min(gg_occs$year)

test<- brm(abundance2~year.num + (year.num|geogr4/geogr),
              data=gg_occs,
              family=cumulative("logit",threshold = "flexible"), 
              control = list(adapt_delta = 0.98),
              cores=4,
              save_all_pars = T,
              warmup = 100, iter = 300, chains = 4)

brms::make_stancode(abundance2~year.num + (year.num|geogr4/geogr),
                    data=gg_occs,
                    family=cumulative(link="logit",threshold = "flexible"), 
                    control = list(adapt_delta = 0.98),
                    cores=4,
                    save_all_pars = T,
                    warmup = 100, iter = 300, chains = 4)


gg_occs_pre10<- subset(gg_occs,year<2010)
X_yr<- matrix(data=c(rep(1,nrow(gg_occs_pre10)),as.numeric(factor(gg_occs_pre10$year))-1),ncol=2,nrow=nrow(gg_occs_pre10))
X<- matrix(data=c(scale(as.numeric(gg_occs_pre10$btime)),scale(as.numeric(gg_occs_pre10$averagedepth)),scale(as.numeric(gg_occs_pre10$visibility)),scale(as.numeric(gg_occs_pre10$current)),gg_occs_pre10$exp_binary),ncol=5,nrow=nrow(gg_occs_pre10))

gg_slope_pre2010<- rstan::stan(model_code = var_slope_test, data = list(y =gg_occs_pre10$abundance2,
                                                           N = nrow(gg_occs_pre10),
                                                           diver=as.numeric(factor(gg_occs_pre10$fish_memberid)),
                                                           N_dv=length(unique(gg_occs_pre10$fish_memberid)),
                                                           dmy=as.numeric(factor(gg_occs_pre10$site_dmy)),
                                                           N_dmy=length(unique(gg_occs_pre10$site_dmy)),
                                                           K=length(unique(gg_occs_pre10$abundance2)),
                                                           X=X,
                                                           X_yr=X_yr,
                                                           Z=ncol(X),
                                                           N_1=length(unique(gg_occs_pre10$geogr)),
                                                           M_1=2,
                                                           J_1=as.numeric(factor(gg_occs_pre10$geogr)),
                                                           Z_1_1=X_yr[,1],
                                                           Z_1_2=X_yr[,2],
                                                           NC_1=1),
                    pars = c('c','sd_dv','sd_dmy','b','r_1','sd_1','L_1'),
                    control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

params_pre2010<- rstan::extract(gg_slope_pre2010)

gg_occs_post10<- subset(gg_occs,year>2009)
X_yr<- matrix(data=c(rep(1,nrow(gg_occs_post10)),as.numeric(factor(gg_occs_post10$year))-1),ncol=2,nrow=nrow(gg_occs_post10))
X<- matrix(data=c(scale(as.numeric(gg_occs_post10$btime)),scale(as.numeric(gg_occs_post10$averagedepth)),scale(as.numeric(gg_occs_post10$visibility)),scale(as.numeric(gg_occs_post10$current)),gg_occs_post10$exp_binary),ncol=5,nrow=nrow(gg_occs_post10))

gg_slope_post2010<- rstan::stan(model_code = var_slope_test, data = list(y =gg_occs_post10$abundance2,
                                                                        N = nrow(gg_occs_post10),
                                                                        diver=as.numeric(factor(gg_occs_post10$fish_memberid)),
                                                                        N_dv=length(unique(gg_occs_post10$fish_memberid)),
                                                                        dmy=as.numeric(factor(gg_occs_post10$site_dmy)),
                                                                        N_dmy=length(unique(gg_occs_post10$site_dmy)),
                                                                        K=length(unique(gg_occs_post10$abundance2)),
                                                                        X=X,
                                                                        X_yr=X_yr,
                                                                        Z=ncol(X),
                                                                        N_1=length(unique(gg_occs_post10$geogr)),
                                                                        M_1=2,
                                                                        J_1=as.numeric(factor(gg_occs_post10$geogr)),
                                                                        Z_1_1=X_yr[,1],
                                                                        Z_1_2=X_yr[,2],
                                                                        NC_1=1),
                               pars = c('c','sd_dv','sd_dmy','b','r_1','sd_1','L_1','cor_1'),
                               control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)


params_post2010<- rstan::extract(gg_slope_post2010)



###Florida (All Regions) Population Trajectory####
X<- matrix(data=c(scale(as.numeric(gg_occs$btime)),scale(as.numeric(gg_occs$averagedepth)),scale(as.numeric(gg_occs$visibility)),scale(as.numeric(gg_occs$current)),gg_occs$exp_binary),ncol=5,nrow=nrow(gg_occs))

gg_SS<- rstan::stan(model_code = SS_trend_ord, data = list(y =gg_occs$abundance2,
                                                                N = nrow(gg_occs),
                                                                site=as.numeric(factor(gg_occs$geogr)),
                                                                N_site=length(unique(gg_occs$geogr)),
                                                                diver=as.numeric(factor(gg_occs$fish_memberid)),
                                                                N_dv=length(unique(gg_occs$fish_memberid)),
                                                                dmy=as.numeric(factor(gg_occs$site_dmy)),
                                                                N_dmy=length(unique(gg_occs$site_dmy)),
                                                                K=length(unique(gg_occs$abundance2)),
                                                                X=X,
                                                                Z=ncol(X),
                                                                TT=27,
                                                                N_yr=length(unique(gg_occs$year)),
                                                                yr_index=sort(unique(as.numeric(factor(gg_occs$year)))),
                                                                year_id=as.numeric(factor(gg_occs$year))),
                   pars = c('c','sd_site','sd_dv','sd_dmy','sd_r','sd_q','a_yr','beta'),
                   control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)


gg_params<- rstan::extract(gg_SS)

diff=(ord_to_n(gg_params$x[,27],gg_params$c)/ord_to_n(gg_params$x[,17],gg_params$c))
median(diff)-1
quantile(diff,0.05)-1
quantile(diff,0.95)-1

TS_stan_state_only_plot_MARSS(sp='Goliath Grouper',GZ='Florida (All Regions)',params1=gg_params,TT=27,ts=gg_ts)
dev.off()

##Florida Keys to Dry Tortugas trajectory #####
gg_occs_fk<- gg_occs[which(substring(gg_occs$geogr4,1,2)==34),]
gg_ts_fk<- ts_reef(gg_occs_fk,sp=goliath)
length(unique(gg_occs_fk$site))

##FK model
X<- matrix(data=c(scale(as.numeric(gg_occs_fk$btime)),scale(as.numeric(gg_occs_fk$averagedepth)),scale(as.numeric(gg_occs_fk$visibility)),scale(as.numeric(gg_occs_fk$current)),gg_occs_fk$exp_binary),ncol=5,nrow=nrow(gg_occs_fk))

gg_SS_fk<- rstan::stan(model_code = SS_trend_ord, data = list(y =gg_occs_fk$abundance2,
                                                           N = nrow(gg_occs_fk),
                                                           site=as.numeric(factor(gg_occs_fk$geogr)),
                                                           N_site=length(unique(gg_occs_fk$geogr)),
                                                           diver=as.numeric(factor(gg_occs_fk$fish_memberid)),
                                                           N_dv=length(unique(gg_occs_fk$fish_memberid)),
                                                           dmy=as.numeric(factor(gg_occs_fk$site_dmy)),
                                                           N_dmy=length(unique(gg_occs_fk$site_dmy)),
                                                           K=length(unique(gg_occs_fk$abundance)),
                                                           X=X,
                                                           Z=ncol(X),
                                                           TT=27,
                                                           N_yr=length(unique(gg_occs_fk$year)),
                                                           yr_index=sort(unique(as.numeric(factor(gg_occs_fk$year)))),
                                                           year_id=as.numeric(factor(gg_occs_fk$year))),
                    pars = c('c','sd_site','sd_dv','sd_dmy','sd_r','sd_q','x','a_yr'),
                    control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)

gg_params_fk<- rstan::extract(gg_SS_fk)

TS_stan_state_only_plot_MARSS(sp='Goliath Grouper',GZ='Florida Keys',params1=gg_params_fk,TT=27,ts=gg_ts_fk)
dev.off()

diff=(ord_to_n(gg_params_fk$x[,27],gg_params_fk$c)/ord_to_n(gg_params_fk$x[,17],gg_params_fk$c))
median(diff)-1
quantile(diff,0.05)-1
quantile(diff,0.95)-1

##East Coast population trajectory ####
gg_occs_ec<- gg_occs[which(substring(gg_occs$geogr4,1,2)==33),]
gg_ts_ec<- ts_reef(gg_occs_ec,sp=goliath)
length(unique(gg_occs_ec$site))

X<- matrix(data=c(scale(as.numeric(gg_occs_ec$btime)),scale(as.numeric(gg_occs_ec$averagedepth)),scale(as.numeric(gg_occs_ec$visibility)),scale(as.numeric(gg_occs_ec$current)),gg_occs_ec$exp_binary),ncol=5,nrow=nrow(gg_occs_ec))

gg_SS_ec<- rstan::stan(model_code = SS_trend_ord, data = list(y =gg_occs_ec$abundance2,
                                                              N = nrow(gg_occs_ec),
                                                              site=as.numeric(factor(gg_occs_ec$geogr)),
                                                              N_site=length(unique(gg_occs_ec$geogr)),
                                                              diver=as.numeric(factor(gg_occs_ec$fish_memberid)),
                                                              N_dv=length(unique(gg_occs_ec$fish_memberid)),
                                                              dmy=as.numeric(factor(gg_occs_ec$site_dmy)),
                                                              N_dmy=length(unique(gg_occs_ec$site_dmy)),
                                                              K=length(unique(gg_occs_ec$abundance)),
                                                              X=X,
                                                              Z=ncol(X),
                                                              TT=27,
                                                              N_yr=length(unique(gg_occs_ec$year)),
                                                              yr_index=sort(unique(as.numeric(factor(gg_occs_ec$year)))),
                                                              year_id=as.numeric(factor(gg_occs_ec$year))),
                       pars = c('c','sd_site','sd_dv','sd_dmy','sd_r','sd_q','x','a_yr','beta'),
                       control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)

gg_params_ec<- rstan::extract(gg_SS_ec)

diff=(ord_to_n(gg_params_ec$x[,27],gg_params_ec$c)/ord_to_n(gg_params_ec$x[,17],gg_params_ec$c))
median(diff)-1
quantile(diff,0.05)-1
quantile(diff,0.95)-1

TS_stan_state_only_plot_MARSS(sp='Goliath Grouper',GZ='Florida (East Coast)',params1=gg_params_ec,TT=27,ts=gg_ts_ec)
dev.off()

##Gulf of Mexico - data summary###
##East Coast
gg_occs_mx<- gg_occs[which(substring(gg_occs$geogr4,1,1)==2),]
gg_occs_mx<- subset(gg_occs_mx,year>=1998)
gg_ts_mx<- ts_reef(gg_occs_ec,sp=goliath)
length(unique(gg_occs_mx$site))


####Binomial occurrence model - Florida ####
##Stan Model - bernoulli occurrence####
SS_trend_binom<-"
data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
}
parameters {
  real x0; //initial popn size

  //deviations from intercept
  vector[Z] beta; //effort coefficients - RVC
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  real<lower = 0> sd_r;
  real<lower = 0> sd_q;
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr] obs_dev; //observation deviations 
}

transformed parameters{
  vector[TT] x;
  vector[N_yr] a_yr;

  x[1] = x0 + pro_dev[1];
   
  for(t in 2:TT){
    x[t] = x[t-1] + u + pro_dev[t];
  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + obs_dev[i]; 
  }
}  

model{
  //priors
  beta ~ normal(0,2); //covariates
  x0 ~ normal(0,5); //initial state

  //variance terms
  sd_q ~inv_gamma(2,0.25);
  sd_r ~ inv_gamma(2,0.25);
  sd_site ~ inv_gamma(2, 1);
  sd_dv ~ inv_gamma(3, 0.5);
  sd_dmy ~ inv_gamma(3, 2);
  
  //varying intercepts
  a_site ~ normal(0, sd_site);
  a_dv ~ normal(0, sd_dv);
  a_dmy ~ normal(0, sd_dmy);
  
  for(t in 1:TT){
    pro_dev[t] ~ normal(0, sd_q);
    obs_dev[t] ~ normal(0,sd_r);
  }
  
  for(1 in r:N_regions){
    a_site ~ normal(beta_region,region_sd)
  }

  y ~ bernoulli_logit(a_yr[year_id]+a_site[site]+a_dv[diver]+a_dmy[dmy]+X*beta);
}
"
