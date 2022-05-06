'%notin%' <- Negate(`%in%`)

reef_filter_sp = function(R,GZ,sp,invert){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% subset(geogr4 %in% GZ) %>% select('formid','speciesid','abundance',everything())
  if(invert==0){
    TempDat<- subset(TempDat,type!=2) # take R and then subset any geogr4 in GZ, don't care about GZ, delete subset(geogr4 in GZ) bit, depending on what your needs are, you can modify this
  }
  if(invert==1){
    TempDat<- subset(TempDat,type!=1)
  }
  if(invert==2){
    TempDat<- subset(TempDat,type!=1) # remove this part if want to do for entire TWA
  }
  #TempDat_sp<- TempDat %>% subset(speciesid==sp)
  TempDat_sp<- TempDat %>% subset(speciesid%in%sp) # taking all data for sp of interest, if the sp matches this, then include it in our subset
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n=n(),n.year=length(unique(year))) %>% subset(n.year>2) # piping that sp data into group_by function, group all data based on geogr which is site, then summarise is saying, summarise by site 1st n = n() # of obs, and n.year= saying take length of unque yrs for each site, will figure out how many diff yrs rep in data based on surveys at site, so will calc total # of yrs that site has data for, final section, give back summarised data table but only if 3 or more years of data
  #summarise n=n will count # sightings for each site and year
  #tempdatsite only showing surveys where you saw the sp
  #to include 0s for total survey effort, separate graph for proportion of surveys where sighted...
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  Zeros<- tidyr::complete(TempDat2,formid,tidyr::nesting(speciesid),
                   fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat2$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat2[m,4:ncol(TempDat2)] #Replicate the survey-level data (except abundance)
  TempDat3<- rbind(TempDat2,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat3<- subset(TempDat3,speciesid==sp)
  if(nrow(TempDat3)==0){return(TempDat3)}
  TempDat3$occ=ifelse(TempDat3$abundance>0,1,0) #Code presence/absence based on abundance
  # turning into binary presence absence
  TempDat3<- TempDat3 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat3$abundance2<- TempDat3$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat3$abund_trans=NA #Transform abundance categories into the minimum counts
  TempDat3<- TempDat3%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
  surveyors<- TempDat3 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
  surveyors_trim<- subset(surveyors,n>=3) #only keep surveys by members with 10 or more dives
  
  TempDat4<- subset(TempDat3, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveyors with less than 5
  site_subsets<-TempDat4 %>% group_by(geogr) %>% summarize(n=n(),n.y=n_distinct(year)) %>% subset(n>=3) #Calculate surveys per site
  
  TempDat5<- TempDat4 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  survs<- TempDat5 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat5<- TempDat5 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  return(TempDat5)

}


reef_filter_sp_nogz = function(R,sp){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>%  select('formid','speciesid','abundance',everything())
  
  #TempDat_sp<- TempDat %>% subset(speciesid==sp)
  TempDat_sp<- TempDat %>% subset(speciesid%in%sp) # taking all data for sp of interest, if the sp matches this, then include it in our subset
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n=n(),n.year=length(unique(year))) %>% subset(n.year>2) # piping that sp data into group_by function, group all data based on geogr which is site, then summarise is saying, summarise by site 1st n = n() # of obs, and n.year= saying take length of unque yrs for each site, will figure out how many diff yrs rep in data based on surveys at site, so will calc total # of yrs that site has data for, final section, give back summarised data table but only if 3 or more years of data
  #summarise n=n will count # sightings for each site and year
  #tempdatsite only showing surveys where you saw the sp
  #to include 0s for total survey effort, separate graph for proportion of surveys where sighted...
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  Zeros<- tidyr::complete(TempDat2,formid,tidyr::nesting(speciesid),
                          fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat2$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat2[m,4:ncol(TempDat2)] #Replicate the survey-level data (except abundance)
  TempDat3<- rbind(TempDat2,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat3<- subset(TempDat3,speciesid %in% sp)
  if(nrow(TempDat3)==0){return(TempDat3)}
  TempDat3$occ=ifelse(TempDat3$abundance>0,1,0) #Code presence/absence based on abundance
  # turning into binary presence absence
  TempDat3<- TempDat3 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat3$abundance2<- TempDat3$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat3$abund_trans=NA #Transform abundance categories into the minimum counts
  TempDat3<- TempDat3%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
  surveyors<- TempDat3 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
  surveyors_trim<- subset(surveyors,n>=3) #only keep surveys by members with 10 or more dives
  
  TempDat4<- subset(TempDat3, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveyors with less than 5
  site_subsets<-TempDat4 %>% group_by(geogr) %>% summarize(n=n(),n.y=n_distinct(year)) %>% subset(n>=3) #Calculate surveys per site
  
  TempDat5<- TempDat4 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  survs<- TempDat5 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat5<- TempDat5 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  return(TempDat5)
  
}




#reef_filter_sp_3(R=R,GZ=3404,sp=c(661),invert=0)

reef_filter_sp_3 = function(R,GZ,sp,invert){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% subset(geogr4 %in% GZ) %>% select('formid','speciesid','abundance',everything())
  # take R and then subset any geogr4 in GZ, don't care about GZ, delete subset(geogr4 in GZ) bit, depending on what your needs are, you can modify this 
  if(invert==0){
    TempDat<- subset(TempDat,type!=2) 
  }
  if(invert==1){
    TempDat<- subset(TempDat,type!=1)
  }
  if(invert==2){
    TempDat<- subset(TempDat,type!=1) # remove this part if want to do for entire TWA
  }
  #TempDat_sp<- TempDat %>% subset(speciesid==sp)
  TempDat_sp<- TempDat %>% subset(speciesid%in%sp) # taking all data for sp of interest, if the sp matches this, then include it in our subset
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n=n(),n.year=length(unique(year))) %>% subset(n.year>2) # piping that sp data into group_by function, group all data based on geogr which is site, then summarise is saying, summarise by site 1st n = n() # of obs, and n.year= saying take length of unque yrs for each site, will figure out how many diff yrs rep in data based on surveys at site, so will calc total # of yrs that site has data for, final section, give back summarised data table but only if 3 or more years of data
  #summarise n=n will count # sightings for each site and year
  #tempdatsite only showing surveys where you saw the sp
  #to include 0s for total survey effort, separate graph for proportion of surveys where sighted...
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  Zeros<- tidyr::complete(TempDat2,formid,tidyr::nesting(speciesid),
                          fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat2$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat2[m,4:ncol(TempDat2)] #Replicate the survey-level data (except abundance)
  TempDat3<- rbind(TempDat2,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat3<- subset(TempDat3,speciesid==sp)
  if(nrow(TempDat3)==0){return(TempDat3)}
  TempDat3$occ=ifelse(TempDat3$abundance>0,1,0) #Code presence/absence based on abundance
  # turning into binary presence absence
  TempDat3<- TempDat3 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat3$abundance2<- TempDat3$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat3$abund_trans=NA #Transform abundance categories into the minimum counts
  TempDat3<- TempDat3%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
  surveyors<- TempDat3 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
  surveyors_trim<- subset(surveyors,n>=3) #only keep surveys by members with 10 or more dives
  
  TempDat4<- subset(TempDat3, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveyors with less than 5
  site_subsets<-TempDat4 %>% group_by(geogr) %>% summarize(n=n(),n.y=n_distinct(year)) %>% subset(n>=3) #Calculate surveys per site
  
  TempDat5<- TempDat4 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  survs<- TempDat5 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n(),propSurv = sum(TempDat3$occ)/n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat5<- TempDat5 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  return(TempDat5)
  
}




# add line after all filtering, then take tempdat5, groupby zone code, year, breaks it down by both, summarise n, or you can summarise, calc proportion sum of occ/n, knowing # of sites is helpful but breaking it down by sites is too fine scale
# total # of sites they consistently occur, otherwise by zone code
# could do survey effort at that scale as well, kind of depends on how many sites there are, if not that many sites then just fine to summarise based on GZ
# easiest way to navigate is to go line by line, but to do thta need to define R=R, GZ=GZ, sp=661, invert=0, then if do that, can start at top, at occ...and just go down and do line by line to SEE what they actually do then can figure out what these labels all mean
#run tempdat line and then CCALL IT and see how it differs
#what happens when go from tempdat2-3
# return(TempDat5) - if want fn to just tell you the # of tempdat_site, chg return to return (TempDatsite) can remove EVERYTHING else and just make it a new function
# very often have to repeat same thing over and over
# can just tell fn what want it to do and just plug in the approp varis
# so can modify varis at top and return to whatever you want
# getting # of surveys over time, could do based on zone codes
#going back to tempdat site...group by geogr, so nice trick: add comma





reef_filter_sp_2 = function(R,sp,invert){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% select('formid','speciesid','abundance',everything())
  if(invert==0){
    TempDat<- subset(TempDat,type!=2) # take R and then subset any geogr4 in GZ, don't care about GZ, delete subset(geogr4 in GZ) bit, depending on what your needs are, you can modify this
  }
  if(invert==1){
    TempDat<- subset(TempDat,type!=1)
  }
  if(invert==2){
    TempDat<- subset(TempDat,type!=1) # remove this part if want to do for entire TWA
  }
  #TempDat_sp<- TempDat %>% subset(speciesid==sp)
  TempDat_sp<- TempDat %>% subset(speciesid%in%sp) # taking all data for sp of interest, if the sp matches this, then include it in our subset
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n=n(),n.year=length(unique(year))) %>% subset(n.year>2) # piping that sp data into group_by function, group all data based on geogr which is site, then summarise is saying, summarise by site 1st n = n() # of obs, and n.year= saying take length of unque yrs for each site, will figure out how many diff yrs rep in data based on surveys at site, so will calc total # of yrs that site has data for, final section, give back summarised data table but only if 3 or more years of data
  #summarise n=n will count # sightings for each site and year
  #tempdatsite only showing surveys where you saw the sp
  #to include 0s for total survey effort, separate graph for proportion of surveys where sighted...
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  Zeros<- tidyr::complete(TempDat2,formid,tidyr::nesting(speciesid),
                          fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat2$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat2[m,4:ncol(TempDat2)] #Replicate the survey-level data (except abundance)
  TempDat3<- rbind(TempDat2,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat3<- subset(TempDat3,speciesid==sp)
  if(nrow(TempDat3)==0){return(TempDat3)}
  TempDat3$occ=ifelse(TempDat3$abundance>0,1,0) #Code presence/absence based on abundance
  # turning into binary presence absence
  TempDat3<- TempDat3 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat3$abundance2<- TempDat3$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat3$abund_trans=NA #Transform abundance categories into the minimum counts
  TempDat3<- TempDat3%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
  surveyors<- TempDat3 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
  surveyors_trim<- subset(surveyors,n>=3) #only keep surveys by members with 10 or more dives
  
  TempDat4<- subset(TempDat3, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveyors with less than 5
  site_subsets<-TempDat4 %>% group_by(geogr) %>% summarize(n=n(),n.y=n_distinct(year)) %>% subset(n>=3) #Calculate surveys per site
  
  TempDat5<- TempDat4 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  survs<- TempDat5 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat5<- TempDat5 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  return(TempDat5)
  
}











reef_filter_geog = function(R,geog,sp,invert){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% subset(geogr %in% geog) %>% select('formid','speciesid','abundance',everything())
  if(invert==0){
    TempDat<- subset(TempDat,type!=2)
  }
  if(invert==1){
    TempDat<- subset(TempDat,type!=1)
  }
  if(invert==2){
    TempDat<- subset(TempDat,type!=1)
  }
  TempDat_sp<- TempDat %>% subset(speciesid==sp)
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n=n(),n.year=length(unique(year))) %>% subset(n.year>2)
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  Zeros<- tidyr::complete(TempDat2,formid,nesting(speciesid),
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
 # surveyors_trim<- subset(surveyors,n>=3) #only keep surveys by members with 5 or more dives
  
#  TempDat4<- subset(TempDat3, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveys with less than 5
  site_subsets<-TempDat3 %>% group_by(geogr) %>% summarize(n=n_distinct(formid)) %>% subset(n>=3) #Calculate surveys per site
  
  TempDat5<- TempDat3 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  survs<- TempDat5 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat5<- TempDat5 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  
  occ_dat<- subset(TempDat5,speciesid==sp) #Subset out each species in the provided dataframe
  return(occ_dat)
}


State_space_timeseries_plot<- function(sp,GZ,params1,TT,ts,cols){
  #  pdf(paste('./figures/',paste(sp,GZ,'state_space',sep='_'),'.pdf',sep=''),width=8,height=6)
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params1$c)))
    
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
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
      reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,i])
      reef_coef[,10]<- 0
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==4){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
      reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,i])-plogis(params1$c[,3]-params1$a_yr[,i])
      reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,i])
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
      reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,i])-plogis(params1$c[,3]-params1$x[,i])
      reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,i])
      
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
  plot(y_mat$median.reef~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(na.omit(y_mat[,2])),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col=cols[1])
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  polygon(x, y1, col = adjustcolor(cols[1], alpha = 0.2), border=NA) # Add uncertainty polygon
  
  lines(y_mat$median.reef~y_mat$year,col=cols[2],lwd=2)
  points(y_mat$median.reef~y_mat$year,col='white',pch=21,bg=cols[2],cex=1.5)
  text(y=rep(min(x_mat),nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  #  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
  #  dev.off()
}

Multi_state_space_timeseries_plot<- function(sp,grp,params1,TT,ts,n.groups,col.palette1,col.palette2,log=0){
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params1$c)))
      
      if(ncol(params1$c)==2){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==3){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==4){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,,r][,i])-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,,r][,i])
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,,r][,i])-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,,r][,i])
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
    }
    lambda_mat[[r]]=ts_estimates
  } 
  
  
  for(r in 1:n.groups){  
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x)
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.95)
    }
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y)
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.95)
    }
  }
  x_mat_full<- do.call(rbind, lapply(x_mat, data.frame, stringsAsFactors=FALSE))
  y_mat_full<- do.call(rbind, lapply(y_mat, data.frame, stringsAsFactors=FALSE))
  par(xpd=T,mar=c(4.1,4.1,3.1,3.1))
  if(log==1){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=sp,log='y')
  }
  if(log==0){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=sp)
  }
  for(r in 1:n.groups){
    lines(x_mat[[r]][,1]~y_mat[[r]]$year,lty=5,lwd=2,col=col.palette1[r])
    x<- c(y_mat[[r]]$year, rev(y_mat[[r]]$year))
    y.90<- c(x_mat[[r]][,4], rev(x_mat[[r]][,5]))
#    y.80<- c(x_mat[[r]][,2], rev(x_mat[[r]][,3]))
    polygon(x, y.90, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
#    polygon(x, y.80, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    
    lines(y_mat[[r]]$median.reef~y_mat[[r]]$year,col=col.palette2[r],lwd=2)
    points(y_mat[[r]]$median.reef~y_mat[[r]]$year,col='white',pch=21,bg=col.palette2[r],cex=1.5)
  }
  text(y=rep(par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.015,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  text(y=par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.045,x=par("usr")[2]+0.2,'Surveys',cex=0.6,pos=1,xpd=T,font=1)
  for(r in 1:n.groups){
    lines(y=rep(par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,2),x=c(par("usr")[2]-0.5,par("usr")[2]+0.5),col=col.palette2[r],xpd=T,lwd=2)
  }
  for(r in 1:n.groups){
    text(y=par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,x=par("usr")[2]-2.5,grp[r],col=col.palette2[r],xpd=T,lwd=2,cex=0.6)
  }
}

Multi_state_space_timeseries_plot_scaled<- function(sp,grp,params1,TT,ts,n.groups,col.palette1,col.palette2,log){
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()
  median.x= numeric(3)
  median.y= numeric(3)
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params$c)))
      
      if(ncol(params1$c)==2){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==3){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==4){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,,r][,i])-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,,r][,i])
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,,r][,i])-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,,r][,i])
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
    }
    lambda_mat[[r]]=ts_estimates
    median.x[r]=median(do.call(rbind, lapply(lambda_mat[[r]], data.frame, stringsAsFactors=FALSE))$lambda.x)
    median.y[r]=median(do.call(rbind, lapply(lambda_mat[[r]], data.frame, stringsAsFactors=FALSE))$lambda.y)
  } 
  
  
  for(r in 1:n.groups){  
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x/median.x[r])
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.95)
    }
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y/median.y[r])
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.95)
    }
  }
  x_mat_full<- do.call(rbind, lapply(x_mat, data.frame, stringsAsFactors=FALSE))
  y_mat_full<- do.call(rbind, lapply(y_mat, data.frame, stringsAsFactors=FALSE))
  par(xpd=T,mar=c(4.1,4.1,3.1,3.1))
  if(log==1){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=sp,log='y')
  }
  if(log==0){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=sp)
  }
  for(r in 1:n.groups){
    lines(x_mat[[r]][,1]~y_mat[[r]]$year,lty=5,lwd=2,col=col.palette1[r])
    x<- c(y_mat[[r]]$year, rev(y_mat[[r]]$year))
    y.90<- c(x_mat[[r]][,4], rev(x_mat[[r]][,5]))
    #    y.80<- c(x_mat[[r]][,2], rev(x_mat[[r]][,3]))
    polygon(x, y.90, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    #    polygon(x, y.80, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    
    lines(y_mat[[r]]$median.reef~y_mat[[r]]$year,col=col.palette2[r],lwd=2)
    points(y_mat[[r]]$median.reef~y_mat[[r]]$year,col='white',pch=21,bg=col.palette2[r],cex=1.5)
  }
  text(y=rep(par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.015,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  text(y=par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.045,x=par("usr")[2]+0.2,'Surveys',cex=0.6,pos=1,xpd=T,font=1)
  for(r in 1:n.groups){
    lines(y=rep(par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,2),x=c(par("usr")[2]-0.5,par("usr")[2]+0.5),col=col.palette2[r],xpd=T,lwd=2)
  }
  for(r in 1:n.groups){
    text(y=par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,x=par("usr")[2]-2.5,grp[r],col=col.palette2[r],xpd=T,lwd=2,cex=0.6)
  }

}

Multi_state_space_timeseries_plot_pdf<- function(sp,grp,params1,TT,ts,n.groups,col.palette1,col.palette2){
  pdf(file.path(path,paste(sprintf("%02d",i),'_',sp,'_',GZ,'.pdf',sep='')),width=8,height=6)
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params$c)))
      
      if(ncol(params1$c)==2){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==3){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==4){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,,r][,i])-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,,r][,i])
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,,r][,i])-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,,r][,i])
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
    }
    lambda_mat[[r]]=ts_estimates
  } 
  
  
  for(r in 1:n.groups){  
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x)
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.95)
    }
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y)
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.95)
    }
  }
  x_mat_full<- do.call(rbind, lapply(x_mat, data.frame, stringsAsFactors=FALSE))
  y_mat_full<- do.call(rbind, lapply(y_mat, data.frame, stringsAsFactors=FALSE))
  par(xpd=T,mar=c(4.1,4.1,3.1,3.1))
  plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=sp)
  for(r in 1:n.groups){
    lines(x_mat[[r]][,1]~y_mat[[r]]$year,lty=5,lwd=2,col=col.palette1[r])
    x<- c(y_mat[[r]]$year, rev(y_mat[[r]]$year))
    y.90<- c(x_mat[[r]][,4], rev(x_mat[[r]][,5]))
    #    y.80<- c(x_mat[[r]][,2], rev(x_mat[[r]][,3]))
    polygon(x, y.90, col = adjustcolor(col.palette1[r], alpha = 0.1), border=NA) # Add uncertainty polygon
    #    polygon(x, y.80, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    
    lines(y_mat[[r]]$median.reef~y_mat[[r]]$year,col=col.palette2[r],lwd=2)
    points(y_mat[[r]]$median.reef~y_mat[[r]]$year,col='white',pch=21,bg=col.palette2[r],cex=1.5)
  }
  text(y=rep(par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.015,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  text(y=par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.045,x=par("usr")[2]+0.2,'Surveys',cex=0.6,pos=1,xpd=T,font=1)
  for(r in 1:n.groups){
    lines(y=rep(par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,2),x=c(par("usr")[2]-0.5,par("usr")[2]+0.5),col=col.palette2[r],xpd=T,lwd=2)
  }
  for(r in 1:n.groups){
    text(y=par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,x=par("usr")[2]-2.5,grp[r],col=col.palette2[r],xpd=T,lwd=2,cex=0.6)
  }
  
  dev.off(file.path(path,paste(sprintf("%02d",i),'_',sp,'_',GZ,'.pdf',sep='')))
  dev.off()
}

State_space_timeseries_plot_pdf<- function(sp,GZ,params1,TT,ts,path,i){
  pdf(file.path(path,paste(sprintf("%02d",i),'_',sp,'_',GZ,'.pdf',sep='')),width=8,height=6)
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params1$c)))
    
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
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
      reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,i])
      reef_coef[,10]<- 0
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==4){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
      reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,i])-plogis(params1$c[,3]-params1$a_yr[,i])
      reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,i])
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
      reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,i])-plogis(params1$c[,3]-params1$x[,i])
      reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,i])
      
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
  
  lines(y_mat$median.reef~y_mat$year,col='dodgerblue4',lwd=2)
  points(y_mat$median.reef~y_mat$year,col='white',pch=21,bg='dodgerblue4',cex=1.5)
  text(y=rep(min(x_mat),nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  dev.off(file.path(path,paste(sprintf("%02d",i),'_',sp,'_',GZ,'.pdf',sep='')))
  dev.off()
}


ts_reef = function(X){
  occ_by_year<- X %>% group_by(year) %>% summarize(n.occ=sum(occ),n.surv=n_distinct(formid),p.occ=n.occ/n.surv)
  abun_by_year<- X %>% group_by(year,geogr) %>% summarize(site_abun=mean(abund_trans),sd_abund=sd(abund_trans),n.surv=n_distinct(formid)) %>% group_by(year) %>% summarize(mean_abund=mean(site_abun),n.survs=sum(n.surv),n.sites=n())
  total_sd<- X %>% group_by(year) %>% summarize(sd=sd(abund_trans))
  
  comb<- left_join(occ_by_year,abun_by_year)
  comb2<- left_join(comb,total_sd)
  ts_dat<- comb2
  return(ts_dat)
}

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

population_growth_rates<- function(params1,TT){
lambda_mat<- list()  
for(i in 1:TT){
  reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params$c)))
  
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
    reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
    reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,i])
    reef_coef[,10]<- 0
    
    reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
    reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
  }
  
  if(ncol(params1$c)==4){
    reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
    reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
    reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
    reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,i])-plogis(params1$c[,3]-params1$a_yr[,i])
    reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,i])
    reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
    reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,i])-plogis(params1$c[,3]-params1$x[,i])
    reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,i])
    
    reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
    reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
  }
  lambda_mat[[i]]=reef_coef
}  


x_mat<- matrix(nrow=nrow(lambda_mat[[1]]),ncol=TT)
for(i in 1:TT){
  x_mat[,i]=lambda_mat[[i]]$lambda.x
}
r_mat<- matrix(nrow=nrow(lambda_mat[[1]]),ncol=TT-1)
for(i in 1:nrow(r_mat)){
  for(t in 1:TT-1)
    r_mat[i,t]=log(x_mat[i,t+1])-log(x_mat[i,t])
}
return(r_mat)
}

gm_mean_r<- function(x){
  prod(exp(x))^(1/length(x))
}

population_growth_rates_multi<- function(params1,TT,n.groups){
lambda_mat<- list()
x_mat<- list()
r_mat<- list()
for(r in 1:n.groups){
  ts_estimates<- list()
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params$c)))
    
    if(ncol(params1$c)==2){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
      reef_coef[,4]<- 0
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,,r][,i])
      reef_coef[,9]<- 0
      reef_coef[,10]<- 0
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==3){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
      reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
      reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,,r][,i])
      reef_coef[,10]<- 0
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==4){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
      reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,,r][,i])-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
      reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,,r][,i])
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
      reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,,r][,i])-plogis(params1$c[,3]-params1$x[,,r][,i])
      reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,,r][,i])
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    ts_estimates[[i]]=reef_coef
  }
  lambda_mat[[r]]=ts_estimates
} 


for(r in 1:n.groups){  
  x_mat[[r]]<- matrix(nrow=nrow(lambda_mat[[1]][[1]]),ncol=TT)
  for(i in 1:TT){
    x_mat[[r]][,i]=lambda_mat[[r]][[i]]$lambda.x
  }
  r_mat[[r]]<- matrix(nrow=nrow(lambda_mat[[1]][[1]]),ncol=TT-1)
  for(i in 1:nrow(r_mat[[r]])){
    for(t in 1:TT-1)
      r_mat[[r]][i,t]=log(x_mat[[r]][i,t+1])-log(x_mat[[r]][i,t])
  }
}
return(r_mat) 
}
