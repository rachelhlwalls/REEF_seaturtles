rm(list=ls())
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Load in REEF functions####
source(here('code','reef_functions.R'))

#Stan model####
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
  int<lower=0> N_region; //number of regions
  int<lower=1,upper=N_region> region[N]; // vector of site identities
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
  vector[N_region] a_region; //deviation between basins
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_site;
  real<lower = 0> sd_reg;
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

  x[1] = x0 + pro_dev[1]*sd_q;
   
  for(t in 2:TT){
    x[t] = x[t-1] + pro_dev[t]*sd_q;
  }
   
  for(i in 1:N_yr){
    a_yr[i] = x[yr_index[i]] + obs_dev[i]*sd_r; 
  }
  
}  
model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  x0 ~ normal(0,5); //initial state

  //variance terms
  sd_q ~ inv_gamma(3,0.5);
  sd_r ~ inv_gamma(3,0.5);
  sd_site ~ inv_gamma(4, 2);
  sd_reg ~ inv_gamma(4, 2);
  sd_dv ~ inv_gamma(4, 2);
  sd_dmy ~ inv_gamma(4, 2);
  
  //varying intercepts
  for(s in 1:N_site){
   a_site[s] ~ normal(a_region[region[s]]*sd_reg,1);
  }
  a_region ~ std_normal();
  a_dv ~ std_normal();
  a_dmy ~ std_normal();
  
  for(t in 1:TT){
    pro_dev[t] ~ std_normal();
    obs_dev[t] ~ std_normal();
  }

  y ~ ordered_logistic(a_yr[year_id]+a_site[site]*sd_site+a_dv[diver]*sd_dv+a_dmy[dmy]*sd_dmy+X*beta,c);
  
}
"

#Batched models in data####
salish_sea_data<- list.files(here('data','species data','Salish Sea'))

for(i in 1:length(salish_sea_data)){
  sp_ss<- read.csv(here('data','species data','Salish Sea',salish_sea_data[i]),na.strings = "NULL")
  sp_ss<- subset(sp_ss,year<2021)
  
  X<- matrix(data=c(scale(as.numeric(sp_ss$btime)),scale(as.numeric(sp_ss$visibility)),scale(as.numeric(sp_ss$current)),sp_ss$exp_binary),ncol=4,nrow=nrow(sp_ss))
  
  stan_ss<- rstan::stan(model_code = SS_trend_ord, data = list(y=sp_ss$abundance2,
                                                              N = nrow(sp_ss),
                                                              site=as.numeric(factor(sp_ss$geogr)),
                                                              N_site=length(unique(sp_ss$geogr)),
                                                              region=as.numeric(factor(sp_ss$basin)),
                                                              N_region=length(unique(sp_ss$basin)),
                                                              diver=as.numeric(factor(sp_ss$fish_memberid)),
                                                              N_dv=length(unique(sp_ss$fish_memberid)),
                                                              dmy=as.numeric(factor(sp_ss$site_dmy)),
                                                              N_dmy=length(unique(sp_ss$site_dmy)),
                                                              K=length(unique(sp_ss$abundance2)),
                                                              TT=max(sp_ss$year)-min(sp_ss$year)+1,
                                                              N_yr=length(unique(sp_ss$year)),
                                                              yr_index=sort(unique(as.numeric(factor(sp_ss$year)))),
                                                              year_id=as.numeric(factor(sp_ss$year)),
                                                              X=X,
                                                              Z=ncol(X)),                                                         
                      pars = c('c','sd_site','sd_dv','sd_reg','sd_dmy','sd_r','sd_q','x','a_yr','a_region','a_site','beta'),
                      control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 600, thin = 1)
  
  params<- rstan::extract(stan_ss)
  
  ts_ss<- ts_reef(sp_ss)
  
  plot_path<- here('outputs','figures','Salish Sea','time series')
  State_space_timeseries_plot(sp=gsub('_ss.csv','',salish_sea_data[i]),GZ='Salish Sea',params1=params,TT=max(sp_ss$year)-min(sp_ss$year)+1,ts=ts_ss)
  State_space_timeseries_plot_pdf(sp=gsub('_ss.csv','',salish_sea_data[i]),GZ='Salish Sea',params1=params,TT=max(sp_ss$year)-min(sp_ss$year)+1,ts=ts_ss,path=plot_path,i)
  
  mod_par_path<- here('outputs','species parameter estimates','Salish Sea')
  write.csv(as.data.frame(params),file.path(mod_par_path,paste(sprintf("%02d",i),'_',gsub('_ss.csv','',salish_sea_data[i]),'.csv',sep='')))
  
  mod_par_sum_path<- here('outputs','species parameter summary','Salish Sea')
  pars_mod<- c('c','sd_site','sd_reg','sd_dv','sd_dmy','sd_r','sd_q','x','a_yr','beta','a_region','a_site')
  write.csv(as.data.frame(summary(stan_ss,pars = pars_mod)$summary),file.path(mod_par_sum_path,paste(sprintf("%02d",i),'_',gsub('_ss.csv','',salish_sea_data[i]),'.csv',sep='')))

}

