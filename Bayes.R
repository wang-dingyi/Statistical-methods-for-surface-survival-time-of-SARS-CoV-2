library(rstan)
library(magrittr)
library(tidybayes)


dat<-read.csv("van_data.csv")
data<-list(
  log_y=dat$log10_titer,
  dll=dat$detection_limit_log10_titer,
  material=dat$trial_unique_id,
  t=dat$time,
  lambda_mu=0.5,
  lambda_sigma=4,
  b_mu=4.5, 
  b_sigma=2.5,
  sigma_mu=0,
  sigma_sigma=2
  )
iter<-1500
chains<-8
init<-list()
for(i in 1:8){
  init1<-list(decay_rate=runif(3, 0, 0.5),sigma=runif(3, 20, 30),intercept=runif(3, 0, 0.2))
  init[i]<-list(init1)
} 


fit<-stan(
  file="Bayes_model.stan",
  data=data,
  iter=iter,
  chains=chains,
  init=init
  )
lambda<-fit%>%spread_draws(lambda[trial_unique_id])
lifetime<-c(0,0,0)
lifetime_L<-c(0,0,0)
lifetime_U<-c(0,0,0)
for(i in 1:3){
  lifetime[i]<-mean(as.numeric(lambda[lambda$trial_unique_id==i,]$lambda))
  lifetime_U[i]<-quantile(as.numeric(lambda[lambda$trial_unique_id==i,]$lambda),c(0.025))
  lifetime_L[i]<-quantile(as.numeric(lambda[lambda$trial_unique_id==i,]$lambda),c(0.975))
}


lifetime_table<-data.frame("0.025"=3/lifetime_L,"lifetime"=3/lifetime,"0.975"=3/lifetime_U)
write.csv(lifetime_table,file="lifetime_table_Bayes.csv")



