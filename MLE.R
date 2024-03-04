library(maxLik)
library(numDeriv)


dat<-read.csv("van_data.csv")
lifetime<-c(0,0,0)
lifetime_L<-c(0,0,0)
lifetime_U<-c(0,0,0)
i=1


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Copper")
logLikFun<-function(par){
  decay_rate<-par[1];intercept<-par[2];sigma<-par[3];
  f=0;
  for(i in 1:nrow(data))
  {
    if(as.numeric(data[i,10])>as.numeric(data[i,9]))
      f<-f+log(pnorm((as.numeric(data[i,10])-intercept+decay_rate*as.numeric(data[i,6]))/sigma)
               -pnorm((as.numeric(data[i,10])-0.25-intercept+decay_rate*as.numeric(data[i,6]))/sigma))
    else
      f<-f+log(pnorm((as.numeric(data[i,9])-intercept+decay_rate*as.numeric(data[i,6]))/sigma))
  }
  return(f)
}
MLE<-maxLik(logLik=logLikFun,start=c(0.1,5,0.5))
par<-coef(MLE)
lifetime[i]<-3/par[1]
hessian<--1*hessian(logLikFun,par)
I=solve(hessian)
a<-qnorm(0.975)*sqrt((3/(par[1])^2)^2*as.numeric(I[1,1]))
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a
i=i+1


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Plastic")
logLikFun<-function(par){
  decay_rate<-par[1];intercept<-par[2];sigma<-par[3];
  f=0;
  for(i in 1:nrow(data))
  {
    if(as.numeric(data[i,10])>as.numeric(data[i,9]))
      f<-f+log(pnorm((as.numeric(data[i,10])-intercept+decay_rate*as.numeric(data[i,6]))/sigma)
               -pnorm((as.numeric(data[i,10])-0.25-intercept+decay_rate*as.numeric(data[i,6]))/sigma))
    else
      f<-f+log(pnorm((as.numeric(data[i,9])-intercept+decay_rate*as.numeric(data[i,6]))/sigma))
  }
  return(f)
}
MLE<-maxLik(logLik=logLikFun,start=c(0.1,5,0.5))
par<-coef(MLE)
lifetime[i]<-3/par[1]
hessian<--1*hessian(logLikFun,par)
I=solve(hessian)
a<-qnorm(0.975)*sqrt((3/(par[1])^2)^2*(as.numeric(I[1,1])))
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a
i=i+1


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Steel")
logLikFun<-function(par){
  decay_rate<-par[1];intercept<-par[2];sigma<-par[3];
  f=0;
  for(i in 1:nrow(data))
  {
    if(as.numeric(data[i,10])>as.numeric(data[i,9]))
      f<-f+log(pnorm((as.numeric(data[i,10])-intercept+decay_rate*as.numeric(data[i,6]))/sigma)
               -pnorm((as.numeric(data[i,10])-0.25-intercept+decay_rate*as.numeric(data[i,6]))/sigma))
    else
      f<-f+log(pnorm((as.numeric(data[i,9])-intercept+decay_rate*as.numeric(data[i,6]))/sigma))
  }
  return(f)
}
MLE<-maxLik(logLik=logLikFun,start=c(0.1,5,0.5))
par<-coef(MLE)
lifetime[i]<-3/par[1]
hessian<--1*hessian(logLikFun,par)
I=solve(hessian)
a<-qnorm(0.975)*sqrt((3/(par[1])^2)^2*(as.numeric(I[1,1])))
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a


lifetime_table<-data.frame("0.025"=lifetime_L,"lifetime"=lifetime,"0.975"=lifetime_U)
write.csv(lifetime_table,file="lifetime_table_MLE.csv")

