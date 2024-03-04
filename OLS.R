dat<-read.csv("van_data.csv")
lifetime<-c(0,0,0)
lifetime_L<-c(0,0,0)
lifetime_U<-c(0,0,0)
i=1


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Copper")
data<-as.data.frame(cbind(data[,6],data[,10]))
for(k in 1:11){
  if(mean(data[k,2]+data[k+11,2]+data[k+22,2]<4.6)){
    data[k,2]<-0
    data[k+11,2]<-0
    data[k+22,2]<-0
  }
}
colnames(data)<-c("time","log10_titer")
olsdata<-data[which(data$log10_titer>0), ]
ols=lm(olsdata$log10_titer~olsdata$time)
par<-coef(ols)
par2<-confint(ols)
lifetime[i]<--3/par[2]
a<-(3/(par[2])^2)*(par2[2,2]-par2[2,1])/2
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a
i=i+1


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Plastic")
data<-as.data.frame(cbind(data[,6],data[,10]))
for(k in 1:11){
  if(mean(data[k,2]+data[k+11,2]+data[k+22,2]<1.6)){
    data[k,2]<-0
    data[k+11,2]<-0
    data[k+22,2]<-0
  }
}
colnames(data)<-c("time","log10_titer")
olsdata<-data[which(data$log10_titer>0), ]
ols=lm(olsdata$log10_titer~olsdata$time)
par<-coef(ols)
par2<-confint(ols)
lifetime[i]<--3/par[2]
a<-(3/(par[2])^2)*(par2[2,2]-par2[2,1])/2
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a
i=i+1


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Steel")
data<-as.data.frame(cbind(data[,6],data[,10]))
for(k in 1:11){
  if(mean(data[k,2]+data[k+11,2]+data[k+22,2]<1.6)){
    data[k,2]<-0
    data[k+11,2]<-0
    data[k+22,2]<-0
  }
}
colnames(data)<-c("time","log10_titer")
olsdata<-data[which(data$log10_titer>0), ]
ols=lm(olsdata$log10_titer~olsdata$time)
par<-coef(ols)
par2<-confint(ols)
lifetime[i]<--3/par[2]
a<-(3/(par[2])^2)*(par2[2,2]-par2[2,1])/2
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a


lifetime_table<-data.frame("0.025"=lifetime_L,"lifetime"=lifetime,"0.975"=lifetime_U)
write.csv(lifetime_table,file="lifetime_table_OLS.csv")
