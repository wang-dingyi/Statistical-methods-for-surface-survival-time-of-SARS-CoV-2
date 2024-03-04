library(icenReg)


dat<-read.csv("van_data.csv")
lifetime<-c(0,0,0)
lifetime_L<-c(0,0,0)
lifetime_U<-c(0,0,0)
i=1
sol<-function(f,a,b,e){
    while(abs(a-b)>e){
      c<-(a+b)/2
      if(f(c)*f(a)>0){a=c}
      else if(f(c)*f(a)==0){return(c)}
           else{b=c}
    }
    return(c)
  }


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Copper")
dll<-data$detection_limit_log10_titer[1]
data<-as.data.frame(cbind(data$time,data$log10_titer))
for(k in 1:11){
  if(mean(data[k,2]+data[k+11,2]+data[k+22,2]<3*dll+0.25)){
    data[k,2]<-0
    data[k+11,2]<-0
    data[k+22,2]<-0
  }
}
colnames(data)<-c("time","y2")
data<-data[which(data$y2>0), ]
y1<-data[,2]-0.25
for(k in 1:nrow(data)){
  if(data[k,2]<dll+0.25) {y1[k]<--50}
}
data<-as.data.frame(cbind(data[,1],y1,data[,2]))
colnames(data)<-c("time","y1","y2")
x<-data$time
y2<-data$y2
Eq<-function(par){
  k<-par
  data_fit<-as.data.frame(cbind(y1-k*x,y2-k*x,1)) 
  np_fit<-ic_np(cbind(y1-k*x,y2-k*x)~V3,data=data_fit)
  error<-np_fit$scurves$`1`$Tbull_ints
  error<-as.numeric(error[,1])
  Fun<-1-np_fit$scurves$`1`$S_curves$baseline
  error[1]<-0
  f1<-function(x){
    f<-0
    for(k in 2:length(error)){
      if(x>error[k]) f<-Fun[k]
    }
    return(f)
  }
  f2<-function(x){
    f<-0
    for(k in 2:length(error)){
      if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*error[k]
    }
    return(f)
  }
  
  xbar<-mean(x)
  h<-0
  for(a in 1:length(x)){
    E1<-f1(y2[a]-k*x[a])-f1(y1[a]-k*x[a])
    E2<-f2(y2[a]-k*x[a])-f2(y1[a]-k*x[a])
    h<-h+(x[a]-xbar)*E2/E1
  }
  return(h)
}
par2<-sol(Eq,-1,0.05,1e-4)
lifetime[i]<--3/par2
k<-par2
data_fit<-as.data.frame(cbind(y1-k*x,y2-k*x,1)) 
np_fit<-ic_np(cbind(y1-k*x,y2-k*x)~V3,data=data_fit)
error<-np_fit$scurves$`1`$Tbull_ints
error<-as.numeric(error[,1])
Fun<-1-np_fit$scurves$`1`$S_curves$baseline
error[1]<-0
f1<-function(x){
  f<-0
  for(k in 2:length(error)){
    if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*error[k]
  }
  return(f)
}

f2<-function(x){
  f<-0
  for(k in 2:length(error)){
    if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*(error[k])^2
  }
  return(f)
}
var<-f2(20)-(f1(20))^2
hessian<--1*(Eq(par2+0.001)-Eq(par2-0.001))/0.002/var
I=1/hessian
a<-qnorm(0.975)*(3/(par2)^2)*sqrt(I)
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a
i=i+1


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Plastic")
dll<-data$detection_limit_log10_titer[1]
data<-as.data.frame(cbind(data$time,data$log10_titer))
for(k in 1:11){
  if(mean(data[k,2]+data[k+11,2]+data[k+22,2]<3*dll+0.25)){
    data[k,2]<-0
    data[k+11,2]<-0
    data[k+22,2]<-0
  }
}
colnames(data)<-c("time","y2")
data<-data[which(data$y2>0), ]
y1<-data[,2]-0.25
for(k in 1:nrow(data)){
  if(data[k,2]<dll+0.25) {y1[k]<--50}
}
data<-as.data.frame(cbind(data[,1],y1,data[,2]))
colnames(data)<-c("time","y1","y2")
x<-data$time
y2<-data$y2
Eq<-function(par){
  k<-par
  data_fit<-as.data.frame(cbind(y1-k*x,y2-k*x,1)) 
  np_fit<-ic_np(cbind(y1-k*x,y2-k*x)~V3,data=data_fit)
  error<-np_fit$scurves$`1`$Tbull_ints
  error<-as.numeric(error[,1])
  Fun<-1-np_fit$scurves$`1`$S_curves$baseline
  error[1]<-0
  f1<-function(x){
    f<-0
    for(k in 2:length(error)){
      if(x>error[k]) f<-Fun[k]
    }
    return(f)
  }
  f2<-function(x){
    f<-0
    for(k in 2:length(error)){
      if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*error[k]
    }
    return(f)
  }
  
  xbar<-mean(x)
  h<-0
  for(a in 1:length(x)){
    E1<-f1(y2[a]-k*x[a])-f1(y1[a]-k*x[a])
    E2<-f2(y2[a]-k*x[a])-f2(y1[a]-k*x[a])
    h<-h+(x[a]-xbar)*E2/E1
  }
  return(h)
}
par2<-sol(Eq,-1,0.05,1e-4)
lifetime[i]<--3/par2
k<-par2
data_fit<-as.data.frame(cbind(y1-k*x,y2-k*x,1)) 
np_fit<-ic_np(cbind(y1-k*x,y2-k*x)~V3,data=data_fit)
error<-np_fit$scurves$`1`$Tbull_ints
error<-as.numeric(error[,1])
Fun<-1-np_fit$scurves$`1`$S_curves$baseline
error[1]<-0
f1<-function(x){
  f<-0
  for(k in 2:length(error)){
    if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*error[k]
  }
  return(f)
}

f2<-function(x){
  f<-0
  for(k in 2:length(error)){
    if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*(error[k])^2
  }
  return(f)
}
var<-f2(20)-(f1(20))^2
hessian<--1*(Eq(par2+0.001)-Eq(par2-0.001))/0.002/var
I=1/hessian
a<-qnorm(0.975)*(3/(par2)^2)*sqrt(I)
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a
i=i+1


data<-subset(dat,virus=="HCoV-19")
data<-subset(data,material=="Steel")
dll<-data$detection_limit_log10_titer[1]
data<-as.data.frame(cbind(data$time,data$log10_titer))
for(k in 1:11){
  if(mean(data[k,2]+data[k+11,2]+data[k+22,2]<3*dll+0.25)){
    data[k,2]<-0
    data[k+11,2]<-0
    data[k+22,2]<-0
  }
}
colnames(data)<-c("time","y2")
data<-data[which(data$y2>0), ]
y1<-data[,2]-0.25
for(k in 1:nrow(data)){
  if(data[k,2]<dll+0.25) {y1[k]<--50}
}
data<-as.data.frame(cbind(data[,1],y1,data[,2]))
colnames(data)<-c("time","y1","y2")
x<-data$time
y2<-data$y2
Eq<-function(par){
  k<-par
  data_fit<-as.data.frame(cbind(y1-k*x,y2-k*x,1)) 
  np_fit<-ic_np(cbind(y1-k*x,y2-k*x)~V3,data=data_fit)
  error<-np_fit$scurves$`1`$Tbull_ints
  error<-as.numeric(error[,1])
  Fun<-1-np_fit$scurves$`1`$S_curves$baseline
  error[1]<-0
  f1<-function(x){
    f<-0
    for(k in 2:length(error)){
      if(x>error[k]) f<-Fun[k]
    }
    return(f)
  }
  f2<-function(x){
    f<-0
    for(k in 2:length(error)){
      if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*error[k]
    }
    return(f)
  }
  
  xbar<-mean(x)
  h<-0
  for(a in 1:length(x)){
    E1<-f1(y2[a]-k*x[a])-f1(y1[a]-k*x[a])
    E2<-f2(y2[a]-k*x[a])-f2(y1[a]-k*x[a])
    h<-h+(x[a]-xbar)*E2/E1
  }
  return(h)
}
par2<-sol(Eq,-1,0.05,1e-4)
lifetime[i]<--3/par2
k<-par2
data_fit<-as.data.frame(cbind(y1-k*x,y2-k*x,1)) 
np_fit<-ic_np(cbind(y1-k*x,y2-k*x)~V3,data=data_fit)
error<-np_fit$scurves$`1`$Tbull_ints
error<-as.numeric(error[,1])
Fun<-1-np_fit$scurves$`1`$S_curves$baseline
error[1]<-0
f1<-function(x){
  f<-0
  for(k in 2:length(error)){
    if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*error[k]
  }
  return(f)
}

f2<-function(x){
  f<-0
  for(k in 2:length(error)){
    if(x>error[k]) f<-f+(Fun[k]-Fun[k-1])*(error[k])^2
  }
  return(f)
}
var<-f2(20)-(f1(20))^2
hessian<--1*(Eq(par2+0.001)-Eq(par2-0.001))/0.002/var
I=1/hessian
a<-qnorm(0.975)*(3/(par2)^2)*sqrt(I)
lifetime_L[i]<-lifetime[i]-a
lifetime_U[i]<-lifetime[i]+a


lifetime_table<-data.frame("0.025"=lifetime_L,"lifetime"=lifetime,"0.975"=lifetime_U)
write.csv(lifetime_table,file="lifetime_table_BJ.csv")

