rm(list=ls())
library(survival)
library(dMod)
library(mvtnorm)
set.seed(1)
#setwd("C:/Users/lingyee.sze/OneDrive - University of Calgary/STAT631")
q_Prior_sample=read.table("prior_sample.txt",header=TRUE)

#Generate time to event data for n=300 subjects
n=300
z1<-c(rep(0, n/2), rep(1, n/2)) #treatment indicator

bio<-runif(n) #normalized biomarker

surv_time_generator<-function(beta1,beta2,beta3,c){
  beta<-c(beta1,beta2,beta3)
  I<-ifelse(bio>c,1,0)
  z = cbind(z1, I, z1*I) #trt I interact
  h0 = 1
  h = h0*exp(z%*%beta)
  t = rexp(n, h)   
  Cens<-runif(n,2,5)
  del<-ifelse(t<Cens,1,0)
  x<-  ifelse(t>Cens, Cens, t)
  data<-data.frame(cbind(x,del,z1,bio))
  names(data)<-c("obs_time","status","treatment","biomarker")
  return(data)
}

#scenario: exp(beta2)=1.5, exp(beta3)=1 , c=0.2
data1<-surv_time_generator(beta1=0,beta2 = log(1.5),beta3 = 0,c=0.2)
x<-data1$obs_time ; del<-data1$status 
mean_data1<-mean(data1$obs_time)

#joint posterior density of q,c,beta given D
lf1<-function(c,beta,q){ #use log density 
  z2<-ifelse(bio>c,1,0) ; z3<-z2*z1
  D<-data.frame(cbind(x,del,z1,z2,z3))
  beta.model<-coxph(Surv(x,del) ~ z1+z2+z3,data=D, init = beta)
  logdens<-log(c)+(q-1)*log(1-c)+log(q-1)+beta.model$loglik[1]
  return(logdens)
}


B=1000
S=3000
epsilon=0.01
c.v<-NULL  
beta.m<-matrix(0,S,3)
q.v<-NULL

q<-sample(q_Prior_sample$q,size=1)
c<-rbeta(1,2,q) ; c<-ifelse(c<0.05,0.05,ifelse(c>0.95,0.95,c))
z2<-ifelse(bio>c,1,0) ; z3<-z1*z2
beta.model<-coxph(Surv(x,del) ~ z1+z2+z3)
beta<-t(rmvnorm(1,beta.model$coefficients,beta.model$var))
D<-surv_time_generator(beta[1],beta[2],beta[3],c)
mean_D<-mean(D$obs_time)
dis<-abs(mean_D- mean_data1)
while (dis>epsilon){
  q<-sample(q_Prior_sample$q,size=1)
  c<-rbeta(1,2,q) ; c<-ifelse(c<0.05,0.05,ifelse(c>0.95,0.95,c))
  z2<-ifelse(bio>c,1,0) ; z3<-z1*z2
  beta.model<-coxph(Surv(x,del) ~ z1+z2+z3)
  beta<-t(rmvnorm(1,beta.model$coefficients,beta.model$var))
  D<-surv_time_generator(beta[1],beta[2],beta[3],c)
  mean_D<-mean(D$obs_time)
  dis<-abs(mean_D- mean_data1)
}
c.v[1]<-c  
beta.m[1,]<-t(beta)
q.v[1]<-q


for ( s in 2:S){
  
  #Propose c:
  lb<-c-0.3  ; ub<-c+0.3
  lb<-max(0.05,lb) ; ub<-min(0.95,ub)
  uk = runif(1,lb,ub)
  
  #Propose beta|c:
  z2<-ifelse(bio>uk,1,0) ; z3<-z1*z2
  D<-data.frame(cbind(x,del,z1,z2,z3))
  beta.model<-coxph(Surv(x,del) ~ z1+z2+z3,data=D)
  beta_prop<-t(rmvnorm(1,beta,beta.model$var))
  
  #Compute q|c,beta:
  lambda<--log(1-uk)
  v<-rgamma(1,2,lambda) ; q_prop<-1+v
  
  #######################################
  Data<-surv_time_generator(beta_prop[1],beta_prop[2],beta_prop[3],uk)
  mean_Data<-mean(Data$obs_time)
  dis_prop<-abs(mean_Data- mean_data1)
  
  ##################################
  while (dis_prop>epsilon){
    #Propose c:
    lb<-c-0.3  ; ub<-c+0.3
    lb<-max(0.05,lb) ; ub<-min(0.95,ub)
    uk = runif(1,lb,ub)
    
    #Propose beta|c:
    z2<-ifelse(bio>uk,1,0) ; z3<-z1*z2
    D<-data.frame(cbind(x,del,z1,z2,z3))
    beta.model<-coxph(Surv(x,del) ~ z1+z2+z3,data=D)
    beta_prop<-t(rmvnorm(1,beta,beta.model$var))
    
    #Compute q|c,beta:
    lambda<--log(1-uk)
    v<-rgamma(1,2,lambda) ; q_prop<-1+v
    
    #######################################
    Data<-surv_time_generator(beta_prop[1],beta_prop[2],beta_prop[3],uk)
    mean_Data<-mean(Data$obs_time)
    dis_prop<-abs(mean_Data- mean_data1)
    
  }
  
  r<-exp(lf1(uk,beta_prop,q_prop)-lf1(c,beta,q))
  u<-runif(1)
  while (u>r){
    #Propose c:
    lb<-c-0.3  ; ub<-c+0.3
    lb<-max(0.05,lb) ; ub<-min(0.95,ub)
    uk = runif(1,lb,ub)
    
    
    #Propose beta|c:
    z2<-ifelse(bio>uk,1,0) ; z3<-z1*z2
    D<-data.frame(cbind(x,del,z1,z2,z3))
    beta.model<-coxph(Surv(x,del) ~ z1+z2+z3,data=D)
    #beta_prop<-t(rmvnorm(1,beta.model$coefficients,beta.model$var))
    beta_prop<-t(rmvnorm(1,beta,beta.model$var))
    
    #Propose q|c,beta:
    lambda<--log(1-uk)
    v<-rgamma(1,2,lambda) ; q_prop<-1+v
    
    #######################################
    Data<-surv_time_generator(beta_prop[1],beta_prop[2],beta_prop[3],uk)
    mean_Data<-mean(Data$obs_time)
    dis_prop<-abs(mean_Data- mean_data1)
    
    ##################################
    while (dis_prop>epsilon){
      #Propose c:
      lb<-c-0.3  ; ub<-c+0.3
      lb<-max(0.05,lb) ; ub<-min(0.95,ub)
      uk = runif(1,lb,ub)
      
      
      #Propose beta|c:
      z2<-ifelse(bio>uk,1,0) ; z3<-z1*z2
      D<-data.frame(cbind(x,del,z1,z2,z3))
      beta.model<-coxph(Surv(x,del) ~ z1+z2+z3,data=D)
      beta_prop<-t(rmvnorm(1,beta,beta.model$var))
      
      #Propose q|c,beta:
      lambda<--log(1-uk)
      v<-rgamma(1,2,lambda) ; q_prop<-1+v
      
      #######################################
      Data<-surv_time_generator(beta_prop[1],beta_prop[2],beta_prop[3],uk)
      mean_Data<-mean(Data$obs_time)
      dis_prop<-abs(mean_Data- mean_data1)
    }
    
    r<-exp(lf1(uk,beta_prop,q_prop)-lf1(c,beta,q))
    u<-runif(1) 
  }
  
  c<-uk ; beta<-beta_prop ; q<-q_prop
  c.v[s]<-c  
  beta.m[s,]<-t(beta)
  q.v[s]<-q  
}


############################################  

par(mfrow=c(2,2))
plot(c(1:length(c.v)),c.v,type="l",main="Trace plot of c at eps=0.01, S=3000",xlab="iteration",ylab="c")
plot(c(1:length(c.v)),beta.m[,1],type="l",main="Trace plot of beta1, S=3000",xlab="iteration",ylab="beta1")
plot(c(1:length(c.v)),beta.m[,2],type="l",main="Trace plot of beta2, S=3000",xlab="iteration",ylab="beta2")
plot(c(1:length(c.v)),beta.m[,3],type="l",main="Trace plot of beta3, S=3000",xlab="iteration",ylab="beta3")

###############################################

c.v<-as.matrix(c.v)
c.v<-as.vector(submatrix(c.v,seq(B+1,S,by=2),1))  
beta.m<-submatrix(beta.m, seq(B+1,S,by=2), c(1:3))
Gibbs_Sample<-data.frame(cbind(c.v,beta.m))
colnames(Gibbs_Sample) <- c("c","beta1","beta2","beta3")

# change the file name BELOW!
write.table(Gibbs_Sample,file="data(0.2,0,0.4,0).txt") 
#D=read.table("data(0.2,1.5,1.0).txt",header=TRUE) # says first column are rownames

#Record the estimates of c:
hat.c<-mean(Gibbs_Sample$c) #0.2060863
cred.int.c<-quantile(Gibbs_Sample$c,c(.025,.975)) # 2.5%     97.5% 
#0.07136967 0.38844485 

#Marginal method for the estimation of betas:
hat.beta<-apply(Gibbs_Sample[,2:4],2,mean) # 0.15236342 0.35791620 0.01279727  vs target: 0  0.4  0
