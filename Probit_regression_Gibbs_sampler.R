#Probit regression
rm(list=ls())

#6.3
Divorce_Data<-read.table("divorce.dat.txt", header=T)
x<-Divorce_Data$X2
y<-Divorce_Data$X0
#Gibbs sampling
tau2.beta<-16
tau2.c<-16
S=100000
n=24
z.m<-matrix(nrow=S,ncol=n)
beta.v<-rep(0,S)
c.v<-rep(0,S)

z.m[1,]<-z<-rnorm(n,0,1)
beta.v[1]<-beta<-0
c.v[1]<-c<-0

set.seed(200)
for (s in 2:S){
  #sample z
  for (i in 1:n){
    if (y[i]==1){
      u<-runif(1,pnorm(c-beta*x[i]),1)
      z[i]<-beta*x[i]+qnorm(u)
    } 
    if (y[i]==0) {
    u<-runif(1,0,pnorm(c-beta*x[i]))
    z[i]<-beta*x[i]+qnorm(u)
    }
  }
 
  z.m[s,]<-z
  
  #sample beta
  mean<-sum(x*z)*((sum(x^2)+tau2.beta^(-2))^(-1))
  var<-(sum(x^2)+tau2.beta^(-2))^(-1)
  beta<-rnorm(1,mean,sqrt(var))
  beta.v[s]<-beta
  
  #sample c
  dat<-data.frame(cbind(z,y))
  a<-max(dat$z[dat$y==0]) 
  b<-min(dat$z[dat$y==1]) 
  v<-runif(1,pnorm(a/sqrt(tau2.c)),pnorm(b/sqrt(tau2.c)))
  c<-sqrt(tau2.c)*qnorm(v)              
  c.v[s]<-c
}
#Plot
par(mfrow=c(2,2))
hist(z.m[,1],main=expression(z[1]),prob=T,xlab="",breaks = 50)
hist(z.m[,24],main=expression(z[24]),prob=T,xlab="",breaks = 50)
hist(beta.v,main=expression(beta),prob=T,xlab="",breaks = 50)
hist(c.v,main=expression(c),prob=T,xlab="",breaks = 50)
#check S.eff
#if effective sample size of 1000 is acheived for 1 z[i], we add 1 to eff.count
library(coda)
s.eff.z<-rep(0,n)
eff.count<-0
for (i in 1:n){
  s.eff.z[i]<-effectiveSize(z.m[,i])
  if (s.eff.z[i]>=1000){
    eff.count<-eff.count+1
  }
}
eff.count #24
effectiveSize(beta.v) #5694.947
effectiveSize(c.v)   #4037.15 

#plot ACF
par(mfrow=c(2,2))
acf(z.m[,1])
acf(z.m[,n])
acf(beta.v)
acf(c.v)

#95% posterior c.i. for beta:
quantile(beta.v,c(.025,.975)) #0.1098507 0.6866497 
#Pr(beta>0|y,x)
mean(beta.v>0) #0.99927

