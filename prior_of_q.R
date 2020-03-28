set.seed(1)
q<-seq(1.01,1000,length=100000)
p_q<-(q-1)/(q*(q+1))

plot(q,p_q,type="l",col="black", ylab = "p(q)",ylim=c(0,0.2))
max(p_q) #0.1715725
abline(h=0,col="red")
abline(h=0.1715725,col="red")
abline(v=1,col="red")
abline(v=500,col="red")

Rej_Samp<-function(n){
  x<-NULL
  for (i in 1:n){
  u<-runif(1,1,500)
  x[i]<-500*u
  v<-runif(1)
  y<-0.1715725*v
  p_q<-(x[i]-1)/(x[i]*(x[i]+1))
  
  while (y>p_q){
    u<-runif(1,1,500)
    x[i]<-500*u
    v<-runif(1)
    y<-0.1715725*v
    p_q<-(x[i]-1)/(x[i]*(x[i]+1))
  }
  }
  return(x)
}
q<-Rej_Samp(10000)
q1<-q/500

hist(q1,breaks = 100,main="Prior sample generated from Rectangular RS",
     xlab="q")
q1<-data.frame(q1)
colnames(q1) <- c("q")
#setwd("C:/Users/lingyee.sze.UC/OneDrive - University of Calgary/STAT619/project")
write.table(q1,file="prior_sample.txt") 
#Prior_sample=read.table("prior_sample.txt",header=TRUE)