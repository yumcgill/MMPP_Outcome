library(parallel)
library(HiddenMarkov)
options(mc.cores = 22)
set.seed(443)
state<-2
ini<-c(0.7,0.3)
Q<-matrix(c(0,1,3,0),byrow=T,nrow=state)
diag(Q)<--rowSums(Q)
be<-matrix(c(-1,1),byrow=F,ncol=state)
## be<-<-matrix(c(0.8,1),byrow=F,ncol=state)
simdata<-NULL
## number of subjects
N.sim<-50
n.obs<-100
## end of the obervation interval
T.obs<-5
lam.rate<-c(4,12)
##lam.rate<-c(8,8)
sigma<-1
for(m in 1:N.sim){
  y<-NULL
  # NULL indicates that we have no data at this point
  x <- mmpp(NULL, Q, delta=ini, lambda=lam.rate)
  x <- simulate(x, nsim=n.obs-1)
  
  xh<-x$ys
  tim.g<-x$tau
  del<-0
  for (j in 1:n.obs){
    y[j]<-rnorm(1,cbind(1)%*%be[,xh[j]],sigma)
    del<-c(del,tim.g[j]-tim.g[j-1])
  }
  
  id<-rep(m,n.obs)
  sub<-cbind(id=id[tim.g<T.obs],y=y[tim.g<T.obs],tim.g=tim.g[tim.g<T.obs],del=del[tim.g<T.obs],xh=xh[tim.g<T.obs])
  simdata<-as.data.frame(rbind(simdata,sub))
}

simdata$int<-1
