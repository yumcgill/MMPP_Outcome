library(expm)
library(gtools)
library(ctmcd)
## Forward backwards for MMPP and outcome process
## subdat is the dataset for one individual
## The window ends at time tend
## nu is the initial distribution
## In reality the first "event" happens at time 0 but since this is
## guaranteed it does not contribute towards the likelihood

forward.backward <- function(Qe, nu, beta,lam, subdat, tend) {
  n.stat <- length(nu)
  #n.cov <- nrow(beta)
  N <- nrow(subdat)
  QmL=Qe-diag(lam)
  beta.mult <- t(as.matrix(subdat[c("int")]) %*% beta)
  def <- mapply(dnorm, subdat$y, as.list(data.frame(beta.mult)))
  expm.N <- lapply(subdat$del, function(x) expm(QmL*x))
  
  
  A=matrix(ncol=n.stat,nrow=N+1)
  
  a = as.vector(nu) * def[,1] 
  a=a/sum(a)
  A[1,]=a
  for (i in 2:N) { ## do a^\top exp[(Q-Lambda)Deltatau] diag(lambda)
    a=as.numeric(t(a)%*%expm.N[[i]])*(def[,i] * lam)
    a=a/sum(a) +1e-300 ## renormalise so doesn't get too small
    A[i,]=a
  }
  
  Deltat=tend-subdat$tim.g[N]
  eta = expm(QmL*Deltat)
  a=as.numeric(t(a)%*%eta) ## no event at end
  a=a/sum(a)
  A[N+1,]=a

  list(a=A)
}

### Simulate the path backwards given A
SampleBackwards<-function(A,taus,tend,Q,lambda) {
  d<-nrow(Q)
  ntau=length(taus)
  QmL=Q-diag(lambda)
  
  states=rep(0,ntau+2)
  thisstate=sample(1:d,1,prob=A[ntau+2,])
  thistau=tend
  states[ntau+2]=thisstate
  
  for (j in ntau:1) {
    prevtau=taus[j]
    Deltat=thistau-prevtau
    probs=A[j+1,]*as.numeric(expm(QmL*Deltat)[,thisstate])
    probs=probs/sum(probs)
    prevstate=sample(1:d,1,prob=probs)
    
    thisstate=prevstate
    thistau=prevtau
    states[j+1]=thisstate
  }
  Deltat=thistau
  probs=A[1,]*as.numeric(expm(QmL*Deltat)[,thisstate])
  states[1]=sample(1:d,1,prob=probs)
  
  return(states)
}

### function to simulate the underlying Markov chain
markov.sim<-function(subdat,a,Qe,lam,tend){
  n.stat<-ncol(Qe)
  Qew<-cbind(Qe-diag(lam),lam)
  Qew<-rbind(Qew,rep(0,n.stat+1))
  new.tim<-c(subdat$del[-1],tend-subdat$tim.g[dim(subdat)[1]])
  pp<-lapply(1:(dim(subdat)[1]),function(j) {N <- array(0, dim=c(n.stat+1,n.stat+1))
  tpm <- expm(Qew*new.tim[j])
  N[a[j],a[j+1]] <- 1
  rNijTRiT_Unif(N, new.tim[j], gm=Qew,tpm)})
  time.R<-Reduce("+",lapply(pp, function(t){t$RiT}))[1:n.stat]
  time.N<-Reduce("+",lapply(pp, function(t){t$NijT}))[1:n.stat,1:n.stat]
  return(list(time.N=time.N,time.R=time.R))
}



### function to run one Gibbs iteration
OneGibbsSampler<-function(Q,lam,nu,bet,tend,split.data) {
  ## Priors
  alambda=1; blambda=alambda/8
  aQ=1; bQ=aQ/8
  
  n.stat <- length(lam)
  a.and.b <- 
    mclapply(split.data, forward.backward, 
             Qe=Q, nu=nu, beta=bet,lam=lam,tend=tend)
  
  split.a <- mclapply(a.and.b, `[[`, "a")
  ss <- mclapply(1:length(split.data.50), function(tt) {SampleBackwards(A=split.a[[tt]], taus=split.data[[tt]]$tim.g[-1],tend=tend,Q=Q,lambda=lam)})
  ss<-mclapply(ss,function(zz) factor(zz,levels=c(1:n.stat)))
  
  ## update initial probability
  ini.prob <- do.call(rbind, mclapply(split.a,function(x) x[1,]))
  sam.ini <- apply(ini.prob,1,function(x) sample(1:n.stat,1,prob=x))
  sam.ini<-factor(sam.ini,levels=c(1:n.stat))
  nu.update<-rdirichlet(1,table(sam.ini)+1)
  
  
  
  ## sample the path given the start and end states 
  time.RQ<-mclapply(1:length(split.data.50),function(tt){markov.sim(split.data[[tt]],a=ss[[tt]],Qe=Q,lam=lam,tend=tend)})
  N.n<-mclapply(time.RQ, function(x) x$time.N)
  R.n<-mclapply(time.RQ, function(x) x$time.R)
  
  
  ###update point rate lambda
  time.R.tot <- Reduce('+',R.n) + blambda
  count.end <- mclapply(mclapply(ss,function(x) x[-c(1,length(x))]),table)
  count.end.tot<-Reduce(`+`,count.end)+alambda
  lam.update<-rgamma(n.stat,shape=count.end.tot,rate=time.R.tot)
  
  ### update Q
  N.n1<-Reduce(`+`, N.n)
  R.n1<-Reduce(`+`, R.n)
  Q.update<-array(0,c(n.stat,n.stat))
  
  for (mmm in 1:(n.stat)){ 
    Q.update[mmm,-mmm]<-rgamma(n.stat-1,shape=N.n1[mmm,-mmm]+aQ,rate=R.n1[mmm]+bQ)
  }
  
  diag(Q.update) <- -rowSums(Q.update)
  
  
  ### update outcome mean
  sss <- Reduce(c,mclapply(ss,function(x) x[-length(x)]))
  data<-do.call(rbind, split.data)
  beta.update<-array(NA,c(n.cov,n.stat))
  for (k in 1:n.stat){ if (nrow(data[sss==k,])!=0) {
    beta.update[,k]<-rnorm(1,sum(data$y[sss==k])/(length(data$y[sss==k])+1),1/(length(data$y[sss==k])+1))
  } else {
    beta.update[,k]<-rnorm(1,be[,k],1)
  }
  }
  
  return(list(Q=Q.update,bet=beta.update,nu=nu.update,lam=lam.update))
}


