source("functions.r")
### starting values

Q0<-matrix(c(0,1.5,1,0),byrow=T,nrow=state)
diag(Q0)<- -rowSums(Q0)


split.data<-function(data){
  data <- data[order(data$id, data$tim.g),]
  split.data <- split(data, data$id)
  split.data <- 
    lapply(split.data, function(x) data.frame(x, obs.num = 1:nrow(x)))
  return(split.data)
}
split.data.50<-split.data(simdata)


mcmc.sam.Q<-list()
mcmc.sam.Q[[1]]<-Q0
mcmc.sam.nu<-list()
mcmc.sam.nu[[1]]<-c(0.6,0.4)
mcmc.sam.beta<-list()
mcmc.sam.beta[[1]]<-matrix(c(-0.5,0.5),nrow=1)
mcmc.sam.lam<-list()
mcmc.sam.lam[[1]]<-c(5,10)

n.cov<-1
N.iter<-10000
##for Ex 1
T.obs <- 5

for(i in 2:N.iter){
  
  one.iter.result<-OneGibbsSampler(split.data.50, Q = mcmc.sam.Q[[i-1]], nu = mcmc.sam.nu[[i-1]],
                                   bet = mcmc.sam.beta[[i-1]], lam = mcmc.sam.lam[[i-1]],tend = T.obs)
  
  
  mcmc.sam.Q[[i]]<-one.iter.result$Q
  mcmc.sam.nu[[i]]<-one.iter.result$nu
  mcmc.sam.beta[[i]]<-one.iter.result$bet
  mcmc.sam.lam[[i]]<-one.iter.result$lam
  
  if (i/10==round(i/10)) {
    print("Iteration")
    print(i)
    print(one.iter.result$lam)
    print(one.iter.result$Q)
  }
  
} 

