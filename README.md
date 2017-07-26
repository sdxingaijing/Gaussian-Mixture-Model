# Gaussian-Mixture-Model of multivariate gaussian(normal distribution)
#This model is estimated by MCMC
#input some packages for random distribution and matrix manipilation

library(matrixStats)
library(bayesm)
library(MCMCpack)

set.seed(100)
#generate mixture model
num.vars=3  #number of variables
n=1000       #number of observations
num.clusters=3 #number of clusters
seed=NULL
if(is.null(seed)){
  set.seed(runif(1,0,100))
} else {
  set.seed(seed)
}
data <- data.frame(matrix(NA, nrow=n, ncol=num.vars+1))

mu <- NULL
for(m in 1:num.vars){
  mu <- cbind(mu,rnorm(num.clusters, runif(1,-10,15), 5))
}
mix<-bayesm::rdirichlet(rep(1,num.clusters))
for (i in 1:n) {
  cluster <- sample(1:num.clusters, size = 1,prob=mix)
  data[i, 1] <- cluster
  for(j in 1:num.vars){
    data[i, j+1] <- rnorm(1, mu[cluster,j], 1)
  }
}
data$X1 <- factor(data$X1)
var.names <- paste("VAR",seq(1,ncol(data)-1), sep="")
names(data) <- c("cluster",var.names)

library(ggplot2)
ggplot(data, aes(x = VAR1, y = VAR2, color = cluster)) + geom_point()
ggplot(data, aes(x = VAR1, y = VAR3, color = cluster)) + geom_point()
ggplot(data, aes(x = VAR2, y = VAR3, color = cluster)) + geom_point()
hist(data$VAR1)
hist(data$VAR3)
# data generation is over


#Gaussian mixture model of MCMC , also called as finite mixture model
#the MCMC inference illustrated in the book of "Finite Mixture and Markov Switching Models" Sylvia, 2006
## Set the initial value

K<-3 #specify the number of clusters
data0<-as.matrix(data[,-1])
I<-dim(data0)[1]
J<-dim(data0)[2]

ppi<-rep(1,K)/K #initial value of pi, prior of mixing propotion
post_prob<-rep(1,K)/K #initial value of posterior classification probability\
alpha<-rep(10,K)#hyperparameter of dirichlet distribution
state<-sample(c(1:K),size = I,prob = ppi,replace = TRUE)#initial state of classification
n_state<-rep(0,K) # initial value of assignment number in each group

  
nburn<-100 # iterations of burn-in
smcmc<-100 # iterations after burn_in
swst<-rep(0,(nburn+smcmc)) #number of switched state in each MCMC iteration. The less number of switched states, the better convergency the model is. 

mu=matrix(0,K,J)            # initial value of mean value of each state
mug=matrix(0,nburn+smcmc,K*J)       # keeping mcmc sample of mu in each MCMC iteration, the check the convergency

sigma=matrix(100,K,J)         # initial value of variance sigma, each row of this matrix refers to the diagonal elements of covariance matrix
#sigmag=matrix(0,smcmc,1)      # keeping mcmc sample  
inverse_diag<-function(x) {
  nn<-nrow(x)
  matrix_x<-list()
  for(i in 1:nn){
    matrix_x[[i]]<-solve(diag(x[i,]))
  }
  return(matrix_x)
}
sigmai=inverse_diag(sigma)  # inverse of sigma

#prior distribution of mu and sigma
##  mu ~ N(u0,v0)   normal distribution
u0=matrix(0,J,1)      # mean
v0=diag(J)*1000     # variance
v0i=solve(v0)
v0iu0=v0i%*%u0

#  sigma ~ IW_m(f0,g0)  inverted wishart distribution
f0=J+2
g0i=diag(1,J)*f0     # scale matrix 

#Log-density of Multivariate-Normal
lndMvn=function(x,mu,sigma){
  sigmai=solve(sigma)
  ll=-log(2*pi)/2-log(det(sigma))/2-(diag((x-mu)%*%sigmai%*%(x-mu)))/2;
  return(ll)
}

GMM<-function(data0,I,J,K,state,mu,sigma,sigmai,alpha,v0i,v0iu0,fn,g0i){
  #step 1 for mu and sigma
  for(k in 1:K){
    n_state[k]<-sum(state==k)
    state_k<-(data0[state==k,])
    #mu ~ N(b_k_S,B_k_S)
    B_k_S<-solve(v0i+n_state[k]*sigmai[[k]]) 
    #print(B_k_S)
    b_k_S<-B_k_S%*%(v0iu0+(sigmai[[k]]*n_state[k])%*%(colMeans(state_k))) #posterior mean of mu
    mu[k,]<- b_k_S+chol(B_k_S)%*%rnorm(J,0,1)
    
    #sigma ~ iWishart(c_k_S,C_k_S)
    c_k_S<-f0+n_state[k]/2
    resid<-state_k-t(matrix(1,J,n_state[k])*mu[k,])
    C_k_S<-g0i+t(resid)%*%resid/2
    sigma[k,]<-diag(MCMCpack::riwish(c_k_S,C_k_S))
    #print(c(k,n_state[k]))
  }
  sigmai=inverse_diag(sigma)
  
  #step 2 estimation of ppi
  ppi<-bayesm::rdirichlet(n_state+alpha)
  
  #step 3 for state assignment estimation
  for(i in 1:I){
    for(k in 1:K){
      post_prob[k]<- ppi[k]*exp(lndMvn(data0[i,],mu[k,],diag(sigma[k,])))
    }
    state[i]<-sample(1:K,size = 1,prob = post_prob) 
  }
  return(list(state=state,mu=mu,sigma=sigma,sigmai=sigmai,ppi=ppi))
}

#MCMC
# burn in period
for(ite in 1:nburn){
  #record the state of previous iteration
  state0<-state
 
  output<-GMM(data0,I,J,K,state,mu,sigma,sigmai,alpha,v0i,v0iu0,fn,g0i)
  state<-output$state
  mu   <-output$mu
  sigma<-output$sigma
  sigmai<-output$sigmai
  ppi  <-output$ppi
  
  swst[ite]<-sum(state!=state0)
  mug[ite,]<-t(vec(t(mu)))
}
# after burn in
for(ite in 1:smcmc){
  #record the state of previous iteration
  state0<-state
  
  output<-GMM(data0,I,J,K,state,mu,sigma,sigmai,alpha,v0i,v0iu0,fn,g0i)
  state<-output$state
  mu   <-output$mu
  sigma<-output$sigma
  sigmai<-output$sigmai
  ppi  <-output$ppi
  
  swst[ite+nburn]<-sum(state!=state0)
  mug[ite+nburn,]<-t(vec(t(mu)))
}

data2<-as.data.frame(cbind(state,data0))
data2$state<-as.factor(data2$state)
library(ggplot2)
#true data
ggplot(data, aes(x=VAR1, y=VAR2, color=cluster)) + 
  geom_point(size=1, alpha=0.6)
#classification output
ggplot(data2, aes(x=VAR1, y=VAR2, color=state)) + 
  geom_point(size=1, alpha=0.6)
#trace of switched states in MCMC iterations
plot.ts(swst)
#trace of mu
plot.ts(mug)
