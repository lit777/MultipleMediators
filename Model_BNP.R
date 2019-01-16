library(truncnorm)
library(mnormt)
library(BART)
library(DPpackage)
library(condMVNorm)
#------ Load R functions
load("Master.RData")
source("Correlation_Alt1.R")
source("MetropolisHastingsAlgorithm.R")

#------ Set Treatment (TRT), Outcome (OUT), Mediators (M) and Covariates (X)
Data <- Master
OUT <- Data$PM.2.5
TRT <- Data$SO2.SC
M <- (cbind(Data$SO2_Annual/10000, Data$NOx_Annual/1000, Data$CO2_Annual/10000000))

XX <- (cbind(Data$NumNOxControls, log(Data$Heat_Input), Data$Barometric_Pressure/100, Data$Temperature,  Data$PctCapacity/100, Data$Sulfur_Content, log(Data$Operating_Time), log(Data$Heat_Rate)))


dim.cov <- dim(XX)[2] #<--------- Num. of Covariates
dimx <- dim.cov+1

#------ Variables by treatments
x0 <- XX[which(TRT==0),]
x1 <- XX[which(TRT==1),]

y0 <- OUT[which(TRT==0)]
y1 <- OUT[which(TRT==1)]

m0 <- (M[which(TRT==0),])
m1 <- (M[which(TRT==1),])

n0 <- dim(x0)[1]
n1 <- dim(x1)[1]

n <- n0+n1

#------- load required libraries
library(mnormt)
library(gtools)
library(numDeriv)
library(matrixcalc)
library(corpcor)
library(rootSolve)

P <- 6 # number of marginal distributions
K <- 8 # numner of clusters

KK <- 15

#-------- Initial Settings
MCMC <- 20000 # Num. of Iterations


para.m11 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov+1+n0+n1) # Y(0)
para.m11[1,] <- para.m11[2,] <- c(rep(1/(K),(K)),rep(coefficients(lm(M[,1]~TRT+XX))[1],K),coefficients(lm(M[,1]~TRT+XX))[-c(1,2)],rep(0.1,K),2,2,2,mean(m1[,1]),500,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K),1, sample(1:K, n0+n1, replace=TRUE))

para.m21 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov+1+n0+n1) # m2(1)
para.m21[1,] <- para.m21[2,] <- c(rep(1/(K),(K)),rep(coefficients(lm(M[,2]~TRT+XX))[1],K),coefficients(lm(M[,2]~TRT+XX))[-c(1,2)],rep(0.1,K),2,2,2,mean(m1[,2]),500,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K),1, sample(1:K, n0+n1, replace=TRUE))

para.m31 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov+1+n0+n1) # m3(1)
para.m31[1,] <- para.m31[2,] <- c(rep(1/(K),(K)),rep(coefficients(lm(M[,3]~TRT+XX))[1],K),coefficients(lm(M[,3]~TRT+XX))[-c(1,2)],rep(0.1,K),2,2,2,mean(m1[,3]),500,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K),1, sample(1:K, n0+n1, replace=TRUE))

para.m10 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov+1+n0+n1) # m1(0)
para.m10[1,] <- para.m10[2,] <- c(rep(1/(K),(K)),rep(coefficients(lm(M[,1]~TRT+XX))[1],K),coefficients(lm(M[,1]~TRT+XX))[-c(1,2)],rep(0.1,K),2,2,2,mean(m0[,1]),500,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K),1, sample(1:K, n0+n1, replace=TRUE))

para.m20 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov+1+n0+n1) # m2(0)
para.m20[1,] <- para.m20[2,] <- c(rep(1/(K),(K)),rep(coefficients(lm(M[,2]~TRT+XX))[1],K),coefficients(lm(M[,2]~TRT+XX))[-c(1,2)],rep(0.1,K),2,2,2,mean(m0[,2]),500,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K),1, sample(1:K, n0+n1, replace=TRUE))

para.m30 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov+1+n0+n1) # m3(0)
para.m30[1,] <- para.m30[2,] <- c(rep(1/(K),(K)),rep(coefficients(lm(M[,3]~TRT+XX))[1],K),coefficients(lm(M[,3]~TRT+XX))[-c(1,2)],rep(0.1,K),2,2,2,mean(m0[,3]),500,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K),1, sample(1:K, n0+n1, replace=TRUE))

para.C <- matrix(nrow = MCMC, ncol = (15*2+1)) # Correlations
para.C[1,] <- para.C[2,] <- c(rep(NA,15),rep(0,15),c(0.001))

# Indices of parameters
ind1 <- (K+1):(2*K) # intercept
ind2 <- (2*K+1):(2*K+dim.cov) # regression coefficient except the intercept
ind3 <- (2*K+dim.cov+1):(3*K+dim.cov) # variance
ind4 <- 3*K+dim.cov+1 # alpha
ind5 <- ind4+1 # alpha_beta0
ind6 <- ind5+1 # alpha_sigma
ind7 <- ind6+1 # mu_beta0
ind8 <- ind7+1 # sigma_beta0

ind9 <- length(para.m11[1,])

# Initial values for Gaussian variables for the Copula model
h <- rbind(cbind(m1[,1:3],0,0,0),
           cbind(0,0,0,m0[,1:3]))

# Initial setting for R
R <- diag(1,P)

# Starting values for variances in adaptive sampler (1)
cov.m11 <- diag(vcov(lm(m1[,1]~x1))[1,1],K)/4
cov.m21 <- diag(vcov(lm(m1[,2]~x1))[1,1],K)/4
cov.m31 <- diag(vcov(lm(m1[,3]~x1))[1,1],K)/4
cov.m10 <- diag(vcov(lm(m0[,1]~x0))[1,1],K)/4
cov.m20 <- diag(vcov(lm(m0[,2]~x0))[1,1],K)/4
cov.m30 <- diag(vcov(lm(m0[,3]~x0))[1,1],K)/4

# Starting values for variances in adaptive sampler (2)
cov2.m11 <- diag(diag(vcov(lm(m1[,1]~x1))[-1,-1]),dim.cov)/8
cov2.m21 <- diag(diag(vcov(lm(m1[,2]~x1))[-1,-1]),dim.cov)/8
cov2.m31 <- diag(diag(vcov(lm(m1[,3]~x1))[-1,-1]),dim.cov)/8
cov2.m10 <- diag(diag(vcov(lm(m0[,1]~x0))[-1,-1]),dim.cov)/8
cov2.m20 <- diag(diag(vcov(lm(m0[,2]~x0))[-1,-1]),dim.cov)/8
cov2.m30 <- diag(diag(vcov(lm(m0[,3]~x0))[-1,-1]),dim.cov)/8

cov3.m11 <- diag(var(log(abs(resid(lm(m1[,1]~x1))))),K)/4
cov3.m21 <- diag(var(log(abs(resid(lm(m1[,2]~x1))))),K)/4
cov3.m31 <- diag(var(log(abs(resid(lm(m1[,3]~x1))))),K)/4
cov3.m10 <- diag(var(log(abs(resid(lm(m0[,1]~x0))))),K)/4
cov3.m20 <- diag(var(log(abs(resid(lm(m0[,2]~x0))))),K)/4
cov3.m30 <- diag(var(log(abs(resid(lm(m0[,3]~x0))))),K)/4


# Initial values for complete data: Y(1),Y(0),M(1,1,1),M(0,0,0)
y1 <- c(y1, rnorm(n0,mean(y1),0.1))
y0 <- c(rnorm(n1,mean(y0),0.1),y0)
m1 <- rbind(m1, cbind(rtruncnorm(n0,a=0, b=Inf,mean(m1[,1]),0.1),rtruncnorm(n0,a=0, b=Inf,mean(m1[,2]),0.1),rtruncnorm(n0,a=0, b=Inf,mean(m1[,3]),0.1)))
m0 <- rbind(cbind(rtruncnorm(n1,a=0, b=Inf,mean(m0[,1]),0.1),rtruncnorm(n1,a=0, b=Inf,mean(m0[,2]),0.1),rtruncnorm(n1,a=0, b=Inf,mean(m0[,3]),0.1)), m0)

ME <- list()

# Initial DPM fitting
w1 <- cbind(y1,m1,m0,rbind(x1,x0))
wbar1 <- apply(w1,2,mean)
wcov1 <- var(w1)
prior1 <- list(a0=10,b0=1,nu1=22,nu2=22,s2=0.5*wcov1,m2=wbar1,psiinv2=2*solve(wcov1),tau1=6.01,tau2=2.01)
state <- NULL
mcmc1 <- list(nburn=500, nsave=100, nskip=0, ndisplay=100)
y1111 <- DPdensity(y=cbind(y1[1:n1],m1[1:n1,],m0[1:n1,],x1),prior=prior1,mcmc=mcmc1,state=state,status=TRUE)


w0 <- cbind(y0,m1,m0,rbind(x1,x0))
wbar0 <- apply(w0,2,mean)
wcov0 <- var(w0)
prior0 <- list(a0=10,b0=1,nu1=32,nu2=32,s2=0.5*wcov0,m2=wbar0,psiinv2=2*solve(wcov0),tau1=6.01,tau2=2.01)
state <- NULL
mcmc0 <- list(nburn=500, nsave=100, nskip=0, ndisplay=100)
y0000 <- DPdensity(y=cbind(y0[(1+n1):(n1+n0)],m1[(1+n1):(n1+n0),],m0[(1+n1):(n1+n0),],x0),prior=prior0,mcmc=mcmc0,state=state,status=TRUE)

Y_1111 <- Y_1000 <- Y_1011 <- Y_1101 <- Y_1110 <- Y_1001 <- Y_1100 <- Y_1010 <- Y_0000 <- ede_1 <- ede_0 <- eae_1 <- eae_0 <- list()

acc1 <- acc0 <-0

cc<-0

#-------- Run MCMC
pb <- txtProgressBar(min = 0, max = MCMC, style = 3)
for (t in 3:MCMC){
  
  w1 <- cbind(y1,m1,m0,rbind(x1,x0))
  wbar1 <- apply(w1,2,mean)
  wcov1 <- var(w1)
  prior1 <- list(a0=10,b0=1,nu1=20,nu2=20,s2=0.5*wcov1,m2=wbar1,psiinv2=2*solve(wcov1),tau1=6.01,tau2=2.01)
  mcmc1 <- list(nburn=0, nsave=1, nskip=0, ndisplay=100)
  
  w0 <- cbind(y0,m1,m0,rbind(x1,x0))
  wbar0 <- apply(w0,2,mean)
  wcov0 <- var(w0)
  prior0 <- list(a0=10,b0=1,nu1=20,nu2=20,s2=0.5*wcov0,m2=wbar0,psiinv2=2*solve(wcov0),tau1=6.01,tau2=2.01)
  mcmc0 <- list(nburn=0, nsave=1, nskip=0, ndisplay=100)
  
  state <- y1111$state;
  y1111 <- DPdensity(y=cbind(y1[1:n1],m1[1:n1,],m0[1:n1,],x1),ngrid=500,prior=prior1,mcmc=mcmc1,state=state,status=FALSE)
  state <- y0000$state;
  y0000 <- DPdensity(y=cbind(y0[(1+n1):(n1+n0)],m1[(1+n1):(n1+n0),],m0[(1+n1):(n1+n0),],x0),ngrid=500,prior=prior0,mcmc=mcmc0,state=state,status=FALSE)
  
  pi <- table(y1111$state$ss)
  Sigma <- list()
  Mu <- list()
  for(i in 1:length(pi)){
    sigma <- diag(length(wbar1))
    sigma[lower.tri(sigma, diag=TRUE)] <- y1111$state$sigmaclus[i,]
    sigma <- sigma + t(sigma) - diag(diag(sigma))
    Sigma[[i]] <- sigma
    Mu[[i]] <- y1111$state$muclus[i,]
  }
  
  pi <- round(table(y1111$state$ss)/n1,2)
  
  pi0 <- table(y0000$state$ss)
  Sigma0 <- list()
  Mu0 <- list()
  for(i in 1:length(pi0)){
    sigma0 <- diag(length(wbar0))
    sigma0[lower.tri(sigma0, diag=TRUE)] <- y0000$state$sigmaclus[i,]
    sigma0 <- sigma0 + t(sigma0) - diag(diag(sigma0))
    Sigma0[[i]] <- sigma0
    Mu0[[i]] <- y0000$state$muclus[i,]
  }
  
  p0 <- round(table(y0000$state$ss)/n0,2)
  
  
  SEQ <- seq(1, 20000, by=5) # Thinning by 5
  
  # Conditional MVN given the joint distribution parameters
  CondMVN<-function(mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE)
  {
    B <- sigma[dependent.ind, dependent.ind]
    C <- sigma[dependent.ind, given.ind, drop = FALSE]
    D <- sigma[given.ind, given.ind]
    CDinv <- C %*% solve(D)
    cMu <- mean[dependent.ind] + CDinv %*% t(X.given - matrix(mean[given.ind],nrow=dim(X.given)[1], ncol=dim(X.given)[2],byrow=TRUE))
    list(condMean = cMu)
  }
  
  if(t %in% SEQ){
    cc <- cc + 1
    den<-sapply(1:length(pi), function(i) dmnorm(cbind(m1,m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i])/
      rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m1,m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_1111[[cc]]<-rowSums(sapply(1:length(pi), function(i) CondMVN(mean=Mu[[i]] , sigma=Sigma[[i]]  , dependent=c(1), given=c(2:15),X.given=cbind(m1,m0,rbind(x1,x0)))$condMean)*(den1))
    
    den<-sapply(1:length(pi), function(i) dmnorm(cbind(m0,m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i])/
      rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m0,m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_1000[[cc]]<-rowSums(sapply(1:length(pi), function(i) CondMVN(mean=Mu[[i]] , sigma=Sigma[[i]]  , dependent=c(1), given=c(2:15),X.given=cbind(m0,m0,rbind(x1,x0)))$condMean)*(den1))
    
    den<-sapply(1:length(pi), function(i) dmnorm(cbind(m0[,1],m1[,2:3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i])/
      rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m0[,1],m1[,2:3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_1011[[cc]]<-rowSums(sapply(1:length(pi), function(i) CondMVN(mean=Mu[[i]] , sigma=Sigma[[i]]  , dependent=c(1), given=c(2:15),X.given=cbind(m0[,1],m1[,2:3],m0,rbind(x1,x0)))$condMean)*(den1))
    
    den<-sapply(1:length(pi), function(i) dmnorm(cbind(m1[,1],m0[,2],m1[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i])/
      rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m1[,1],m0[,2],m1[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_1101[[cc]]<-rowSums(sapply(1:length(pi), function(i) CondMVN(mean=Mu[[i]] , sigma=Sigma[[i]]  , dependent=c(1), given=c(2:15),X.given=cbind(m1[,1],m0[,2],m1[,3],m0,rbind(x1,x0)))$condMean)*(den1))
    
    den<-sapply(1:length(pi), function(i) dmnorm(cbind(m1[,1],m1[,2],m0[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i])/
      rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m1[,1],m1[,2],m0[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_1110[[cc]]<-rowSums(sapply(1:length(pi), function(i) CondMVN(mean=Mu[[i]] , sigma=Sigma[[i]]  , dependent=c(1), given=c(2:15),X.given=cbind(m1[,1],m1[,2],m0[,3],m0,rbind(x1,x0)))$condMean)*(den1))
    
    den<-sapply(1:length(pi), function(i) dmnorm(cbind(m0[,1],m1[,2],m0[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i])/
      rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m0[,1],m1[,2],m0[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_1010[[cc]]<-rowSums(sapply(1:length(pi), function(i) CondMVN(mean=Mu[[i]] , sigma=Sigma[[i]]  , dependent=c(1), given=c(2:15),X.given=cbind(m0[,1],m1[,2],m0[,3],m0,rbind(x1,x0)))$condMean)*(den1))
    
    den<-sapply(1:length(pi), function(i) dmnorm(cbind(m0[,1],m0[,2],m1[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i])/
      rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m0[,1],m0[,2],m1[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_1001[[cc]]<-rowSums(sapply(1:length(pi), function(i) CondMVN(mean=Mu[[i]] , sigma=Sigma[[i]]  , dependent=c(1), given=c(2:15),X.given=cbind(m0[,1],m0[,2],m1[,3],m0,rbind(x1,x0)))$condMean)*(den1))
    
    den<-sapply(1:length(pi), function(i) dmnorm(cbind(m1[,1],m0[,2],m0[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i])/
      rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m1[,1],m0[,2],m0[,3],m0,rbind(x1,x0)), Mu[[i]][-c(1)], Sigma[[i]][-c(1),-c(1)])*pi[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_1100[[cc]]<-rowSums(sapply(1:length(pi), function(i) CondMVN(mean=Mu[[i]] , sigma=Sigma[[i]]  , dependent=c(1), given=c(2:15),X.given=cbind(m1[,1],m0[,2],m0[,3],m0,rbind(x1,x0)))$condMean)*(den1))
    
    den<-sapply(1:length(pi0), function(i) dmnorm(cbind(m1,m0,rbind(x1,x0)), Mu0[[i]][-c(1)], Sigma0[[i]][-c(1),-c(1)])*pi0[i])/
      rowSums(sapply(1:length(pi0), function(i) dmnorm(cbind(m1,m0,rbind(x1,x0)), Mu0[[i]][-c(1)], Sigma0[[i]][-c(1),-c(1)])*pi0[i]))
    den1<-ifelse(is.na(den), 0, den)
    Y_0000[[cc]]<-rowSums(sapply(1:length(pi0), function(i) CondMVN(mean=Mu0[[i]], sigma=Sigma0[[i]], dependent=c(1), given=c(2:15),X.given=cbind(m1,m0,rbind(x1,x0)))$condMean)*(den1))
    
    
    ME[[cc]] <- list(m1=m1, m0=m0)
    
    points(cc,mean(unlist(Y_1111[[cc]]))-mean(unlist(Y_0000[[cc]])),col=1)
    points(cc,mean(unlist(Y_1111[[cc]]))-mean(unlist(Y_1000[[cc]])),col=2)
    points(cc,mean(unlist(Y_1111[[cc]]))-mean(unlist(Y_1011[[cc]])),col=3)
  }
  
  #### M1(0) ####
  if(t > 200){
    cov.m10 <- cov(para.m10[3:(t-1),(K+1):(2*K)])
  }
  
  if(t > 200){
    cov2.m10 <- cov(para.m10[3:(t-1),(2*K+1):(2*K+dim.cov)])
  }
  
  para.m10[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m0[,1], R=R, w_pre=para.m10[t-1,1:K],
                             beta0_pre=para.m10[t-1,ind1], beta_pre=para.m10[t-1,ind2], sigma_pre=para.m10[t-1,ind3],
                             alpha_pre=para.m10[t-1,ind4], alpha_beta0_pre=para.m10[t-1,ind5], alpha_sigma_pre=para.m10[t-1,ind6],
                             mu_beta0_pre=para.m10[t-1,ind7], sigma_beta0_pre=para.m10[t-1,ind8], K=K, eps=0.1,
                             del1=15, del2=10, del3=20, index=4, cov1=cov.m10, cov2=cov2.m10,cov3=cov3.m10, zz=0, gamma=para.m10[t-1,(ind9-(n0+n1))], Z=para.m10[t-1,(ind9-(n0+n1)+1):ind9])
  
  Zm10 <- para.m10[t,(ind9-(n0+n1)+1):ind9][(1+n1):(n1+n0)]
  h[(1+n1):(n1+n0),4] <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(m0[(1+n1):(n1+n0),1], a=0, b=Inf,para.m10[t,ind1][Zm10]+rbind(x0)%*%para.m10[t,ind2] ,sqrt(para.m10[t,ind3])[Zm10]))))
  
  
  #### M2(0) ####
  if(t > 200){
    cov.m20 <- cov(para.m20[3:(t-1),(K+1):(2*K)])
  }
  
  if(t > 200){
    cov2.m20 <- cov(para.m20[3:(t-1),(2*K+1):(2*K+dim.cov)])
  }
  
  para.m20[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m0[,2], R=R, w_pre=para.m20[t-1,1:K],
                             beta0_pre=para.m20[t-1,ind1], beta_pre=para.m20[t-1,ind2], sigma_pre=para.m20[t-1,ind3],
                             alpha_pre=para.m20[t-1,ind4], alpha_beta0_pre=para.m20[t-1,ind5], alpha_sigma_pre=para.m20[t-1,ind6],
                             mu_beta0_pre=para.m20[t-1,ind7], sigma_beta0_pre=para.m20[t-1,ind8], K=K, eps=0.1,
                             del1=15, del2=10, del3=20, index=5, cov1=cov.m20, cov2=cov2.m20,cov3=cov3.m20, zz=0, gamma=para.m20[t-1,(ind9-(n0+n1))],Z= para.m20[t-1,(ind9-(n0+n1)+1):ind9])
  
  Zm20 <- para.m20[t,(ind9-(n0+n1)+1):ind9][(1+n1):(n1+n0)]
  h[(1+n1):(n1+n0),5] <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(m0[(1+n1):(n1+n0),2], a=0, b=Inf,para.m20[t,ind1][Zm20]+rbind(x0)%*%para.m20[t,ind2] ,sqrt(para.m20[t,ind3])[Zm20]))))
  
  #### M3(0) ####
  if(t > 200){
    cov.m30 <- cov(para.m30[3:(t-1),(K+1):(2*K)])
  }
  
  if(t > 200){
    cov2.m30 <- cov(para.m30[3:(t-1),(2*K+1):(2*K+dim.cov)])
  }
  
  para.m30[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m0[,3], R=R, w_pre=para.m30[t-1,1:K],
                             beta0_pre=para.m30[t-1,ind1], beta_pre=para.m30[t-1,ind2], sigma_pre=para.m30[t-1,ind3],
                             alpha_pre=para.m30[t-1,ind4], alpha_beta0_pre=para.m30[t-1,ind5], alpha_sigma_pre=para.m30[t-1,ind6],
                             mu_beta0_pre=para.m30[t-1,ind7], sigma_beta0_pre=para.m30[t-1,ind8], K=K, eps=0.1,
                             del1=15, del2=10, del3=20, index=6, cov1=cov.m30, cov2=cov2.m30,cov3=cov3.m30, zz=0, gamma=para.m30[t-1,(ind9-(n0+n1))],  Z=para.m30[t-1,(ind9-(n0+n1)+1):ind9])
  
  Zm30 <- para.m30[t,(ind9-(n0+n1)+1):ind9][(1+n1):(n1+n0)]
  h[(1+n1):(n1+n0),6] <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(m0[(1+n1):(n1+n0),3], a=0, b=Inf,para.m30[t,ind1][Zm30]+rbind(x0)%*%para.m30[t,ind2] ,sqrt(para.m30[t,ind3])[Zm30]))))
  
  #### M1(1) ####
  if(t > 200){
    cov.m11 <- cov(para.m11[3:(t-1),(K+1):(2*K)])
  }
  
  if(t > 200){
    cov2.m11 <- cov(para.m11[3:(t-1),(2*K+1):(2*K+dim.cov)])
  }
  
  para.m11[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m1[,1], R=R, w_pre=para.m11[t-1,1:K],
                             beta0_pre=para.m11[t-1,ind1], beta_pre=para.m11[t-1,ind2], sigma_pre=para.m11[t-1,ind3],
                             alpha_pre=para.m11[t-1,ind4], alpha_beta0_pre=para.m11[t-1,ind5], alpha_sigma_pre=para.m11[t-1,ind6],
                             mu_beta0_pre=para.m11[t-1,ind7], sigma_beta0_pre=para.m11[t-1,ind8], K=K, eps=0.1,
                             del1=15, del2=10, del3=20, index=1, cov1=cov.m11, cov2=cov2.m11, cov3=cov3.m11, zz=1, gamma=para.m11[t-1,(ind9-(n0+n1))],  Z=para.m11[t-1,(ind9-(n0+n1)+1):ind9])
  
  
  Zm11 <- para.m11[t,(ind9-(n0+n1)+1):ind9][(1):(n1)]
  h[(1):(n1),1] <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(m1[(1):(n1),1], a=0, b=Inf,para.m11[t,ind1][Zm11]+rbind(x1)%*%para.m11[t,ind2] ,sqrt(para.m11[t,ind3])[Zm11]))))
  
  
  #### M2(1) ####
  if(t > 200){
    cov.m21 <- cov(para.m21[3:(t-1),(K+1):(2*K)])
  }
  
  if(t > 200){
    cov2.m21 <- cov(para.m21[3:(t-1),(2*K+1):(2*K+dim.cov)])
  }
  
  para.m21[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m1[,2], R=R, w_pre=para.m21[t-1,1:K],
                             beta0_pre=para.m21[t-1,ind1], beta_pre=para.m21[t-1,ind2], sigma_pre=para.m21[t-1,ind3],
                             alpha_pre=para.m21[t-1,ind4], alpha_beta0_pre=para.m21[t-1,ind5], alpha_sigma_pre=para.m21[t-1,ind6],
                             mu_beta0_pre=para.m21[t-1,ind7], sigma_beta0_pre=para.m21[t-1,ind8], K=K, eps=0.1,
                             del1=15, del2=10, del3=20, index=2, cov1=cov.m21, cov2=cov2.m21,cov3=cov3.m21,zz=1, gamma=para.m21[t-1,(ind9-(n0+n1))],  Z=para.m21[t-1,(ind9-(n0+n1)+1):ind9])
  
  
  Zm21 <- para.m21[t,(ind9-(n0+n1)+1):ind9][(1):(n1)]
  h[(1):(n1),2] <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(m1[(1):(n1),2], a=0, b=Inf, para.m21[t,ind1][Zm21]+rbind(x1)%*%para.m21[t,ind2] ,sqrt(para.m21[t,ind3])[Zm21]))))
  
  if(t > 200){
    cov.m31 <- cov(para.m31[3:(t-1),(K+1):(2*K)])
  }
  
  if(t > 200){
    cov2.m31 <- cov(para.m31[3:(t-1),(2*K+1):(2*K+dim.cov)])
  }
  
  para.m31[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m1[,3], R=R, w_pre=para.m31[t-1,1:K],
                             beta0_pre=para.m31[t-1,ind1], beta_pre=para.m31[t-1,ind2], sigma_pre=para.m31[t-1,ind3],
                             alpha_pre=para.m31[t-1,ind4], alpha_beta0_pre=para.m31[t-1,ind5], alpha_sigma_pre=para.m31[t-1,ind6],
                             mu_beta0_pre=para.m31[t-1,ind7], sigma_beta0_pre=para.m31[t-1,ind8], K=K, eps=0.1,
                             del1=15, del2=10, del3=20, index=3, cov1=cov.m31, cov2=cov2.m31,cov3=cov3.m31, zz=1, gamma=para.m31[t-1,(ind9-(n0+n1))], Z= para.m31[t-1,(ind9-(n0+n1)+1):ind9])
  
  
  Zm31 <- para.m31[t,(ind9-(n0+n1)+1):ind9][(1):(n1)]
  h[(1):(n1),3] <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(m1[(1):(n1),3], a=0, b=Inf, para.m31[t,ind1][Zm31]+rbind(x1)%*%para.m31[t,ind2] ,sqrt(para.m31[t,ind3])[Zm31]))))
  
  
  # Correlation parameters
  para.C[t,1:31] <- metropolisC(h=h, rho=para.C[t-1,16:30], prho=para.C[t-1,31])
  prop1 <- para.C[t,16:30]
  
  # Update the correlation matrix R
  R <- matrix(c(1,prop1[1:5],prop1[1],1,prop1[6:9],prop1[2],prop1[6],1,prop1[10:12],prop1[3],
                prop1[7],prop1[10],1,prop1[13:14],prop1[4],prop1[8],prop1[11],prop1[13],1,
                prop1[15],prop1[5],prop1[9],prop1[12],prop1[14],prop1[15],1),6,6,byrow=TRUE)
  
  
  # Impute missing mediators
  
  h_prop <- h
  m1_prop <- m1
  m0_prop <- m0
  
  MNORM <- function(n = 1, mean = rep(0, d), varcov){
    sqrt.varcov <-  chol(varcov)
    d <- ncol(sqrt.varcov)
    return(drop(mean + t(matrix(rnorm(n * d), d, n)) %*% sqrt.varcov))
  }
  
  # Proposal distribution
  m1_prop[(n1+1):(n1+n0),] <- MNORM(n0,m1[(n1+1):(n1+n0),], cov(m1)/5)
  
  clus.m11 <- para.m11[t,(ind9-(n0+n1)+1):ind9][(n1+1):(n1+n0)]
  h_prop[(1+n1):(n1+n0),1] <- qnorm(pmin(1-0.1^16,pmax(0.1^100, ptruncnorm(m1_prop[(1+n1):(n1+n0),1], a=0, b=Inf,para.m11[t,ind1][clus.m11]+rbind(x0)%*%para.m11[t,ind2] ,sqrt(para.m11[t,ind3])[clus.m11]))))
  clus.m21 <- para.m21[t,(ind9-(n0+n1)+1):ind9][(n1+1):(n1+n0)]
  h_prop[(1+n1):(n1+n0),2] <- qnorm(pmin(1-0.1^16,pmax(0.1^100, ptruncnorm(m1_prop[(1+n1):(n1+n0),2], a=0, b=Inf,para.m21[t,ind1][clus.m21]+rbind(x0)%*%para.m21[t,ind2] ,sqrt(para.m21[t,ind3])[clus.m21]))))
  clus.m31 <- para.m31[t,(ind9-(n0+n1)+1):ind9][(n1+1):(n1+n0)]
  h_prop[(1+n1):(n1+n0),3] <- qnorm(pmin(1-0.1^16,pmax(0.1^100, ptruncnorm(m1_prop[(1+n1):(n1+n0),3], a=0, b=Inf,para.m31[t,ind1][clus.m31]+rbind(x0)%*%para.m31[t,ind2] ,sqrt(para.m31[t,ind3])[clus.m31]))))
  
  rat1 <- dmnorm(h_prop[(n1+1):(n1+n0),1:3], rep(0,3), R[1:3,1:3], log=TRUE)+
    log(pmax(0.1^100,dtruncnorm(m1_prop[(n1+1):(n1+n0),1], a=0, b=Inf, para.m11[t,ind1][clus.m11]+rbind(x0)%*%para.m11[t,ind2], sqrt(para.m11[t,ind3][clus.m11]))))+
    log(pmax(0.1^100,dtruncnorm(m1_prop[(n1+1):(n1+n0),2], a=0, b=Inf, para.m21[t,ind1][clus.m21]+rbind(x0)%*%para.m21[t,ind2], sqrt(para.m21[t,ind3][clus.m21]))))+
    log(pmax(0.1^100,dtruncnorm(m1_prop[(n1+1):(n1+n0),3], a=0, b=Inf, para.m31[t,ind1][clus.m31]+rbind(x0)%*%para.m31[t,ind2], sqrt(para.m31[t,ind3][clus.m31]))))+
    log(pmax(0.1^300,rowSums(sapply(1:length(pi0), function(i) dmnorm(cbind(y0[(n1+1):(n1+n0)],m1_prop[(n1+1):(n1+n0),],x0), Mu0[[i]][-c(5:7)], Sigma0[[i]][-c(5:7),-c(5:7)])*pi0[i]))/
               rowSums(sapply(1:length(pi0), function(i) dmnorm(cbind(m1_prop[(n1+1):(n1+n0),],x0), Mu0[[i]][-c(1,5:7)], Sigma0[[i]][-c(1,5:7),-c(1,5:7)])*pi0[i]))))+
    sapply((1+n1):(n1+n0), function(c) dmnorm(m1[c,1:3],m1_prop[(c),1:3], cov(m1)/5, log=TRUE))
  
  rat2 <- dmnorm(h[(n1+1):(n1+n0),1:3], rep(0,3), R[1:3,1:3], log=TRUE)+
    log(pmax(0.1^100,dtruncnorm(m1[(n1+1):(n1+n0),1], a=0, b=Inf, para.m11[t,ind1][clus.m11]+rbind(x0)%*%para.m11[t,ind2], sqrt(para.m11[t,ind3][clus.m11]))))+
    log(pmax(0.1^100,dtruncnorm(m1[(n1+1):(n1+n0),2], a=0, b=Inf, para.m21[t,ind1][clus.m21]+rbind(x0)%*%para.m21[t,ind2], sqrt(para.m21[t,ind3][clus.m21]))))+
    log(pmax(0.1^100,dtruncnorm(m1[(n1+1):(n1+n0),3], a=0, b=Inf, para.m31[t,ind1][clus.m31]+rbind(x0)%*%para.m31[t,ind2], sqrt(para.m31[t,ind3][clus.m31]))))+
    log(pmax(0.1^300,rowSums(sapply(1:length(pi0), function(i) dmnorm(cbind(y0[(n1+1):(n1+n0)],m1[(n1+1):(n1+n0),],x0), Mu0[[i]][-c(5:7)], Sigma0[[i]][-c(5:7),-c(5:7)])*pi0[i]))/
               rowSums(sapply(1:length(pi0), function(i) dmnorm(cbind(m1[(n1+1):(n1+n0),],x0), Mu0[[i]][-c(1,5:7)], Sigma0[[i]][-c(1,5:7),-c(1,5:7)])*pi0[i]))))+
    sapply((1+n1):(n1+n0), function(c) dmnorm(m1_prop[c,1:3],m1[(c),1:3], cov(m1)/5, log=TRUE))
  
  rat <- rat1 - rat2
  
  acc1 <- rep(0,n0)
  
  for(ii in (1):(n0)){    
    if(log(runif(1))>rat[ii] | is.na(rat[ii])){
      m1_prop[(ii+n1),] <- m1[(ii+n1),]
      h_prop[(ii+n1),1:3] <- h[(ii+n1),1:3]
    }else{
      m1[(ii+n1),] <- m1_prop[(ii+n1),]
      h[(ii+n1),1:3] <- h_prop[(ii+n1),1:3]
      acc1[ii] <- 1
    }
  }
  
  
  # Proposal Distribution
  m0_prop[(1):(n1),] <- MNORM(n1,m0[(1):(n1),], cov(m0)/5)
  
  clus.m10 <- para.m10[t,(ind9-(n0+n1)+1):ind9][(1):(n1)]
  h_prop[(1):(n1),4] <- qnorm(pmin(1-0.1^16,pmax(0.1^100, ptruncnorm(m0_prop[(1):(n1),1], a=0, b=Inf,para.m10[t,ind1][clus.m10]+rbind(x1)%*%para.m10[t,ind2] ,sqrt(para.m10[t,ind3])[clus.m10]))))
  clus.m20 <- para.m20[t,(ind9-(n0+n1)+1):ind9][(1):(n1)]
  h_prop[(1):(n1),5] <- qnorm(pmin(1-0.1^16,pmax(0.1^100, ptruncnorm(m0_prop[(1):(n1),2], a=0, b=Inf,para.m20[t,ind1][clus.m20]+rbind(x1)%*%para.m20[t,ind2] ,sqrt(para.m20[t,ind3])[clus.m20]))))
  clus.m30 <- para.m30[t,(ind9-(n0+n1)+1):ind9][(1):(n1)]
  h_prop[(1):(n1),6] <- qnorm(pmin(1-0.1^16,pmax(0.1^100, ptruncnorm(m0_prop[(1):(n1),3], a=0, b=Inf,para.m30[t,ind1][clus.m30]+rbind(x1)%*%para.m30[t,ind2] ,sqrt(para.m30[t,ind3])[clus.m30]))))
  
  rat1 <- dmnorm(h_prop[(1):(n1),4:6], rep(0,3), R[4:6,4:6], log=TRUE)+
    log(pmax(0.1^100,dtruncnorm(m0_prop[(1):(n1),1], a=0, b=Inf, para.m10[t,ind1][clus.m10]+rbind(x1)%*%para.m10[t,ind2], sqrt(para.m10[t,ind3][clus.m10]))))+
    log(pmax(0.1^100,dtruncnorm(m0_prop[(1):(n1),2], a=0, b=Inf, para.m20[t,ind1][clus.m20]+rbind(x1)%*%para.m20[t,ind2], sqrt(para.m20[t,ind3][clus.m20]))))+
    log(pmax(0.1^100,dtruncnorm(m0_prop[(1):(n1),3], a=0, b=Inf, para.m30[t,ind1][clus.m30]+rbind(x1)%*%para.m30[t,ind2], sqrt(para.m30[t,ind3][clus.m30]))))+
    log(pmax(0.1^300,rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(y1[(1):(n1)],m0_prop[(1):(n1),],x1), Mu[[i]][-c(2:4)], Sigma[[i]][-c(2:4),-c(2:4)])*pi[i]))/
               rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m0_prop[(1):(n1),],x1), Mu[[i]][-c(1:4)], Sigma[[i]][-c(1:4),-c(1:4)])*pi[i]))))+
    sapply((1):(n1), function(c) dmnorm(m0[c,1:3],m0_prop[c,1:3],cov(m0)/5, log=TRUE))
  
  rat2 <- dmnorm(h[(1):(n1),4:6], rep(0,3), R[4:6,4:6], log=TRUE)+
    log(pmax(0.1^100,dtruncnorm(m0[(1):(n1),1], a=0, b=Inf, para.m10[t,ind1][clus.m10]+rbind(x1)%*%para.m10[t,ind2], sqrt(para.m10[t,ind3][clus.m10]))))+
    log(pmax(0.1^100,dtruncnorm(m0[(1):(n1),2], a=0, b=Inf, para.m20[t,ind1][clus.m20]+rbind(x1)%*%para.m20[t,ind2], sqrt(para.m20[t,ind3][clus.m20]))))+
    log(pmax(0.1^100,dtruncnorm(m0[(1):(n1),3], a=0, b=Inf, para.m30[t,ind1][clus.m30]+rbind(x1)%*%para.m30[t,ind2], sqrt(para.m30[t,ind3][clus.m30]))))+
    log(pmax(0.1^300,rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(y1[(1):(n1)],m0[(1):(n1),],x1), Mu[[i]][-c(2:4)], Sigma[[i]][-c(2:4),-c(2:4)])*pi[i]))/
               rowSums(sapply(1:length(pi), function(i) dmnorm(cbind(m0[(1):(n1),],x1), Mu[[i]][-c(1:4)], Sigma[[i]][-c(1:4),-c(1:4)])*pi[i]))))+
    sapply((1):(n1), function(c) dmnorm(m0_prop[c,1:3],m0[c,1:3], cov(m0)/5, log=TRUE))
  
  rat <- rat1 - rat2
  
  acc0 <- rep(0,n1)
  for(ii in (1):(n1)){    
    if(log(runif(1))>rat[ii] | is.na(rat[ii])){
      m0_prop[(ii),] <- m0[(ii),]
      h_prop[(ii),4:6] <- h[(ii),4:6]
    }else{
      m0[(ii),] <- m0_prop[(ii),]
      h[(ii),4:6] <- h_prop[(ii),4:6]
      acc0[ii] <- 1
    }
  }

  Sys.sleep(0.001)
  setTxtProgressBar(pb, t)
}

save(XX, para.m10,para.m20,para.m30,para.m11,para.m21,para.m31,para.C,Y_1111,Y_1000,Y_1011,Y_1101,Y_1110,Y_0000,Y_1010,Y_1001,Y_1100, ME,file="MCMC.RData")
