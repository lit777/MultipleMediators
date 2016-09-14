#------ Load the main dataset
load("Master.RData")

#------ Load R functions
source("Correlation_Alt.R")
source("MetropolisHastingsAlgorithm.R")

#------ Set Treatment (TRT), Outcome (OUT), Mediators (M) and Covariates (X)
Data <- Master
OUT <- Data$PM.2.5
TRT <- Data$SO2.SC
M <- cbind(Data$SO2_Annual, Data$NOx_Annual, Data$CO2_Annual)
XX <- cbind(Data$S_n_CR, Data$NumNOxControls, Data$Heat_Input/100000, Data$Barometric_Pressure, Data$Temperature,  Data$PctCapacity, Data$sulfur_Content, Data$Phase2_Indicator, Data$Operating_Time/1000)

dim.cov <- dim(XX)[2] #<--------- Num. of Covariates
dimx <- dim.cov+1

#------ Variables by treatments
x0 <- XX[which(TRT==0),]
x1 <- XX[which(TRT==1),]

y0 <- OUT[which(TRT==0)]
y1 <- OUT[which(TRT==1)]

m0 <- log(M[which(TRT==0),])
m1 <- log(M[which(TRT==1),])

n0 <- dim(x0)[1]
n1 <- dim(x1)[1]

#------- load required libraries
library(mnormt)
library(gtools)
library(numDeriv)
library(matrixcalc)
library(corpcor)
library(rootSolve)

P <- 8 # number of marginal distributions
K <- 9 # numner of clusters



#-------- Initial Settings
MCMC <- 100000 # Num. of Iterations

### Set matrices to place posteriors ###
para.y1 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov) # Y(1)
# Set initial values
para.y1[1,] <- para.y1[2,] <- c(NA,rep(0.5,(K-1)),rep(-6,K),coefficients(lm(y1~x1))[-1],rep(var(y1)/4,K),2,2,2,mean(y1),1,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K))

para.y0 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov) # Y(0)
para.y0[1,] <- para.y0[2,] <- c(NA,rep(0.5,(K-1)),rep(coefficients(lm(y0~x0))[1],K),coefficients(lm(y0~x0))[-1],rep(var(y0)/4,K),2,2,2,mean(y0),1,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K))

para.m11 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov) # Y(0)
para.m11[1,] <- para.m11[2,] <- c(NA,rep(0.5,(K-1)),rep(coefficients(lm(m1[,1]~x1))[1],K),coefficients(lm(m1[,1]~x1))[-1],rep(var(m1[,1])/4,K),2,2,2,mean(m1[,1]),1,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K))

para.m21 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov) # m2(1)
para.m21[1,] <- para.m21[2,] <- c(NA,rep(0.5,(K-1)),rep(coefficients(lm(m1[,2]~x1))[1],K),coefficients(lm(m1[,2]~x1))[-1],rep(var(m1[,2])/4,K),2,2,2,mean(m1[,2]),1,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K))

para.m31 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov) # m3(1)
para.m31[1,] <- para.m31[2,] <- c(NA,rep(0.5,(K-1)),rep(coefficients(lm(m1[,3]~x1))[1],K),coefficients(lm(m1[,3]~x1))[-1],rep(var(m1[,3])/4,K),2,2,2,mean(m1[,3]),1,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K))

para.m10 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov) # m1(0)
para.m10[1,] <- para.m10[2,] <- c(NA,rep(0.5,(K-1)),rep(coefficients(lm(m0[,1]~x0))[1],K),coefficients(lm(m0[,1]~x0))[-1],rep(var(m0[,1])/4,K),2,2,2,mean(m0[,1]),1,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K))

para.m20 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov) # m2(0)
para.m20[1,] <- para.m20[2,] <- c(NA,rep(0.5,(K-1)),rep(coefficients(lm(m0[,2]~x0))[1],K),coefficients(lm(m0[,2]~x0))[-1],rep(var(m0[,2])/4,K),2,2,2,mean(m0[,2]),1,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K))

para.m30 <- matrix(nrow = MCMC, ncol = 6+6*K+2*dim.cov) # m3(0)
para.m30[1,] <- para.m30[2,] <- c(NA,rep(0.5,(K-1)),rep(coefficients(lm(m0[,3]~x0))[1],K),coefficients(lm(m0[,3]~x0))[-1],rep(var(m0[,3])/4,K),2,2,2,mean(m0[,3]),1,NA,rep(NA,dim.cov),rep(NA,K),rep(NA,K),rep(NA,K))

para.C <- matrix(nrow = MCMC, ncol = (28*2+3)) # Correlations 
para.C[1,] <- para.C[2,] <- c(rep(NA,28),c(0.130123426029935, -0.0854405048744152, -0.0837929602888095,
                    0, 0, 0, 0, 0.528898227751186, 0.517242870241654, 0, 0, 0, 0,
                    0.804514579489972, 0, 0, 0, 0, 0, 0, 0, 0, 0.125181447538096,
                    -0.0769838292735945, 0.00400666063036421, 0.689859929372401,
                    0.701074907893295, 0.835553284184042),c(0.001,0.001,0.001))

# Indices of parameters
ind1 <- (K+1):(2*K) # intercept
ind2 <- (2*K+1):(2*K+dim.cov) # regression coefficient except the intercept
ind3 <- (2*K+dim.cov+1):(3*K+dim.cov) # variance 
ind4 <- 3*K+dim.cov+1 # alpha
ind5 <- ind4+1 # alpha_beta0
ind6 <- ind5+1 # alpha_sigma
ind7 <- ind6+1 # mu_beta0
ind8 <- ind7+1 # sigma_beta0

# Initial values for Gaussian variables for the Copula model
h <- rbind(cbind(y1,m1[,1:3],0,0,0,0),
           cbind(0,0,0,0,y0,m0[,1:3]))

# Initial setting for R
R <- diag(1,P)

# Starting values for variances in adaptive sampler (1)
cov.y1 <- log(sqrt(diag(vcov(lm(y1~x1))[1,1],K)/4))
cov.m11 <- log(sqrt(diag(vcov(lm(m1[,1]~x1))[1,1],K)/4))
cov.m21 <- log(sqrt(diag(vcov(lm(m1[,2]~x1))[1,1],K)/4))
cov.m31 <- log(sqrt(diag(vcov(lm(m1[,3]~x1))[1,1],K)/4))
cov.y0 <- log(sqrt(diag(vcov(lm(y0~x0))[1,1],K)/4))
cov.m10 <- log(sqrt(diag(vcov(lm(m0[,1]~x0))[1,1],K)/4))
cov.m20 <- log(sqrt(diag(vcov(lm(m0[,2]~x0))[1,1],K)/4))
cov.m30 <- log(sqrt(diag(vcov(lm(m0[,3]~x0))[1,1],K)/4))

# Starting values for variances in adaptive sampler (2)
cov2.y1 <- log(sqrt(diag(diag(vcov(lm(y1~x1))[-1,-1]),K)/8))
cov2.m11 <- log(sqrt(diag(diag(vcov(lm(m1[,1]~x1))[-1,-1]),K)/8))
cov2.m21 <- log(sqrt(diag(diag(vcov(lm(m1[,2]~x1))[-1,-1]),K)/8))
cov2.m31 <- log(sqrt(diag(diag(vcov(lm(m1[,3]~x1))[-1,-1]),K)/8))
cov2.y0 <- log(sqrt(diag(diag(vcov(lm(y0~x0))[-1,-1]),K)/8))
cov2.m10 <- log(sqrt(diag(diag(vcov(lm(m0[,1]~x0))[-1,-1]),K)/8))
cov2.m20 <- log(sqrt(diag(diag(vcov(lm(m0[,2]~x0))[-1,-1]),K)/8))
cov2.m30 <- log(sqrt(diag(diag(vcov(lm(m0[,3]~x0))[-1,-1]),K)/8))

# Initial values for complete data: Y(1),Y(0),M(1,1,1),M(0,0,0)
y1 <- c(y1, rnorm(n0, mean(y1), sd(y1)))
y0 <- c(rnorm(n1, mean(y0), sd(y0)), y0)
m1 <- rbind(m1, rmnorm(n0, apply(m1, 2, mean), var(m1)))
m0 <- rbind(rmnorm(n1, apply(m0, 2, mean), var(m0)), m0)



#-------- Run MCMC
pb <- txtProgressBar(min = 0, max = MCMC, style = 3)
for (t in 3:MCMC){
  
  # Break up the MCMC run into several batches (50 iterations each) 
  # to monitor and manipulate the acceptance rates for the adaptive samplers
  SEQ <- seq(54, MCMC, by=50)

  #### Y(1) ####
  if(t %in% SEQ){
    for(c in 1:K){
      if(mean(para.y1[(t-51):(t-1),(dim(para.y1)[2]-3*K+c)]) < 0.44 ){
        cov.y1[c,c] <- cov.y1[c,c]-min(0.01, 1/sqrt(t)) # reduce the variance by min(0.01, 1/sqrt(t))
      }else{
        cov.y1[c,c] <- cov.y1[c,c]+min(0.01, 1/sqrt(t)) # increase the variance by min(0.01, 1/sqrt(t))
      }
    }
  }
  
  if(t %in% SEQ){
    for(c in 1:K){
      if(mean(para.y1[(t-51):(t-1),(dim(para.y1)[2]-1*K+c)]) < 0.44 ){
        cov2.y1[c,c] <- cov2.y1[c,c]-min(0.01, 1/sqrt(t)) # reduce the variance by min(0.01, 1/sqrt(t))
      }else{
        cov2.y1[c,c] <- cov2.y1[c,c]+min(0.01, 1/sqrt(t)) # increase the variance by min(0.01, 1/sqrt(t))
      }
    }
  }
  
  para.y1[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=y1, R=R, w_pre=para.y1[t-1,2:K],
                         beta0_pre=para.y1[t-1,ind1], beta_pre=para.y1[t-1,ind2],sigma_pre=para.y1[t-1,ind3],
                         alpha_pre=para.y1[t-1,ind4], alpha_beta0_pre=para.y1[t-1,ind5], alpha_sigma_pre=para.y1[t-1,ind6],
                         mu_beta0_pre=para.y1[t-1,ind7], sigma_beta0_pre=para.y1[t-1,ind8],K=K, eps=0.1,
                         del1=15, del2=10, del3=30, index=1, cov1=cov.y1,cov2=cov2.y1, zz=1)

  # Update a Gaussian variable for the Copula model
  prop1 <- para.y1[t,2:K]
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])
  pprop.y1 <- pprop1
    
  h[,1] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(c(y1),rbind(x1,x0)), 1, function(z) 
    sum(pprop1*pnorm(z[1],para.y1[t,ind1]+sum(para.y1[t,ind2]*z[2:(dim(x0)[2]+1)]), sqrt(para.y1[t,ind3])))))))
  

  #### Y(0) ####
  if(t %in% SEQ){
    for(c in 1:K){
      if(mean(para.y0[(t-51):(t-1),(dim(para.y0)[2]-3*K+c)]) < 0.44 ){
        cov.y0[c,c] <- cov.y0[c,c]-min(0.01, 1/sqrt(t))
      }else{
        cov.y0[c,c] <- cov.y0[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  
  if(t %in% SEQ){
    for(c in 1:K){
      if(mean(para.y0[(t-51):(t-1),(dim(para.y0)[2]-1*K+c)]) < 0.44 ){
        cov2.y0[c,c] <- cov2.y0[c,c]-min(0.01, 1/sqrt(t))
      }else{
        cov2.y0[c,c] <- cov2.y0[c,c]+min(0.01, 1/sqrt(t)) }
    }
  }
  
  para.y0[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=y0, R=R, w_pre=para.y0[t-1,2:K], 
                         beta0_pre=para.y0[t-1,ind1], beta_pre=para.y0[t-1,ind2], sigma_pre=para.y0[t-1,ind3],
                         alpha_pre=para.y0[t-1,ind4], alpha_beta0_pre=para.y0[t-1,ind5], alpha_sigma_pre=para.y0[t-1,ind6],
                         mu_beta0_pre=para.y0[t-1,ind7], sigma_beta0_pre=para.y0[t-1,ind8], K=K, eps=0.1,
                         del1=15, del2=10, del3=30, index=5, cov1=cov.y0, cov2=cov2.y0, zz=0)
  
  prop1 <- para.y0[t,2:K]
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])
  pprop.y0 <- pprop1

  h[,5] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(c(y0),rbind(x1,x0)), 1, function(z) 
    sum(pprop1*pnorm(z[1],para.y0[t,ind1]+sum(para.y0[t,ind2]*z[2:(dim(x0)[2]+1)]),sqrt(para.y0[t,ind3])))))))
  
  
  #### M1(0) ####  
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m10[(t-51):(t-1),(dim(para.m10)[2]-3*K+c)]) < 0.44 ){
        cov.m10[c,c] <- cov.m10[c,c]-min(0.01, 1/sqrt(t))
      }else{
        cov.m10[c,c] <- cov.m10[c,c]+min(0.01, 1/sqrt(t))
      }
    }
  }
  
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m10[(t-51):(t-1),(dim(para.m10)[2]-1*K+c)]) < 0.44 ){
        cov2.m10[c,c] <- cov2.m10[c,c]-min(0.01, 1/sqrt(t)) 
      }else{
        cov2.m10[c,c] <- cov2.m10[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  
  para.m10[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m0[,1], R=R, w_pre=para.m10[t-1,2:K],
                          beta0_pre=para.m10[t-1,ind1], beta_pre=para.m10[t-1,ind2], sigma_pre=para.m10[t-1,ind3],
                          alpha_pre=para.m10[t-1,ind4], alpha_beta0_pre=para.m10[t-1,ind5], alpha_sigma_pre=para.m10[t-1,ind6],
                          mu_beta0_pre=para.m10[t-1,ind7], sigma_beta0_pre=para.m10[t-1,ind8], K=K, eps=0.1,
                          del1=15, del2=10, del3=20, index=6, cov1=cov.m10, cov2=cov2.m10, zz=0)
  
  prop1 <- para.m10[t,2:K]
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])
  pprop.m10 <- pprop1

  h[ ,6]<-   qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(c(m0[,1]),rbind(x1,x0)), 1, function(z) 
    sum(pprop1*pnorm(z[1],para.m10[t,ind1]+sum(para.m10[t,ind2]*z[2:(dim(x0)[2]+1)]),sqrt(para.m10[t,ind3])))))))
  
  
  #### M2(0) ####    
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m20[(t-51):(t-1),(dim(para.m20)[2]-3*K+c)]) < 0.44 ){
        cov.m20[c,c] <- cov.m20[c,c]-min(0.01, 1/sqrt(t))
      }else{
        cov.m20[c,c] <- cov.m20[c,c]+min(0.01, 1/sqrt(t))
      }
    }
  }
  
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m20[(t-51):(t-1),(dim(para.m20)[2]-1*K+c)]) < 0.44 ){
        cov2.m20[c,c] <- cov2.m20[c,c]-min(0.01, 1/sqrt(t))
      }else{
        cov2.m20[c,c] <- cov2.m20[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  
  para.m20[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m0[,2], R=R, w_pre=para.m20[t-1,2:K],
                          beta0_pre=para.m20[t-1,ind1], beta_pre=para.m20[t-1,ind2], sigma_pre=para.m20[t-1,ind3],
                          alpha_pre=para.m20[t-1,ind4], alpha_beta0_pre=para.m20[t-1,ind5], alpha_sigma_pre=para.m20[t-1,ind6],
                          mu_beta0_pre=para.m20[t-1,ind7], sigma_beta0_pre=para.m20[t-1,ind8], K=K, eps=0.1,
                          del1=15, del2=10, del3=20, index=7, cov1=cov.m20, cov2=cov2.m20, zz=0)
  
  prop1 <- para.m20[t,2:K]
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])
  pprop.m20 <- pprop1

  h[ ,7] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(c(m0[,2]),rbind(x1,x0)), 1, function(z)
    sum(pprop1*pnorm(z[1],para.m20[t,ind1]+sum(para.m20[t,ind2]*z[2:(dim(x0)[2]+1)]),sqrt(para.m20[t,ind3])))))))
  
  #### M3(0) ####      
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m30[(t-51):(t-1),(dim(para.m30)[2]-3*K+c)]) < 0.44 ){
        cov.m30[c,c] <- cov.m30[c,c]-min(0.01, 1/sqrt(t))
      }else{
        cov.m30[c,c] <- cov.m30[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m30[(t-51):(t-1),(dim(para.m30)[2]-1*K+c)]) < 0.44 ){
        cov2.m30[c,c] <- cov2.m30[c,c]-min(0.01, 1/sqrt(t)) 
      }else{
        cov2.m30[c,c] <- cov2.m30[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  
  para.m30[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m0[,3], R=R, w_pre=para.m30[t-1,2:K],
                          beta0_pre=para.m30[t-1,ind1], beta_pre=para.m30[t-1,ind2], sigma_pre=para.m30[t-1,ind3],
                          alpha_pre=para.m30[t-1,ind4], alpha_beta0_pre=para.m30[t-1,ind5], alpha_sigma_pre=para.m30[t-1,ind6],
                          mu_beta0_pre=para.m30[t-1,ind7], sigma_beta0_pre=para.m30[t-1,ind8], K=K, eps=0.1,
                          del1=15, del2=10, del3=20, index=8, cov1=cov.m30, cov2=cov2.m30, zz=0)
  
  prop1 <- para.m30[t,2:K]
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])
  pprop.m30 <- pprop1

  h[ ,8] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(c(m0[,3]),rbind(x1,x0)), 1, function(z) 
    sum(pprop1*pnorm(z[1],para.m30[t,ind1]+sum(para.m30[t,ind2]*z[2:(dim(x0)[2]+1)]),sqrt(para.m30[t,ind3])))))))
  
  
  
  #### M1(1) ####        
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m11[(t-51):(t-1),(dim(para.m11)[2]-3*K+c)]) < 0.44 ){
        cov.m11[c,c] <- cov.m11[c,c]-min(0.01, 1/sqrt(t)) 
      }else{
        cov.m11[c,c] <- cov.m11[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m11[(t-51):(t-1),(dim(para.m11)[2]-1*K+c)]) < 0.44 ){
        cov2.m11[c,c] <- cov2.m11[c,c]-min(0.01, 1/sqrt(t)) 
      }else{
        cov2.m11[c,c] <- cov2.m11[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  
  para.m11[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m1[,1], R=R, w_pre=para.m11[t-1,2:K],
                          beta0_pre=para.m11[t-1,ind1], beta_pre=para.m11[t-1,ind2], sigma_pre=para.m11[t-1,ind3],
                          alpha_pre=para.m11[t-1,ind4], alpha_beta0_pre=para.m11[t-1,ind5], alpha_sigma_pre=para.m11[t-1,ind6],
                          mu_beta0_pre=para.m11[t-1,ind7], sigma_beta0_pre=para.m11[t-1,ind8], K=K, eps=0.1,
                          del1=15, del2=10, del3=20, index=2, cov1=cov.m11, cov2=cov2.m11, zz=1)
  
  prop1 <- para.m11[t,2:K]
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])
  pprop.m11 <- pprop1

  h[ ,2] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(c(m1[,1]),rbind(x1,x0)), 1, function(z) 
    sum(pprop1*pnorm(z[1],para.m11[t,ind1]+sum(para.m11[t,ind2]*z[2:(dim(x1)[2]+1)]),sqrt(para.m11[t,ind3])))))))
  
  
  
  #### M2(1) ####          
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m21[(t-51):(t-1),(dim(para.m21)[2]-3*K+c)]) < 0.44 ){
        cov.m21[c,c] <- cov.m21[c,c]-min(0.01, 1/sqrt(t))
      }else{
        cov.m21[c,c] <- cov.m21[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m21[(t-51):(t-1),(dim(para.m21)[2]-1*K+c)]) < 0.44 ){
        cov2.m21[c,c] <- cov2.m21[c,c]-min(0.01, 1/sqrt(t)) 
      }else{
        cov2.m21[c,c] <- cov2.m21[c,c]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  
  para.m21[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m1[,2], R=R, w_pre=para.m21[t-1,2:K],
                          beta0_pre=para.m21[t-1,ind1], beta_pre=para.m21[t-1,ind2], sigma_pre=para.m21[t-1,ind3],
                          alpha_pre=para.m21[t-1,ind4], alpha_beta0_pre=para.m21[t-1,ind5], alpha_sigma_pre=para.m21[t-1,ind6],
                          mu_beta0_pre=para.m21[t-1,ind7], sigma_beta0_pre=para.m21[t-1,ind8], K=K, eps=0.1,
                          del1=15, del2=10, del3=20, index=3, cov1=cov.m21, cov2=cov2.m21,zz=1)
  
  prop1 <- para.m21[t,2:K]
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])
  pprop.m21 <- pprop1
  
  h[ ,3] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(c(m1[,2]),rbind(x1,x0)), 1, function(z) 
    sum(pprop1*pnorm(z[1],para.m21[t,ind1]+sum(para.m21[t,ind2]*z[2:(dim(x1)[2]+1)]),sqrt(para.m21[t,ind3])))))))
  
  
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m31[(t-51):(t-1),(dim(para.m31)[2]-3*K+c)]) < 0.44 ){
        cov.m31[c,c] <- cov.m31[c,c]-min(0.01, 1/sqrt(t)) 
      }else{
        cov.m31[c,c] <- cov.m31[c,c]+min(0.01, 1/sqrt(t))
      }
    }
  }
  if(t %in% SEQ){
    for(c in 1:K){
      if( mean(para.m31[(t-51):(t-1),(dim(para.m31)[2]-1*K+c)]) < 0.44 ){
        cov2.m31[tk,tk] <- cov2.m31[tk,tk]-min(0.01, 1/sqrt(t))
      }else{
        cov2.m31[tk,tk] <- cov2.m31[tk,tk]+min(0.01, 1/sqrt(t)) 
      }
    }
  }
  
  para.m31[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m1[,3], R=R, w_pre=para.m31[t-1,2:K],
                          beta0_pre=para.m31[t-1,ind1], beta_pre=para.m31[t-1,ind2], sigma_pre=para.m31[t-1,ind3],
                          alpha_pre=para.m31[t-1,ind4], alpha_beta0_pre=para.m31[t-1,ind5], alpha_sigma_pre=para.m31[t-1,ind6],
                          mu_beta0_pre=para.m31[t-1,ind7], sigma_beta0_pre=para.m31[t-1,ind8], K=K, eps=0.1,
                          del1=15, del2=10, del3=20, index=4, cov1=cov.m31, cov2=cov2.m31, zz=1)
  
  prop1 <- para.m31[t,2:K]
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])
  pprop.m31 <- pprop1

  h[,4] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(c(m1[,3]),rbind(x1,x0)), 1, function(z) 
    sum(pprop1*pnorm(z[1],para.m31[t,ind1]+sum(para.m31[t,ind2]*z[2:(dim(x1)[2]+1)]),sqrt(para.m31[t,ind3])))))))
  
  # Correlation parameters
  para.C[t,1:59] <- metropolisC(h=h, rho=para.C[t-1,29:56], prho=para.C[t-1,57:59])
  prop1 <- para.C[t,29:56]
  
  # Update the correlation matrix R
  R <- matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],
                prop1[9],prop1[14],1,prop1[19:22],prop1[4],prop1[10],prop1[15],prop1[19],1,
                prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],
                prop1[6],prop1[12],prop1[17],prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],
                prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE)
  
  
  # Impute missing part of the Copula model based on R
  h[(1+n1):(n1+n0),1:4] <- t(apply(h[(1+n1):(n1+n0),5:8], 1, function(x) 
    rmnorm(1, R[1:4,5:8]%*%solve(R[5:8,5:8])%*%c(x[1],x[2],x[3],x[4]), R[1:4,1:4]-R[1:4,5:8]%*%solve(R[5:8,5:8])%*%t(R[1:4,5:8]))))
  h[1:n1,5:8] <- t(apply(h[1:n1,1:4], 1, function(x) 
    rmnorm(1, R[5:8,1:4]%*%solve(R[1:4,1:4])%*%c(x[1],x[2],x[3],x[4]), R[5:8,5:8]-R[5:8,1:4]%*%solve(R[1:4,1:4])%*%t(R[5:8,1:4]))))
 
 
  # Update missing part of Y(1), Y(0), M(1,1,1) and M(0,0,0) from the above
  clus.y1 <- apply(rmultinom(n0, 1, pprop.y1), 2, function(x) which(x==1)) # cluster membership
  y1[(n1+1):(n1+n0)] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(h[(1+n1):(n1+n0),1]))), mean = para.y1[t,ind1][clus.y1]+x0%*%para.y1[t,ind2], sd=sqrt(para.y1[t,ind3][clus.y1]) )
  clus.y0 <- apply(rmultinom(n1, 1, pprop.y0), 2, function(x) which(x==1)) # cluster membership
  y0[(1):(n1)] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(h[(1):(n1),5]))), mean = para.y0[t,ind1][clus.y0]+x1%*%para.y0[t,ind2], sd=sqrt(para.y0[t,ind3][clus.y0]) )

  clus.m11 <- apply(rmultinom(n0, 1, pprop.m11), 2, function(x) which(x==1)) # cluster membership
  m1[(n1+1):(n1+n0),1] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(h[(1+n1):(n1+n0),2]))), mean = para.m11[t,ind1][clus.m11]+x0%*%para.m11[t,ind2], sd=sqrt(para.m11[t,ind3][clus.m11]) )
  clus.m21 <- apply(rmultinom(n0, 1, pprop.m21), 2, function(x) which(x==1)) # cluster membership
  m1[(n1+1):(n1+n0),2] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(h[(1+n1):(n1+n0),3]))), mean = para.m21[t,ind1][clus.m21]+x0%*%para.m21[t,ind2], sd=sqrt(para.m21[t,ind3][clus.m21]) )
  clus.m31 <- apply(rmultinom(n0, 1, pprop.m31), 2, function(x) which(x==1)) # cluster membership
  m1[(n1+1):(n1+n0),3] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(h[(1+n1):(n1+n0),4]))), mean = para.m31[t,ind1][clus.m31]+x0%*%para.m31[t,ind2], sd=sqrt(para.m31[t,ind3][clus.m31]) )

  clus.m10 <- apply(rmultinom(n1, 1, pprop.m10), 2, function(x) which(x==1)) # cluster membership
  m0[(1):(n1),1] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(h[(1):(n1),6]))), mean = para.m10[t,ind1][clus.m10]+x1%*%para.m10[t,ind2], sd=sqrt(para.m10[t,ind3][clus.m10]) )
  clus.m20 <- apply(rmultinom(n1, 1, pprop.m20), 2, function(x) which(x==1)) # cluster membership
  m0[(1):(n1),2] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(h[(1):(n1),7]))), mean = para.m20[t,ind1][clus.m20]+x1%*%para.m20[t,ind2], sd=sqrt(para.m20[t,ind3][clus.m20]) )
  clus.m30 <- apply(rmultinom(n1, 1, pprop.m30), 2, function(x) which(x==1)) # cluster membership
  m0[(1):(n1),3] <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(h[(1):(n1),8]))), mean = para.m30[t,ind1][clus.m30]+x1%*%para.m30[t,ind2], sd=sqrt(para.m30[t,ind3][clus.m30]) )
  
  Sys.sleep(0.001)
  setTxtProgressBar(pb, t)
}

save.image("MCMCsamples.RData")


