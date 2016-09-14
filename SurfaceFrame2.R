
##---------------------------------------------------------------
## Required libraries
##---------------------------------------------------------------

#----- Parallel computing
library(doParallel)

#----- Make clusters based on the number of CPU cores
cl<-makeCluster(64)
registerDoParallel(cl)
getDoParWorkers()

#----- Support parallel excution
library(foreach)



##---------------------------------------------------------------
## Extract MCMC samples from Stan outputs
##---------------------------------------------------------------

#----- Load MCMC samples and Data
load("Master.RData")
load("MCMCsamples.RData")
load("SurfaceFrame1.RData")

#------ Set Treatment (TRT), Outcome (OUT), Mediators (M) and Covariates (X)
Data <- Master
OUT <- Data$PM.2.5
TRT <- Data$SO2.SC
M <- cbind(Data$SO2_Annual, Data$NOx_Annual, Data$CO2_Annual)
X <- cbind(Data$S_n_CR, Data$NumNOxControls, Data$Heat_Input/100000, Data$Barometric_Pressure, Data$Temperature,  Data$PctCapacity, Data$sulfur_Content, Data$Phase2_Indicator, Data$Operating_Time/1000)


dim.cov <- dim(X)[2] #<--------- Num. of Covariates

#------ Variables by treatments
x0 <- X[which(TRT==0),]
x1 <- X[which(TRT==1),]

y0 <- OUT[which(TRT==0)]
y1 <- OUT[which(TRT==1)]

m0 <- log(M[which(TRT==0),])
m1 <- log(M[which(TRT==1),])

n0 <- dim(x0)[1]
n1 <- dim(x1)[1]


#----- Extract MCMC samples from each marginal distribution
data.y0 <- para.y0[,2:35]
data.y1 <- para.y1[,2:35]
data.m10 <- para.m10[,2:35]
data.m11 <- para.m11[,2:35]
data.m20 <- para.m20[,2:35]
data.m21 <- para.m21[,2:35]
data.m30 <- para.m30[,2:35]
data.m31 <- para.m31[,2:35]

COR <- para.C[,29:56]

#----- Parameters of the models under Z=1
gamma10 <- data.y1[,9:17]
gamma11 <- data.y1[,18:25]
w1 <- data.y1[,26:34]
beta1_10 <- data.m11[,9:17]
beta1_11 <- data.m11[,18:25]
beta2_10 <- data.m21[,9:17]
beta2_11 <- data.m21[,18:25]
beta3_10 <- data.m31[,9:17]
beta3_11 <- data.m31[,18:25]
s1_1 <- data.m11[,26:34]
s2_1 <- data.m21[,26:34]
s3_1 <- data.m31[,26:34]
psi1 <- data.y1[,1:8]
psi1_1 <- data.m11[,1:8]
psi2_1 <- data.m21[,1:8]
psi3_1 <- data.m31[,1:8]


#----- Parameters of the models under Z=0
gamma00 <- data.y0[,9:17]
gamma01 <- data.y0[,18:25]
w0 <- data.y0[,26:34]
beta1_00 <- data.m10[,9:17]
beta1_01 <- data.m10[,18:25]
beta2_00 <- data.m20[,9:17]
beta2_01 <- data.m20[,18:25]
beta3_00 <- data.m30[,9:17]
beta3_01 <- data.m30[,18:25]
s1_0 <- data.m10[,26:34]
s2_0 <- data.m20[,26:34]
s3_0 <- data.m30[,26:34]
psi0 <- data.y0[,1:8]
psi1_0 <- data.m10[,1:8]
psi2_0 <- data.m20[,1:8]
psi3_0 <- data.m30[,1:8]

#----- Setting Thinning and Burn-in's
Thin <- 5
Burn0 <- 90000  # Extra burn-in periods for Models under Z=0
Burn1 <- 90000   # Extra burn-in periods for Models under Z=1


#----- Set the number of post-processing steps, N
n.iter <- 2000


##---------------------------------------------------------------
## 'Main' function on each cluster (parallel)
##---------------------------------------------------------------


  main <- function(temp){
      
      
      #----- Require multivariate normal distribuion on each cluster
      library(mnormt)
      
      #----- Index of iterations (on parallel)
      j <- temp
      
      #----- The number of covariates samples (with replacement) for the empirical distribution
      size <- 120000
      s.covariate <- data.frame(X[sample(seq(1,dim(X)[1]), size = size, replace = TRUE),])  #----- Covariates samples
      
      #----- Index of posterior samples
      index0 <- Thin * j + Burn0; index1 <- Thin * j + Burn1
      
      #----- Construct the correlation matrix
      COR1 <- matrix(c(1,COR[index0,c(1,2,3,5,6,7,1)],1,COR[index0,c(8,9,11,12,13,2,8)],1,COR[index0,c(14,16,17,18,3,9,14)],1,COR[index0,c(20,21,22,5,11,16,20)],1,COR[index0,c(26,27,6,12,17,21,26)],1,COR[index0,c(28,7,13,18,22,27,28)],1),7,7,byrow=TRUE)
      
      COR0 <- matrix(c(1,COR[index0,c(8,9,10,11,12,13,8)],1,COR[index0,c(14,15,16,17,18,9,14)],1,COR[index0,c(19,20,21,22,10,15,19)],1,COR[index0,c(23,24,25,11,16,20,23)],1,COR[index0,c(26,27,12,17,21,24,26)],1,COR[index0,c(28,13,18,22,25,27,28)],1),7,7,byrow=TRUE)
      
      #----- mixing parameters
      ppi.y0 <- pmax(psi0[index0,], 0)
      ppi.y1 <- pmax(psi1[index1,], 0)
      ppi.m10 <- pmax(psi1_0[index0,], 0)
      ppi.m20 <- pmax(psi2_0[index0,], 0)
      ppi.m30 <- pmax(psi3_0[index0,], 0)
      ppi.m11 <- pmax(psi1_1[index1,], 0)
      ppi.m21 <- pmax(psi2_1[index1,], 0)
      ppi.m31 <- pmax(psi3_1[index1,], 0)
      
      pi.y1 <- NULL
      pi.y1[1] <- ppi.y1[1]
      pi.y1[2:8] <- sapply(2:8, function(i) ppi.y1[i] * prod(1 - ppi.y1[1:(i-1)]))
      pi.y1[9] <- prod(1-ppi.y1[1:8])
      
      pi.y0 <- NULL
      pi.y0[1] <- ppi.y0[1]
      pi.y0[2:8] <- sapply(2:8, function(i) ppi.y0[i] * prod(1 - ppi.y0[1:(i-1)]))
      pi.y0[9] <- prod(1-ppi.y0[1:8])
      
      pi.m10 <- NULL
      pi.m10[1] <- ppi.m10[1]
      pi.m10[2:8] <- sapply(2:8, function(i) ppi.m10[i] * prod(1 - ppi.m10[1:(i-1)]))
      pi.m10[9] <- prod(1-ppi.m10[1:8])
      
      pi.m20 <- NULL
      pi.m20[1] <- ppi.m20[1]
      pi.m20[2:8] <- sapply(2:8, function(i) ppi.m20[i] * prod(1 - ppi.m20[1:(i-1)]))
      pi.m20[9] <- prod(1-ppi.m20[1:8])
      
      pi.m30 <- NULL
      pi.m30[1] <- ppi.m30[1]
      pi.m30[2:8] <- sapply(2:8, function(i) ppi.m30[i] * prod(1 - ppi.m30[1:(i-1)]))
      pi.m30[9] <- prod(1-ppi.m30[1:8])
      
      pi.m11 <- NULL
      pi.m11[1] <- ppi.m11[1]
      pi.m11[2:8] <- sapply(2:8, function(i) ppi.m11[i] * prod(1 - ppi.m11[1:(i-1)]))
      pi.m11[9] <- prod(1-ppi.m11[1:8])
      
      pi.m21 <- NULL
      pi.m21[1] <- ppi.m21[1]
      pi.m21[2:8] <- sapply(2:8, function(i) ppi.m21[i] * prod(1 - ppi.m21[1:(i-1)]))
      pi.m21[9] <- prod(1-ppi.m21[1:8])
      
      pi.m31 <- NULL
      pi.m31[1] <- ppi.m31[1]
      pi.m31[2:8] <- sapply(2:8, function(i) ppi.m31[i] * prod(1 - ppi.m31[1:(i-1)]))
      pi.m31[9] <- prod(1-ppi.m31[1:8])
      

      
      ##----- Sampling from the joint distribution of [Y(1;M(1,1,1)), M(0,0,0), M(1,1,1)]
      f1.M0M1 <- function(k){
          
          #----- Random samples from the Copula and take the standard normal CDF on them
          F <- pnorm(rmnorm(1, mean=rep(0, 7), COR1), 0, 1)
          X.temp <- s.covariate[k,]
          
          #----- Draw a cluster for each marginal distribution
          clus.y1 <- which(rmultinom(1, 1, pi.y1) == 1)
          clus.m10 <- which(rmultinom(1, 1, pi.m10) == 1)
          clus.m20 <- which(rmultinom(1, 1, pi.m20) == 1)
          clus.m30 <- which(rmultinom(1, 1, pi.m30) == 1)
          clus.m11 <- which(rmultinom(1, 1, pi.m11) == 1)
          clus.m21 <- which(rmultinom(1, 1, pi.m21) == 1)
          clus.m31 <- which(rmultinom(1, 1, pi.m31) == 1)
          
          #----- Samples from marginal distributions
          s.y1 <- qnorm(F[1], mean = gamma10[index1, clus.y1] + gamma11[index1, ] %*% t(X.temp), sd = sqrt(w1[index1, clus.y1]))
          s.m1_0 <- qnorm(F[5], mean = beta1_00[index0, clus.m10] + beta1_01[index0, ] %*% t(X.temp), sd = sqrt(s1_0[index0, clus.m10]))
          s.m2_0 <- qnorm(F[6], mean = beta2_00[index0, clus.m20] + beta2_01[index0, ] %*% t(X.temp), sd = sqrt(s2_0[index0, clus.m20]))
          s.m3_0 <- qnorm(F[7], mean = beta3_00[index0, clus.m30] + beta3_01[index0, ] %*% t(X.temp), sd = sqrt(s3_0[index0, clus.m30]))
          
          s.m1_1 <- qnorm(F[2], mean = beta1_10[index1, clus.m11] + beta1_11[index1, ] %*% t(X.temp), sd = sqrt(s1_1[index1, clus.m11]))
          s.m2_1 <- qnorm(F[3], mean = beta2_10[index1, clus.m21] + beta2_11[index1, ] %*% t(X.temp), sd = sqrt(s2_1[index1, clus.m21]))
          s.m3_1 <- qnorm(F[4], mean = beta3_10[index1, clus.m31] + beta3_11[index1, ] %*% t(X.temp), sd = sqrt(s3_1[index1, clus.m31]))
          return(c(s.y1, s.m1_0, s.m2_0, s.m3_0, s.m1_1, s.m2_1, s.m3_1))
      }
      
      ##----- Sampling from the joint distribution of [Y(0;M(0,0,0)), M(0,0,0), M(1,1,1)]
      f0.M0M1 <- function(k){
          
          #----- Random samples from the Copula and take the standard normal CDF on them
          F <- pnorm(rmnorm(1, mean = rep(0, 7), COR0), 0, 1)
          X.temp <- s.covariate[k, ]
          
          #----- Draw a cluster for each marginal distribution
          clus.y0 <- which(rmultinom(1, 1, pi.y0) == 1)
          clus.m10 <- which(rmultinom(1, 1, pi.m10) == 1)
          clus.m20 <- which(rmultinom(1, 1, pi.m20) == 1)
          clus.m30 <- which(rmultinom(1, 1, pi.m30) == 1)
          clus.m11 <- which(rmultinom(1, 1, pi.m11) == 1)
          clus.m21 <- which(rmultinom(1, 1, pi.m21) == 1)
          clus.m31 <- which(rmultinom(1, 1, pi.m31) == 1)
          
          #----- Samples from marginal distributions
          s.y0 <- qnorm(F[4], mean = gamma00[index0, clus.y0] + gamma01[index0, ] %*% t(X.temp), sd = sqrt(w0[index0, clus.y0]))
          s.m1_0 <- qnorm(F[5], mean = beta1_00[index0, clus.m10] + beta1_01[index0, ] %*% t(X.temp), sd = sqrt(s1_0[index0, clus.m10]))
          s.m2_0 <- qnorm(F[6], mean = beta2_00[index0, clus.m20] + beta2_01[index0, ] %*% t(X.temp), sd = sqrt(s2_0[index0, clus.m20]))
          s.m3_0 <- qnorm(F[7], mean = beta3_00[index0, clus.m30] + beta3_01[index0, ] %*% t(X.temp), sd = sqrt(s3_0[index0, clus.m30]))
          
          s.m1_1 <- qnorm(F[1], mean = beta1_10[index1, clus.m11] + beta1_11[index1, ] %*% t(X.temp), sd = sqrt(s1_1[index1, clus.m11]))
          s.m2_1 <- qnorm(F[2], mean = beta2_10[index1, clus.m21] + beta2_11[index1, ] %*% t(X.temp), sd = sqrt(s2_1[index1, clus.m21]))
          s.m3_1 <- qnorm(F[3], mean = beta3_10[index1, clus.m31] + beta3_11[index1, ] %*% t(X.temp), sd = sqrt(s3_1[index1, clus.m31]))
          return(c(s.y0, s.m1_0, s.m2_0, s.m3_0, s.m1_1, s.m2_1, s.m3_1))
      }
      
      
      #----- Estimating Y1 and Y0 conditional on all mediators
      y1.sample <- sapply(seq(1, dim(s.covariate)[1], by=1), function(x) f1.M0M1(x ))
      y0.sample <- sapply(seq(1, dim(s.covariate)[1], by=1), function(x) f0.M0M1(x ))
      


      #----- Estimating E[Y1-Y0] given M samples
      y1_1 <- apply(ind1, 1, function(x) mean(y1.sample[1,which(abs(y1.sample[2,]-x[1]) < 0.2 & abs(y1.sample[5,]-x[2]) < 0.2)]))
      y0_1 <- apply(ind1, 1, function(x) mean(y0.sample[1,which(abs(y0.sample[2,]-x[1]) < 0.2 & abs(y0.sample[5,]-x[2]) < 0.2)]))
      d_1 <- y1_1-y0_1

      y1_2 <- apply(ind2, 1, function(x) mean(y1.sample[1,which(abs(y1.sample[3,]-x[1]) < 0.2 & abs(y1.sample[6,]-x[2]) < 0.2)]))
      y0_2 <- apply(ind2, 1, function(x) mean(y0.sample[1,which(abs(y0.sample[3,]-x[1]) < 0.2 & abs(y0.sample[6,]-x[2]) < 0.2)]))
      d_2 <- y1_2-y0_2

      y1_3 <- apply(ind3, 1, function(x) mean(y1.sample[1,which(abs(y1.sample[4,]-x[1]) < 0.2 & abs(y1.sample[7,]-x[2]) < 0.2)]))
      y0_3 <- apply(ind3, 1, function(x) mean(y0.sample[1,which(abs(y0.sample[4,]-x[1]) < 0.2 & abs(y0.sample[7,]-x[2]) < 0.2)]))
      d_3 <- y1_3-y0_3

      return(c(d_1, d_2, d_3))

  }
  
  result<-foreach(temp = 1:n.iter, .combine = rbind) %dopar% main(temp)
  
  save(result, file="SurfaceFrame2.RData")
  stopCluster(cl)
  
