
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
      size <- 15000
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
      

    
      
      #----- Sampling for Y(1;M(0,0,0),M(1,1,1))
      f1.M0M1 <- function(m1, m2, m3, H){
          Cor01 <- COR1
          s.m1_0 <- qnorm(sum(apply(cbind(1:9,pi.m10), 1, function(x) x[2]*pnorm(m1, mean=beta1_00[index0,x[1]]+beta1_01[index0,]%*%H, sd=sqrt(s1_0[index0,x[1]])))),0,1)
          s.m2_0 <- qnorm(sum(apply(cbind(1:9,pi.m20), 1, function(x) x[2]*pnorm(m2, mean=beta2_00[index0,x[1]]+beta2_01[index0,]%*%H, sd=sqrt(s2_0[index0,x[1]])))),0,1)
          s.m3_0 <- qnorm(sum(apply(cbind(1:9,pi.m30), 1, function(x) x[2]*pnorm(m3, mean=beta3_00[index0,x[1]]+beta3_01[index0,]%*%H, sd=sqrt(s3_0[index0,x[1]])))),0,1)
          
          F <- pnorm(rmnorm(1, mean=Cor01[1:4,5:7]%*%solve(Cor01[5:7,5:7])%*%c(s.m1_0,s.m2_0,s.m3_0), Cor01[1:4,1:4]-Cor01[1:4,5:7]%*%solve(Cor01[5:7,5:7])%*%Cor01[5:7,1:4]),0,1)
          
          Fy1 <- function(hh,u){abs(sum(apply(cbind(1:9,pi.y1), 1, function(x) x[2]*pnorm(hh, mean=gamma10[index1,x[1]]+gamma11[index1,]%*%H, sd=sqrt(w1[index1,x[1]]))))-u)}
          Fm11 <- function(hh,u){abs(sum(apply(cbind(1:9,pi.m11), 1, function(x) x[2]*pnorm(hh, mean=beta1_10[index1,x[1]]+beta1_11[index1,]%*%H, sd=sqrt(s1_1[index1,x[1]]))))-u)}
          Fm21 <- function(hh,u){abs(sum(apply(cbind(1:9,pi.m21), 1, function(x) x[2]*pnorm(hh, mean=beta2_10[index1,x[1]]+beta2_11[index1,]%*%H, sd=sqrt(s2_1[index1,x[1]]))))-u)}
          Fm31 <- function(hh,u){abs(sum(apply(cbind(1:9,pi.m31), 1, function(x) x[2]*pnorm(hh, mean=beta3_10[index1,x[1]]+beta3_11[index1,]%*%H, sd=sqrt(s3_1[index1,x[1]]))))-u)}
        
          s.y1 <- optimize(Fy1, u=F[1], interval=c(2, 20))$minimum
          s.m1_1 <- optimize(Fm11,  u=F[2], interval=c(0, 10))$minimum
          s.m2_1 <- optimize(Fm21,  u=F[3], interval=c(0, 10))$minimum
          s.m3_1 <- optimize(Fm31,  u=F[4], interval=c(4, 15))$minimum
          
          return(c(s.y1,s.m1_1,s.m2_1,s.m3_1))
      }

      
      
      #----- Sampling for Y(0;M(0,0,0),M(1,1,1))
      f0.M0M1 <- function(m1, m2, m3, H){
          Cor01 <- COR0

          s.m1_1 <- qnorm(sum(apply(cbind(1:9,pi.m11), 1, function(x) x[2]*pnorm(m1, mean=beta1_10[index1,x[1]]+beta1_11[index1,]%*%H, sd=sqrt(s1_1[index1,x[1]])))),0,1)
          s.m2_1 <- qnorm(sum(apply(cbind(1:9,pi.m21), 1, function(x) x[2]*pnorm(m2, mean=beta2_10[index1,x[1]]+beta2_11[index1,]%*%H, sd=sqrt(s2_1[index1,x[1]])))),0,1)
          s.m3_1 <- qnorm(sum(apply(cbind(1:9,pi.m31), 1, function(x) x[2]*pnorm(m3, mean=beta3_10[index1,x[1]]+beta3_11[index1,]%*%H, sd=sqrt(s3_1[index1,x[1]])))),0,1)
          
          F <- pnorm(rmnorm(1, mean=Cor01[4:7,1:3]%*%solve(Cor01[1:3,1:3])%*%c(s.m1_1,s.m2_1,s.m3_1), Cor01[4:7,4:7]-Cor01[4:7,1:3]%*%solve(Cor01[1:3,1:3])%*%Cor01[1:3,4:7]),0,1)
          
          Fy0 <- function(hh,u){abs(sum(apply(cbind(1:9,pi.y0), 1, function(x) x[2]*pnorm(hh, mean=gamma00[index0,x[1]]+gamma01[index0,]%*%H, sd=sqrt(w0[index0,x[1]]))))-u)}
          Fm10 <- function(hh,u){abs(sum(apply(cbind(1:9,pi.m10), 1, function(x) x[2]*pnorm(hh, mean=beta1_00[index0,x[1]]+beta1_01[index0,]%*%H, sd=sqrt(s1_0[index0,x[1]]))))-u)}
          Fm20 <- function(hh,u){abs(sum(apply(cbind(1:9,pi.m20), 1, function(x) x[2]*pnorm(hh, mean=beta2_00[index0,x[1]]+beta2_01[index0,]%*%H, sd=sqrt(s2_0[index0,x[1]]))))-u)}
          Fm30 <- function(hh,u){abs(sum(apply(cbind(1:9,pi.m30), 1, function(x) x[2]*pnorm(hh, mean=beta3_00[index0,x[1]]+beta3_01[index0,]%*%H, sd=sqrt(s3_0[index0,x[1]]))))-u)}
          
          s.y0 <- optimize(Fy0, u=F[1], interval=c(2, 20))$minimum
          s.m1_0 <- optimize(Fm10,  u=F[2], interval=c(0, 10))$minimum
          s.m2_0 <- optimize(Fm20,  u=F[3], interval=c(0, 10))$minimum
          s.m3_0 <- optimize(Fm30,  u=F[4], interval=c(4, 15))$minimum
          return(c(s.y0,s.m1_0,s.m2_0,s.m3_0))
      }

      
      
      #----- Impute Y samples
      new.data <- cbind(TRT,OUT,log(M),X)
      new.data1 <- array(dim=c(dim(new.data)[1],4,40))
      for(j in 1:40){
          new.data1[,,j] <- t(apply(new.data, 1, function(x) x[1]*f0.M0M1(x[3],x[4],x[5],x[6:13])+(1-x[1])*f1.M0M1(x[3],x[4],x[5],x[6:13])))
      }
      data.list <- list(new.data1[,,1],new.data1[,,2],new.data1[,,3],new.data1[,,4],new.data1[,,5],new.data1[,,6],new.data1[,,7],new.data1[,,8],new.data1[,,9],new.data1[,,10],new.data1[,,11],new.data1[,,12],new.data1[,,13],new.data1[,,14],new.data1[,,15],new.data1[,,16],new.data1[,,17],new.data1[,,18],new.data1[,,19],new.data1[,,20],new.data1[,,21],new.data1[,,22],new.data1[,,23],new.data1[,,24],new.data1[,,25],new.data1[,,26],new.data1[,,27],new.data1[,,28],new.data1[,,29],new.data1[,,30],new.data1[,,31],new.data1[,,32],new.data1[,,33],new.data1[,,34],new.data1[,,35],new.data1[,,36],new.data1[,,37],new.data1[,,38],new.data1[,,39],new.data1[,,40])
      new.data2<-cbind(new.data, Reduce("+", data.list) / length(data.list))
      
      #----- Complete data
      y1.sample <- new.data2[which(new.data2[,1]==1),c(2,14,15,16,17,3,4,5)]
      y0.sample <- new.data2[which(new.data2[,1]==0),c(14,2,3,4,5,15,16,17)]
      y.sample <- rbind(y1.sample, y0.sample)
      
      
      
      return(y.sample)
      
  }
  
  result<-foreach(temp = 1:n.iter, .combine = rbind) %dopar% main(temp)
  
  save(result, file="Ysample.RData")
  stopCluster(cl)
  
