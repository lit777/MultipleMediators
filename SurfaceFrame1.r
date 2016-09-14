#----- Load the main dataset
load("Master.RData")

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



#----- Y samples from one iteration of the finite PS model
load("Ysample.RData")
C.sample <- t(result[1:258,3:8])


##---------------------------------------------------------------
## Preparing Plot1
##---------------------------------------------------------------

#----- Define grid of M1
x <- seq(min(log(M[,1])),max(log(M[,1])),length=30)[1:30]
y <- seq(min(log(M[,1])),max(log(M[,1])),length=30)[1:30]

#----- Find cells which include M1 samples
mat_1 <- matrix(0, ncol=30, nrow=30)
for(i in 1:258){
  for(j in 1:29){
    for(k in 1:29){
      if(C.sample[1,i] >= x[j] && C.sample[1,i] <= x[j+1] && C.sample[4,i] >= y[k] && C.sample[4,i] <= y[k+1]) mat_1[j,k] <- 1
    }
  }
}

#----- Find bounds of M1 samples on 2D surface
ind.min <- seq(1,30, by=1)
ind.max <- seq(1,30, by=1)
for(i in 1:30){
  if(sum(mat_1[,i])!=0){
    ind.min[i] <- min(which(mat_1[,i]==1));
    ind.max[i] <- max(which(mat_1[,i]==1));
    mat_1[max(1,ind.min[i]-1):min(ind.max[i]+1,30),i] <- 1}
}

ind.min<-seq(1,30, by=1)
ind.max<-seq(1,30, by=1)
for(i in 1:30){
  if(sum(mat_1[i,])!=0){
    ind.min[i] <- min(which(mat_1[i,]==1));
    ind.max[i] <- max(which(mat_1[i,]==1));
    mat_1[i,max(1,ind.min[i]-1):min(ind.max[i]+1,30)] <- 1}
}

pre.mat <- 0
for(i in 1:30){
  pre.mat <- pre.mat + length(which(mat_1[,i]==1))
}
ind1 <- matrix(ncol=2, nrow=pre.mat)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_1[i,j]==1){index<-1+index; ind1[index,1] <- x[i]; ind1[index,2] <- y[j]}
  }
}



##---------------------------------------------------------------
## Preparing Plot2
##---------------------------------------------------------------

#----- Define grid of M2
x <- seq(min(log(M[,2])),max(log(M[,2])),length=30)[1:30]
y <- seq(min(log(M[,2])),max(log(M[,2])),length=30)[1:30]

#----- Find cells which include M2 samples
mat_2 <- matrix(0, ncol=30, nrow=30)
for(i in 1:258){
  for(j in 1:29){
    for(k in 1:29){
      if(C.sample[2,i] >= x[j] && C.sample[2,i] <= x[j+1] && C.sample[5,i] >= y[k] && C.sample[5,i] <= y[k+1]) mat_2[j,k] <- 1
    }
  }
}

#----- Find bounds of M2 samples on 2D surface
ind.min <- seq(1,30, by=1)
ind.max <- seq(1,30, by=1)
for(i in 1:30){
  if(sum(mat_2[,i])!=0){
    ind.min[i] <- min(which(mat_2[,i]==1));
    ind.max[i] <- max(which(mat_2[,i]==1));
    mat_2[max(1,ind.min[i]-1):min(ind.max[i]+1,30),i] <- 1}
}

ind.min<-seq(1,30, by=1)
ind.max<-seq(1,30, by=1)
for(i in 1:30){
  if(sum(mat_2[i,])!=0){
    ind.min[i] <- min(which(mat_2[i,]==1));
    ind.max[i] <- max(which(mat_2[i,]==1));
    mat_2[i,max(1,ind.min[i]-1):min(ind.max[i]+1,30)] <- 1}
}

pre.mat <- 0
for(i in 1:30){
  pre.mat <- pre.mat + length(which(mat_2[,i]==1))
}
ind2 <- matrix(ncol=2, nrow=pre.mat)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_2[i,j]==1){index<-1+index; ind2[index,1] <- x[i]; ind2[index,2] <- y[j]}
  }
}


##---------------------------------------------------------------
## Preparing Plot3
##---------------------------------------------------------------

#----- Define grid of M1
x <- seq(min(log(M[,3])),max(log(M[,3])),length=30)[1:30]
y <- seq(min(log(M[,3])),max(log(M[,3])),length=30)[1:30]

#----- Find cells which include M1 samples
mat_3 <- matrix(0, ncol=30, nrow=30)
for(i in 1:258){
  for(j in 1:29){
    for(k in 1:29){
      if(C.sample[3,i] >= x[j] && C.sample[3,i] <= x[j+1] && C.sample[6,i] >= y[k] && C.sample[6,i] <= y[k+1]) mat_3[j,k] <- 1
    }
  }
}

#----- Find bounds of M1 samples on 2D surface
ind.min <- seq(1,30, by=1)
ind.max <- seq(1,30, by=1)
for(i in 1:30){
  if(sum(mat_3[,i])!=0){
    ind.min[i] <- min(which(mat_3[,i]==1));
    ind.max[i] <- max(which(mat_3[,i]==1));
    mat_3[max(1,ind.min[i]-1):min(ind.max[i]+1,30),i] <- 1}
}

ind.min<-seq(1,30, by=1)
ind.max<-seq(1,30, by=1)
for(i in 1:30){
  if(sum(mat_3[i,])!=0){
    ind.min[i] <- min(which(mat_3[i,]==1));
    ind.max[i] <- max(which(mat_3[i,]==1));
    mat_3[i,max(1,ind.min[i]-1):min(ind.max[i]+1,30)] <- 1}
}

pre.mat <- 0
for(i in 1:30){
  pre.mat <- pre.mat + length(which(mat_3[,i]==1))
}
ind3 <- matrix(ncol=2, nrow=pre.mat)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_3[i,j]==1){index<-1+index; ind3[index,1] <- x[i]; ind3[index,2] <- y[j]}
  }
}

save(ind1, ind2, ind3, mat_1, mat_2, mat_3, C.sample, file = "SurfaceFrame1.RData")
