##---------------------------------------------------------------
## Required libraries
##---------------------------------------------------------------

#----- Parallel computing
library(doParallel)

#----- Make clusters based on the number of CPU cores
cl<-makeCluster(4) 
registerDoParallel(cl)
getDoParWorkers()

#----- Support parallel excution
library(foreach)



##---------------------------------------------------------------
## Extract MCMC samples from Stan outputs
##---------------------------------------------------------------

#----- Load MCMC samples and Data
load("MCMC.RData")
load("SurfaceFrame1.RData")

##---------------------------------------------------------------
## 'Main' function on each cluster (parallel)
##---------------------------------------------------------------

main <- function(temp){
  i <- temp
  #  eval(parse(text=paste0("load('Surface_", i-1, ".RData')")))
  d_1 <- unlist(apply(ind1, 1, function(x) mean((Y_1111[[i]]-Y_0000[[i]])[which(abs(log(ME[[i]]$m0[,1]*10000)-(x[1])) < 0.4 & abs(log(ME[[i]]$m1[,1]*10000)-(x[2])) < 0.4)],na.rm=TRUE)))
  d_2 <- unlist(apply(ind2, 1, function(x) mean((Y_1111[[i]]-Y_0000[[i]])[which(abs(log(1000*ME[[i]]$m0[,2])-(x[1])) < 0.4 & abs(log(ME[[i]]$m1[,2]*1000)-(x[2])) < 0.4)],na.rm=TRUE)))
  d_3 <- unlist(apply(ind3, 1, function(x) mean((Y_1111[[i]]-Y_0000[[i]])[which(abs(log(ME[[i]]$m0[,3]*10000000)-x[1]) < 0.4 & abs(log(ME[[i]]$m1[,3]*10000000)-x[2]) < 0.4)],na.rm=TRUE)))
  
  return(c(d_1, d_2, d_3))
  
}

result<-foreach(temp = 501:2000, .combine = rbind) %dopar% main(temp)

save(result, file="SurfaceFrame2.RData")
stopCluster(cl)

