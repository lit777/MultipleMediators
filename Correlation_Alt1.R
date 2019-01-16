metropolisC = function(h,rho,prho){
  
  acc1 <- rep(NA,15)  # define AR variables
  prop1 <- rho  # define proposal variables 
  
  # Sampling for the identifiable parameters
  seq <- c(1,2,6,13,14,15)
    
  mem <- function(prop1){
    mempty<-matrix(c(1,prop1[1:5],prop1[1],1,prop1[6:9],prop1[2],prop1[6],1,prop1[10:12],prop1[3],
                     prop1[7],prop1[10],1,prop1[13:14],prop1[4],prop1[8],prop1[11],prop1[13],1,
                     prop1[15],prop1[5],prop1[9],prop1[12],prop1[14],prop1[15],1),6,6,byrow=TRUE)
    return(mempty)
  }
  
  # Use an Independence Sampler for each element with an uniform distribution over positive definite constrained interval. 
  # Each interval is calculated as in Barnard, McCulloch and Meng (Statistica Sinica, 2000)
  for(q in seq){
    prop1[q] <- 1
    RHO1 <- det(mem(prop1))
    prop1[q] <- 0
    RHO0 <- det(mem(prop1))
    prop1[q] <- -1
    RHO_1 <- det(mem(prop1))
    Fun <- function(r){(RHO1+RHO_1-2*RHO0)/2*r^2 + (RHO1-RHO_1)/2*r+RHO0}
    interval <- multiroot(Fun, start=c(-1.1,1.1), maxiter = 200)$root  # Interval satisfying positive definite matrix restriction
    
    # Propose a value from the interval (independence sampler)
    prop1[q] <- runif(1, rho[q]-0.3,rho[q]+0.3 )      
    
    # Proposed Correlation Matrix
    COR <-  mem(prop1)
    
    # Current Correlation matrix
    RHO <-  mem(rho)
    
    # Accept (or reject) a proposed value
    if(is.positive.definite(COR)){        
      rat <- -(n1+n0)/2*log(det(COR))-
        0.5*sum(diag(h%*%(solve(COR))%*%t(h)))+
        dunif(prop1[q],min(interval),max(interval),log=TRUE)+
        dunif(rho[q],prop1[q]-0.3,prop1[q]+0.3,log=TRUE)+
        (n1+n0)/2*log(det(RHO))+0.5*sum(diag(h%*%(solve(RHO))%*%t(h)))-
        dunif(rho[q],min(interval),max(interval),log=TRUE)-
        dunif(prop1[q],rho[q]-0.3,rho[q]+0.3,log=TRUE)
      if(is.na(rat)){
        prop1[q] <- rho[q]
        acc1[q] <- 0
      }else{
        if(log(runif(1))>rat) {
          prop1[q] <- rho[q]
          acc1[q] <- 0
        }else{
          rho[q] <- prop1[q]}
      }
    }else{
      prop1[q] <- rho[q]
      acc1[q] <- 0
    }    
  }
  
  rrho <- NULL # define 3 new parameters for an alternative specification of (partially-identified) parameters  
  seq <- c(3,4,5,7,8,9,10,11,12); 
  
  # New parametrization with rho1 and finding the interval giving the positive definite matrix
  rrho[1] <- 1
  prop1[3] <- rrho[1]
  prop1[4] <- rrho[1]*(prop1[1]+prop1[13])/2
  prop1[5] <- rrho[1]*(prop1[2]+prop1[14])/2
  prop1[7] <- rrho[1]*(prop1[1]+prop1[13])/2
  prop1[8] <- rrho[1]
  prop1[9] <- rrho[1]*(prop1[6]+prop1[15])/2
  prop1[10] <- rrho[1]*(prop1[2]+prop1[14])/2
  prop1[11] <- rrho[1]*(prop1[6]+prop1[15])/2
  prop1[12] <- rrho[1]  
  RHO1 <- det(mem(prop1))
  
  rrho[1] <- 0
  prop1[3] <- rrho[1]
  prop1[4] <- rrho[1]*(prop1[1]+prop1[13])/2
  prop1[5] <- rrho[1]*(prop1[2]+prop1[14])/2
  prop1[7] <- rrho[1]*(prop1[1]+prop1[13])/2
  prop1[8] <- rrho[1]
  prop1[9] <- rrho[1]*(prop1[6]+prop1[15])/2
  prop1[10] <- rrho[1]*(prop1[2]+prop1[14])/2
  prop1[11] <- rrho[1]*(prop1[6]+prop1[15])/2
  prop1[12] <- rrho[1]  
  RHO0 <- det(mem(prop1))
  
  rrho[1] <- -1
  prop1[3] <- rrho[1]
  prop1[4] <- rrho[1]*(prop1[1]+prop1[13])/2
  prop1[5] <- rrho[1]*(prop1[2]+prop1[14])/2
  prop1[7] <- rrho[1]*(prop1[1]+prop1[13])/2
  prop1[8] <- rrho[1]
  prop1[9] <- rrho[1]*(prop1[6]+prop1[15])/2
  prop1[10] <- rrho[1]*(prop1[2]+prop1[14])/2
  prop1[11] <- rrho[1]*(prop1[6]+prop1[15])/2
  prop1[12] <- rrho[1]  
  RHO_1 <- det(mem(prop1))
  
  Fun <- function(r){(RHO1+RHO_1-2*RHO0)/2*r^2 + (RHO1-RHO_1)/2*r+RHO0}
  interval <- multiroot(Fun, start=c(-1.1,1.1), maxiter = 200)$root  # Interval for the positive definite matrix
  
  # Propose a value from the interval (independence sampler)  
  rrho[1] <- runif(1, prho[1]-0.3,prho[1]+0.3 )    
  
  # Proposed values for the original elements related to rho1
  prop1[3] <- rrho[1]
  prop1[4] <- rrho[1]*(prop1[1]+prop1[13])/2
  prop1[5] <- rrho[1]*(prop1[2]+prop1[14])/2
  prop1[7] <- rrho[1]*(prop1[1]+prop1[13])/2
  prop1[8] <- rrho[1]
  prop1[9] <- rrho[1]*(prop1[6]+prop1[15])/2
  prop1[10] <- rrho[1]*(prop1[2]+prop1[14])/2
  prop1[11] <- rrho[1]*(prop1[6]+prop1[15])/2
  prop1[12] <- rrho[1]
  
  # Proposed Correlation Matrix
  COR <- mem(prop1)
  
  # Current Correlation Matrix
  RHO <- mem(rho)
  
  
  if(is.positive.definite(COR)){        
    rat <- -(n1+n0)/2*log(det(COR))-
      0.5*sum(diag(h%*%(solve(COR))%*%t(h)))+
      dunif(rrho[1],max(0,min(interval)),max(interval),log=TRUE)+
      dunif(prho[1],rrho[1]-0.3,rrho[1]+0.3,log=TRUE)+
      (n1+n0)/2*log(det(RHO))+0.5*sum(diag(h%*%(solve(RHO))%*%t(h)))-
      dunif(prho[1],max(0,min(interval)),max(interval),log=TRUE)-
      dunif(rrho[1],prho[1]-0.3,prho[1]+0.3,log=TRUE)
    
    if(is.na(rat)){
      prop1[seq] <- rho[seq]
      acc1[seq] <- 0
      rrho[1] <- prho[1]
    }else{
      if(log(runif(1))>rat){
        prop1[seq] <- rho[seq]
        acc1[seq] <- 0
        rrho[1] <- prho[1]
      }else{rho[seq] <- prop1[seq]
      }
    }
  }else{
    prop1[seq] <- rho[seq]
    acc1[seq] <- 0
    rrho[1] <- prho[1]
  }
  
 
  # Return AR and accepted values and 3 new parameters
  return(c(acc1,prop1,rrho))
}

