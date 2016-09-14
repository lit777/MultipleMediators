metropolisC = function(h,rho,prho){

  acc1 <- rep(NA,28)  # define AR variables
  prop1 <- rho  # define proposal variables 

  # Sampling for the identifiable parameters
  seq <- c(1,2,3,8,9,14,26,27,28)
  
  # Use an Independence Sampler for each element with an uniform distribution over positive definite constrained interval. 
  # Each interval is calculated as in Barnard, McCulloch and Meng (Statistica Sinica, 2000)
  for(q in seq){
    prop1[q] <- 1
    RHO1 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],prop1[4],prop1[10],
                         prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],prop1[6],prop1[12],prop1[17],
                         prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))
    prop1[q] <- 0
    RHO0 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],prop1[4],prop1[10],
                         prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],prop1[6],prop1[12],prop1[17],
                         prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))
    prop1[q] <- -1
    RHO_1 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],prop1[4],prop1[10],
                          prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],prop1[6],prop1[12],prop1[17],
                          prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))
    Fun <- function(r){(RHO1+RHO_1-2*RHO0)/2*r^2 + (RHO1-RHO_1)/2*r+RHO0}
    interval <- multiroot(Fun, start=c(-1.1,1.1), maxiter = 200)$root  # Interval satisfying positive definite matrix restriction

    # Propose a value from the interval (independence sampler)
    prop1[q] <- runif(1, min(interval),max(interval) )    
    
    # Proposed Correlation Matrix
    COR <- matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],prop1[4],prop1[10],
                 prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],prop1[6],prop1[12],prop1[17],prop1[21],
                 prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE)

    # Current Correlation matrix
    RHO <- matrix(c(1,rho[1:7],rho[1],1,rho[8:13],rho[2],rho[8],1,rho[14:18],rho[3],rho[9],rho[14],1,rho[19:22],rho[4],rho[10],rho[15],
                 rho[19],1,rho[1:3],rho[5],rho[11],rho[16],rho[20],rho[1],1,rho[26:27],rho[6],rho[12],rho[17],rho[21],rho[2],rho[26],1,
                 rho[28],rho[7],rho[13],rho[18],rho[22],rho[3],rho[27],rho[28],1),8,8,byrow=TRUE)

    # Accept (or reject) a proposed value
    if(is.positive.definite(COR)){        
        rat <- -(n1+n0)/2*log(det(COR))-
              0.5*sum(diag(h%*%(solve(COR))%*%t(h)))+
              dunif(prop1[q],min(interval),max(interval),log=TRUE)+
              dunif(rho[q],min(interval),max(interval),log=TRUE)+
              (n1+n0)/2*log(det(RHO))+0.5*sum(diag(h%*%(solve(RHO))%*%t(h)))-
              dunif(prop1[q],min(interval),max(interval),log=TRUE)-
              dunif(rho[q],min(interval),max(interval),log=TRUE)
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
  seq <- c(11,12,13,16,17,18,20,21,22); 
  
  # New parametrization with rho1 and finding the interval giving the positive definite matrix
  rrho[1] <- 1
  prop1[11] <- rrho[1]
  prop1[12] <- rrho[1]*(prop1[8]+prop1[26])/2
  prop1[13] <- rrho[1]*(prop1[9]+prop1[27])/2
  prop1[16] <- rrho[1]*(prop1[8]+prop1[26])/2
  prop1[17] <- rrho[1]
  prop1[18] <- rrho[1]*(prop1[14]+prop1[28])/2
  prop1[20] <- rrho[1]*(prop1[9]+prop1[27])/2
  prop1[21] <- rrho[1]*(prop1[14]+prop1[28])/2
  prop1[22] <- rrho[1]  
  RHO1 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],
                       prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],
                       prop1[6],prop1[12],prop1[17],prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],
                       prop1[27],prop1[28],1),8,8,byrow=TRUE))

  rrho[1] <- 0
  prop1[11] <- rrho[1]
  prop1[12] <- rrho[1]*(prop1[8]+prop1[26])/2
  prop1[13] <- rrho[1]*(prop1[9]+prop1[27])/2
  prop1[16] <- rrho[1]*(prop1[8]+prop1[26])/2
  prop1[17] <- rrho[1]
  prop1[18] <- rrho[1]*(prop1[14]+prop1[28])/2
  prop1[20] <- rrho[1]*(prop1[9]+prop1[27])/2
  prop1[21] <- rrho[1]*(prop1[14]+prop1[28])/2
  prop1[22] <- rrho[1]  
  RHO0 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],
                       prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],
                       prop1[6],prop1[12],prop1[17],prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],
                       prop1[27],prop1[28],1),8,8,byrow=TRUE))

  rrho[1] <- -1
  prop1[11] <- rrho[1]
  prop1[12] <- rrho[1]*(prop1[8]+prop1[26])/2
  prop1[13] <- rrho[1]*(prop1[9]+prop1[27])/2
  prop1[16] <- rrho[1]*(prop1[8]+prop1[26])/2
  prop1[17] <- rrho[1]
  prop1[18] <- rrho[1]*(prop1[14]+prop1[28])/2
  prop1[20] <- rrho[1]*(prop1[9]+prop1[27])/2
  prop1[21] <- rrho[1]*(prop1[14]+prop1[28])/2
  prop1[22] <- rrho[1]  
  RHO_1 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],
                        prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],
                        prop1[6],prop1[12],prop1[17],prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],
                        prop1[27],prop1[28],1),8,8,byrow=TRUE))
  
  Fun <- function(r){(RHO1+RHO_1-2*RHO0)/2*r^2 + (RHO1-RHO_1)/2*r+RHO0}
  interval <- multiroot(Fun, start=c(-1.1,1.1), maxiter = 200)$root  # Interval for the positive definite matrix
    
  # Propose a value from the interval (independence sampler)  
  rrho[1] <- runif(1, min(interval), max(interval) )
  
  # Proposed values for the original elements related to rho1
  prop1[11] <- rrho[1]
  prop1[12] <- rrho[1]*(prop1[8]+prop1[26])/2
  prop1[13] <- rrho[1]*(prop1[9]+prop1[27])/2
  prop1[16] <- rrho[1]*(prop1[8]+prop1[26])/2
  prop1[17] <- rrho[1]
  prop1[18] <- rrho[1]*(prop1[14]+prop1[28])/2
  prop1[20] <- rrho[1]*(prop1[9]+prop1[27])/2
  prop1[21] <- rrho[1]*(prop1[14]+prop1[28])/2
  prop1[22] <- rrho[1]  
  
  # Proposed Correlation Matrix
  COR <- matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,
               prop1[19:22],prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[23:25],prop1[5],prop1[11],prop1[16],prop1[20],
               prop1[23],1,prop1[26:27],prop1[6],prop1[12],prop1[17],prop1[21],prop1[24],prop1[26],1,prop1[28],prop1[7],
               prop1[13],prop1[18],prop1[22],prop1[25],prop1[27],prop1[28],1),8,8,byrow=TRUE)

  # Current Correlation Matrix
  RHO <- matrix(c(1,rho[1:7],rho[1],1,rho[8:13],rho[2],rho[8],1,rho[14:18],rho[3],rho[9],rho[14],1,rho[19:22],rho[4],
               rho[10],rho[15],rho[19],1,rho[23:25],rho[5],rho[11],rho[16],rho[20],rho[23],1,rho[26:27],rho[6],rho[12],
               rho[17],rho[21],rho[24],rho[26],1,rho[28],rho[7],rho[13],rho[18],rho[22],rho[25],rho[27],rho[28],1),8,8,byrow=TRUE)
  
  
  if(is.positive.definite(COR)){        
    rat <- -(n1+n0)/2*log(det(COR))-
      0.5*sum(diag(h%*%(solve(COR))%*%t(h)))+
      dunif(rrho[1],max(0,min(interval)),max(interval),log=TRUE)+
      dunif(prho[1],min(interval),max(interval),log=TRUE)+
      (n1+n0)/2*log(det(RHO))+0.5*sum(diag(h%*%(solve(RHO))%*%t(h)))-
      dunif(prho[1],max(0,min(interval)),max(interval),log=TRUE)-
      dunif(rrho[1],min(interval),max(interval),log=TRUE)

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
  
  
  seq <- c(5,6,7,10,15,19)
  
  # New parametrization with rho2 and finding the interval giving the positive definite matrix  
  rrho[2] <- 1
  prop1[5] <- rrho[2]*(prop1[1]+prop1[1])/2
  prop1[6] <- rrho[2]*(prop1[2]+prop1[2])/2
  prop1[7] <- rrho[2]*(prop1[3]+prop1[3])/2
  prop1[10] <- rrho[2]*(prop1[1]+prop1[1])/2
  prop1[15] <- rrho[2]*(prop1[2]+prop1[2])/2
  prop1[19] <- rrho[2]*(prop1[3]+prop1[3])/2
  RHO1 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],
                       prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],
                       prop1[6],prop1[12],prop1[17],prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))

  rrho[2] <- 0
  prop1[5] <- rrho[2]*(prop1[1]+prop1[1])/2
  prop1[6] <- rrho[2]*(prop1[2]+prop1[2])/2
  prop1[7] <- rrho[2]*(prop1[3]+prop1[3])/2
  prop1[10] <- rrho[2]*(prop1[1]+prop1[1])/2
  prop1[15] <- rrho[2]*(prop1[2]+prop1[2])/2
  prop1[19] <- rrho[2]*(prop1[3]+prop1[3])/2
  RHO0 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],
                       prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],
                       prop1[6],prop1[12],prop1[17],prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))

  rrho[2] <- -1
  prop1[5] <- rrho[2]*(prop1[1]+prop1[1])/2
  prop1[6] <- rrho[2]*(prop1[2]+prop1[2])/2
  prop1[7] <- rrho[2]*(prop1[3]+prop1[3])/2
  prop1[10] <- rrho[2]*(prop1[1]+prop1[1])/2
  prop1[15] <- rrho[2]*(prop1[2]+prop1[2])/2
  prop1[19] <- rrho[2]*(prop1[3]+prop1[3])/2
  RHO_1 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],
                        prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],
                        prop1[6],prop1[12],prop1[17],prop1[21],prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))
  
  Fun <- function(r){(RHO1+RHO_1-2*RHO0)/2*r^2 + (RHO1-RHO_1)/2*r+RHO0}
  interval <- multiroot(Fun, start=c(-1.1,1.1), maxiter = 200)$root # Interval for the positive definite matrix
  
  # Propose a value from the interval (independence sampler)  
  rrho[2] <- runif(1, min(interval),max(interval))
  
  # Proposed values for the original elements related to rho1
  prop1[5] <- rrho[2]*(prop1[1]+prop1[1])/2
  prop1[6] <- rrho[2]*(prop1[2]+prop1[2])/2
  prop1[7] <- rrho[2]*(prop1[3]+prop1[3])/2
  prop1[10] <- rrho[2]*(prop1[1]+prop1[1])/2
  prop1[15] <- rrho[2]*(prop1[2]+prop1[2])/2
  prop1[19] <- rrho[2]*(prop1[3]+prop1[3])/2
  
  # Proposed Correlation Matrix
  COR <- matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],
                  prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[23:25],prop1[5],prop1[11],prop1[16],prop1[20],prop1[23],1,prop1[26:27],
                  prop1[6],prop1[12],prop1[17],prop1[21],prop1[24],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[25],prop1[27],prop1[28],1),8,8,byrow=TRUE)
  # Current Correlation Matrix
  RHO <- matrix(c(1,rho[1:7],rho[1],1,rho[8:13],rho[2],rho[8],1,rho[14:18],rho[3],rho[9],rho[14],1,rho[19:22],rho[4],rho[10],rho[15],
                  rho[19],1,rho[23:25],rho[5],rho[11],rho[16],rho[20],rho[23],1,rho[26:27],rho[6],rho[12],rho[17],rho[21],rho[24],rho[26],1,
                  rho[28],rho[7],rho[13],rho[18],rho[22],rho[25],rho[27],rho[28],1),8,8,byrow=TRUE)
  
  
  if(is.positive.definite(COR)){        
    rat= -(n1+n0)/2*log(det(COR))-0.5*sum(diag(h%*%(solve(COR))%*%t(h)))+
      dunif(rrho[2],max(0,min(interval)),max(interval),log=TRUE)+
      dunif(prho[2],min(interval),max(interval),log=TRUE)+
      (n1+n0)/2*log(det(RHO))+0.5*sum(diag(h%*%(solve(RHO))%*%t(h)))-
      dunif(prho[2],max(0,min(interval)),max(interval),log=TRUE)-
      dunif(rrho[2],min(interval),max(interval),log=TRUE)
    
    if(is.na(rat)){
      prop1[seq] <- rho[seq]
      acc1[seq] <- 0
      rrho[2] <- prho[2]
    }else{
      if(log(runif(1))>rat){
        prop1[seq] <- rho[seq]
        acc1[seq] <- 0
        rrho[2] <- prho[2]
      }else{rho[seq] <- prop1[seq]
            }
      }
    }else{
      prop1[seq] <- rho[seq]
      acc1[seq] <- 0
      rrho[2] <- prho[2]
    }
  
  
  # New parametrization with rho3 and finding the interval giving the positive definite matrix    
  seq <- 4
  rrho[3] <- 1
  prop1[4] <- rrho[3]
  RHO1 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],
                       prop1[14],1,prop1[19:22],prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],
                       prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],prop1[6],prop1[12],prop1[17],prop1[21],
                       prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))

  rrho[3] <- 0
  prop1[4] <- rrho[3]
  RHO0 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],
                       prop1[14],1,prop1[19:22],prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],
                       prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],prop1[6],prop1[12],prop1[17],prop1[21],
                       prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))

  rrho[3] <- -1
  prop1[4] <- rrho[3]
  RHO_1 <- det(matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],
                        prop1[14],1,prop1[19:22],prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[1:3],prop1[5],
                        prop1[11],prop1[16],prop1[20],prop1[1],1,prop1[26:27],prop1[6],prop1[12],prop1[17],prop1[21],
                        prop1[2],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[3],prop1[27],prop1[28],1),8,8,byrow=TRUE))
  
  Fun <- function(r){(RHO1+RHO_1-2*RHO0)/2*r^2 + (RHO1-RHO_1)/2*r+RHO0}
  interval <- multiroot(Fun, start=c(-1.1,1.1), maxiter = 200)$root # Interval for the positive definite matrix

  # Propose a value from the interval (independence sampler)  
  rrho[3] <- runif(1, min(interval),max(interval) )
  
  # Proposed value for the original element related to rho1
  prop1[4] <- rrho[3]
  
  # Proposed Correlation Matrix
  COR <- matrix(c(1,prop1[1:7],prop1[1],1,prop1[8:13],prop1[2],prop1[8],1,prop1[14:18],prop1[3],prop1[9],prop1[14],1,prop1[19:22],
                  prop1[4],prop1[10],prop1[15],prop1[19],1,prop1[23:25],prop1[5],prop1[11],prop1[16],prop1[20],prop1[23],1,prop1[26:27],
                  prop1[6],prop1[12],prop1[17],prop1[21],prop1[24],prop1[26],1,prop1[28],prop1[7],prop1[13],prop1[18],prop1[22],prop1[25],prop1[27],prop1[28],1),8,8,byrow=TRUE)

  # Current Correlation Matrix
  RHO=matrix(c(1,rho[1:7],rho[1],1,rho[8:13],rho[2],rho[8],1,rho[14:18],rho[3],rho[9],rho[14],1,rho[19:22],rho[4],rho[10],rho[15],
               rho[19],1,rho[23:25],rho[5],rho[11],rho[16],rho[20],rho[23],1,rho[26:27],rho[6],rho[12],rho[17],rho[21],rho[24],rho[26],1,
               rho[28],rho[7],rho[13],rho[18],rho[22],rho[25],rho[27],rho[28],1),8,8,byrow=TRUE)
  
  if(is.positive.definite(COR)){        
    rat= -(n1+n0)/2*log(det(COR))-0.5*sum(diag(h%*%(solve(COR))%*%t(h)))+dunif(rrho[3],max(0,min(interval)),max(interval),log=TRUE)+dunif(prho[3],min(interval),max(interval),log=TRUE)+(n1+n0)/2*log(det(RHO))+0.5*sum(diag(h%*%(solve(RHO))%*%t(h)))-dunif(prho[3],max(0,min(interval)),max(interval),log=TRUE)-dunif(rrho[3],min(interval),max(interval),log=TRUE)

    if(is.na(rat)){
      prop1[seq] <- rho[seq]
      acc1[seq] <- 0
      rrho[3] <- prho[3]
    }else{
      if(log(runif(1))>rat){
        prop1[seq] <- rho[seq]
        acc1[seq] <- 0
        rrho[3] <- prho[3]
      }else{
        rho[seq] <- prop1[seq]
        }
      }
    }else{
      prop1[seq] <- rho[seq]
      acc1[seq] <- 0
      rrho[3] <- prho[3]
    }
  
  # As an additional assumption to sharpen posterior inference, 
  # we assume that Cor(M_k(1),Y(1;M(1,1,1))) = Cor(M_k(0),Y(0;M(0,0,0))) for all k = 1,2,3.
  prop1[23:25] <- prop1[1:3]

  # Return AR and accepted values and 3 new parameters
  return(c(acc1,prop1,rrho))
}


