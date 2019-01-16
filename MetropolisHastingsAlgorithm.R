


#------- log-likelihood function of parameters
loglik <- function(Y, X, w, beta0, beta, sigma,Z){
  lik <- sum(log(pmax(0.1^100,dtruncnorm(Y,a=0, b=Inf,beta0[Z]+ X%*%beta, sd=sqrt(sigma[Z])))))
  
  return(lik)
}

#------- log priors of hyper-parameters
logprior <- function(Y, X, w, beta0, beta, sigma, alpha,
                     mu_beta0, sigma_beta0, alpha_beta0, alpha_sigma, mu_pop){
  
  beta0_part <- sum(log(dnorm(beta0, mu_beta0, sqrt(sigma_beta0))))+
    log(dnorm(mu_beta0, mu_pop, sqrt(50)))+
    sum(log(dgamma(1/sigma_beta0,alpha_beta0, alpha_beta0*100)))+
    log(dunif(alpha_beta0,1,5)) # Intercept parameters in the outcome models
  beta_part <- log(dmnorm(beta, rep(0,dim(X)[2]), diag(100,dim(X)[2]))) # Coefficients (other than the intercept) parameters
  sigma_part <- sum(dgamma(sigma, 1, 1, log=TRUE)) # Variance parameters
  
  return(beta0_part + beta_part + sigma_part)
}

#------- log density of Gaussian kernel part of Copula model
logCopula <- function(h, index, h.cur, R){
  h[,index] <- h.cur # update h variable with the current values
  hR <- h%*%chol(R)
  return(0.5*sum(rowSums(h*h)-rowSums(hR*hR)) )
}



#------- metropolis-hastings Algorithm
metropolis <- function(h, X, Y, R, w_pre, beta0_pre, beta_pre, sigma_pre, alpha_pre, mu_beta0_pre,
                    sigma_beta0_pre, alpha_beta0_pre, alpha_sigma_pre, K, eps, del1, del2, del3,
                    index, cov1, cov2,cov3, zz,Z, gamma){
  # if zz=1, this is for the data under Z=1; if zz=0, this is for the data under Z=0   
    if(zz==1){
    YY <- Y[(1:n1)]; XX <- X[(1:n1),]; nn <- n1; nind <- 1:n1; 
  }else{
    YY <- Y[(n1+1):(n1+n0)]; XX <- X[(n1+1):(n1+n0),]; nn <- n0; nind <- (n1+1):(n1+n0)
  }
  
  H <- h
  
  piK_pre <- matrix(ncol=K, nrow=n0+n1)
  ww_pre <- w_pre
  for(k in 1:K){
      piK_pre[,k] <- dtruncnorm(Y, a=0, b=Inf,beta0_pre[k]+X%*%beta_pre, sqrt(sigma_pre[k]))*ww_pre[k]
  }
  piK <- sweep(piK_pre, 1, rowSums(piK_pre), "/")
  
  zZ <- NULL
  for(n in 1:(n0+n1)){
      piK[n,] <- ifelse(is.na(piK[n,]), 0, piK[n,])
      if(sum(piK[n,])==0){piK[n,] <- rep(1/K,K)}
      zZ[n] <- sample(1:K, size=1, prob=piK[n,], replace=TRUE)
  }


  # Estimate of the intercept
  mu_pop <- coef(lm(Y~cbind(c(rep(1,n1),rep(0,n0)), X)))[1]
  sigma_pop <- var(resid(lm(Y~cbind(c(rep(1,n1),rep(0,n0)), X))))
  
  # Define AR variables
  acc1 <- acc2 <- 1
  acc3 <- rep(1,dim(X)[2])
  acc4 <- acc5 <- acc6 <- rep(1,K)
  
  
  ###### Proposal distributions for hyper-parameters
  prop5 <- max(rgamma(1,del1*alpha_pre, del1),0.1^10) # alpha
  prop6 <- runif(1, 1, 5)
  prop7 <- alpha_sigma_pre # alpha_sigma
  prop8 <- rnorm(1,mean=mu_beta0_pre, sqrt(var(YY)/2)) # mu_beta0
  prop9 <- runif(1, sigma_beta0_pre-0.1, sigma_beta0_pre+0.1) # sigma_beta0
  
  rat <- logprior(Y, X, w_pre, beta0_pre, beta_pre, sigma_pre, alpha=prop5, mu_beta0=prop8, sigma_beta0=prop9, alpha_beta0=prop6, alpha_sigma=prop7, mu_pop)+
    dgamma(alpha_pre,del1*prop5, del1, log=TRUE)+
    dunif(alpha_beta0_pre, 1, 5, log=TRUE)+
    dnorm(mu_beta0_pre, prop8, sqrt(var(YY)/2),log=TRUE)+
    dunif(sigma_beta0_pre, prop9-0.1, prop9+0.1, log=TRUE)-
    logprior(Y, X, w_pre, beta0_pre, beta_pre, sigma_pre, alpha_pre, mu_beta0_pre, sigma_beta0_pre, alpha_beta0_pre, alpha_sigma_pre, mu_pop)-
    dgamma(prop5,del1*alpha_pre, del1, log=TRUE)-
    dunif(prop6,1,5, log=TRUE)-
    dnorm(prop8, mu_beta0_pre, sqrt(var(YY)/2),log=TRUE)-
    dunif(prop9, sigma_beta0_pre-0.1,sigma_beta0_pre+0.1, log=TRUE)
  

  if(log(runif(1))>rat | abs(rat)==Inf){
    prop5 <- alpha_pre
    prop6 <- alpha_beta0_pre
    prop7 <- alpha_sigma_pre
    prop8 <- mu_beta0_pre
    prop9 <- sigma_beta0_pre
    acc1=0
  }
  
  W <- rep(0,K-1)
  for(c in 1:(K-1)){
      W[c] <- pmin(pmax(rbeta(1, 1+length(which(zZ==c)), gamma+sum(length(which(zZ>c)))),0.00001),0.99999)
  }
  prop1=NULL
  prop1_temp=NULL
  prop1_temp[1]=1
  prop1[1]=W[1]
  for(j in 2:(K-1)){
      prop1_temp[j]=prop1_temp[j-1]*(1-W[j-1])
      prop1[j]=W[j]*prop1_temp[j]
  }
  prop1[K]=prop1_temp[K-1]*(1-W[K-1])
  

  ## update gamma (when gamma follows Gamma(1, 1))
  gamma <- rgamma(1, 1+K-1, 1-sum(log(1-W)))
  
    
    ###### Proposal distribution for the intercept parameter
    prop2 <- pprop2 <- beta0_pre # define variables
    
    sR <- solve(R)
    
    # Adaptive Metropolis-Within-Gibbs 
    # The variables are updated one at a time ('Optimal Proposal Distributions and Adaptive MCMC'; Rosenthal (2010))
    for(k in 1:K){
    pprop2[k] <- rnorm(1,mean=beta0_pre[k], 0.1)

    # Gaussian variable for the Copula evaluated at the outcome model with the current parameter values
    h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(Y, a=0, b=Inf, prop2[zZ]+X%*%beta_pre, sqrt(sigma_pre[zZ])))))
                                                     

    # Gaussian variable for the Copula evaluated at the outcome model with the proposed parameter values
    h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(Y, a=0, b=Inf,pprop2[zZ]+X%*%beta_pre, sqrt(sigma_pre[zZ])))))
    
    h.cur1 <- ifelse(is.na(h.cur1), qnorm(1-0.1^15), h.cur1)
    h.cur2 <- ifelse(is.na(h.cur2), qnorm(1-0.1^15), h.cur2)
    

    rat1 <- logCopula(H,index,h.cur2,sR) + loglik(Y,X,prop1,pprop2,beta_pre,sigma_pre,zZ) +
      logprior(Y,X,prop1,pprop2[k],beta_pre,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) +
      dnorm(prop2[k], pprop2[k], 0.1,log=TRUE)
    
    rat2 <- logCopula(H,index,h.cur1,sR) + loglik(Y,X,prop1,prop2,beta_pre,sigma_pre,zZ) +
      logprior(Y,X,prop1,prop2[k],beta_pre,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) +
      dnorm(pprop2[k], prop2[k], 0.1,log=TRUE)

    rat <- rat1 - rat2
    
    if(log(runif(1))>rat){
      pprop2[k] <- beta0_pre[k]
      H[,index] <- h.cur1
    }else{
      prop2[k] <- pprop2[k]
      H[,index] <- h.cur2
    }
    
    }

  ###### Proposal distributions for the regression coefficients (except the intercept)
  prop3 <- pprop3 <- beta_pre
  
  Cov2 = 2.38^2/((dim.cov)*10)*(cov2/2+diag(0.1^10,dim.cov))

    # Adaptive Metropolis-Within-Gibbs 
    # The variables are updated one at a time ('Optimal Proposal Distributionsand Adaptive MCMC'; Rosenthal (2010))  
    pprop3 <- rmnorm(1,mean=beta_pre, Cov2)
    
    
    # Gaussian variable for the Copula evaluated at the outcome model with the current parameter values
    h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(Y, a=0, b=Inf, prop2[zZ]+X%*%prop3, sqrt(sigma_pre[zZ])))))
    
    # Gaussian variable for the Copula evaluated at the outcome model with the proposed parameter values
    h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(Y, a=0, b=Inf, prop2[zZ]+X%*%pprop3, sqrt(sigma_pre[zZ])))))
    
    h.cur1 <- ifelse(is.na(h.cur1), qnorm(1-0.1^15), h.cur1)
    h.cur2 <- ifelse(is.na(h.cur2), qnorm(1-0.1^15), h.cur2)
    
    rat <- logCopula(H,index,h.cur2,sR) + loglik(Y,X,prop1,prop2,pprop3,sigma_pre,zZ) +
      logprior(Y,X,prop1,prop2,pprop3,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) +
      dmnorm(prop3, pprop3, Cov2,log=TRUE) -
      logCopula(H,index,h.cur1,sR) - loglik(Y,X,prop1,prop2,prop3,sigma_pre,zZ) -
      logprior(Y,X,prop1,prop2,prop3,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) -
      dmnorm(pprop3, prop3, Cov2,log=TRUE)


    if(log(runif(1))>rat){
      pprop3 <- beta_pre
      H[, index] <- h.cur1
    }else{
      prop3 <- pprop3
      H[, index] <- h.cur2
    }
  
  
  
  ###### Proposal distribution for the variance
    
  prop4 <- pprop4 <- sigma_pre
  
  for(k in 1:K){
  pprop4[k] <- rgamma(1, sigma_pre[k]*1000, 1000)
  

  # Gaussian variable for the Copula evaluated at the outcome model with the current parameter values
  h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(Y, a=0, b=Inf, prop2[zZ]+X%*%prop3, sqrt(prop4[zZ])))))
  

  # Gaussian variable for the Copula evaluated at the outcome model with the proposed parameter values
  h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15, ptruncnorm(Y, a=0, b=Inf, prop2[zZ]+X%*%prop3, sqrt(pprop4[zZ])))))
  
  h.cur1 <- ifelse(is.na(h.cur1), qnorm(1-0.1^15), h.cur1)
  h.cur2 <- ifelse(is.na(h.cur2), qnorm(1-0.1^15), h.cur2)
  
  rat1 <- logCopula(H,index,h.cur2,sR) + loglik(Y,X,prop1,prop2,prop3,pprop4,zZ) +
    sum(dgamma(pprop4[k], 1, 1, log=TRUE))+dgamma(prop4[k], pprop4[k]*1000,1000,log=TRUE)


  rat2 <- logCopula(H,index,h.cur1,sR) + loglik(Y,X,prop1,prop2,prop3,prop4,zZ) +
    sum(dgamma(prop4[k], 1, 1, log=TRUE))+dgamma(pprop4[k],prop4[k]*1000,1000,log=TRUE)

  rat <- rat1 - rat2  


  if(log(runif(1))>rat){
      pprop4[k] <- sigma_pre[k]
      H[,index] <- h.cur1
  }else{
      prop4[k] <- pprop4[k]
      H[,index] <- h.cur2
  }
}


  
return(c(prop1,prop2,prop3,prop4,prop5,prop6,prop7,prop8,prop9,acc2,acc3,acc4,acc5,acc6,gamma,zZ))
}

