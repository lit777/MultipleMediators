#------- log-likelihood function of parameters
loglik <- function(Y, X, w, beta0, beta, sigma){
  dimx <- dim(X)[2]+1
  lik <- sum(apply(cbind(Y,X), 1, function(z) log(sum(w*dnorm(z[1],beta0+sum(beta*z[2:dimx]),sqrt(sigma))))))
  return(lik)
}

#------- log priors of hyper-parameters
logprior <- function(Y, X, w, beta0, beta, sigma, alpha, 
                     mu_beta0, sigma_beta0, alpha_beta0, alpha_sigma, mu_pop){
  
  w_part <- sum(log(dbeta(w,1,alpha))) + log(dgamma(alpha,0.1,0.1)) # theta_k^\prime parameters
  beta0_part <- sum(log(dnorm(beta0, mu_beta0, sqrt(sigma_beta0))))+ 
    log(dnorm(mu_beta0, mu_pop, sqrt(var(Y)/2)))+
    sum(log(dgamma(1/sigma_beta0,alpha_beta0, alpha_beta0*var(Y)/2)))+
    log(dunif(alpha_beta0,1,5)) # Intercept parameters in the outcome models
  beta_part <- log(dmnorm(beta, rep(0,dim(X)[2]), diag(100,dim(X)[2]))) # Coefficients (other than the intercept) parameters
  sigma_part <- sum(log(dgamma(1/sigma,alpha_sigma, alpha_sigma*var(Y)/2)))+log(dunif(alpha_sigma, 0.5, 10)) # Variance parameters
  
  return(w_part + beta0_part + beta_part + sigma_part)
}

#------- log density of Gaussian kernel part of Copula model
logCopula <- function(h, index, h.cur, R){
  h[1:(n0+n1),index] <- h.cur # update h variable with the current values
  return(0.5*sum(diag(h%*%(diag(P)-solve(R))%*%t(h))))
}

#------- metropolis-hastings Algorithm
metropolis <- function(h, X, Y, R, w_pre, beta0_pre, beta_pre, sigma_pre, alpha_pre, mu_beta0_pre,
                    sigma_beta0_pre, alpha_beta0_pre, alpha_sigma_pre, K, eps, del1, del2, del3,
                    index, cov1, cov2, zz){
  # if zz=1, this is for the data under Z=1; if zz=0, this is for the data under Z=0   
  if(zz==1){
    YY <- Y[(1:n1)]; XX <- X[(1:n1),]
  }else{
    YY <- Y[(n1+1):(n1+n0)]; XX <- X[(n1+1):(n1+n0),]
  }
  
  # Estimate of the intercept
  mu_pop <- coef(lm(YY~XX))[1]
  
  # Define AR variables
  acc1 <- acc2 <- 1
  acc3 <- rep(1,dim(X)[2])
  acc4 <- acc5 <- acc6 <- rep(1,K)
  
  
  ###### Proposal distributions for hyper-parameters
  prop5 <- max(rgamma(1,del1*alpha_pre, del1),0.1^10) # alpha
  prop6 <- max(rgamma(1,del2*alpha_beta0_pre, del2),0.1^10) # alpha_beta0
  prop7 <- max(rgamma(1,del3*alpha_sigma_pre, del3),0.1^10) # alpha_sigma
  prop8 <- rnorm(1,mean=mu_beta0_pre, sqrt(var(YY)/2)) # mu_beta0
  prop9 <- runif(1, 0, var(YY)) # sigma_beta0
  
  rat <- logprior(YY, X, w_pre, beta0_pre, beta_pre, sigma_pre, alpha=prop5, mu_beta0=prop8, sigma_beta0=prop9, alpha_beta0=prop6, alpha_sigma=prop7, mu_pop)+
    dgamma(alpha_pre,del1*prop5, del1, log=TRUE)+
    dgamma(alpha_beta0_pre, del2*prop6, del2, log=TRUE)+
    dnorm(mu_beta0_pre, prop8, sqrt(var(YY)/2),log=TRUE)+
    dgamma(alpha_sigma_pre, del3*prop7, del3, log=TRUE)+
    dunif(sigma_beta0_pre, 0, var(YY), log=TRUE)-
    logprior(YY, X, w_pre, beta0_pre, beta_pre, sigma_pre, alpha_pre, mu_beta0_pre, sigma_beta0_pre, alpha_beta0_pre, alpha_sigma_pre, mu_pop)-
    dgamma(prop5,del1*alpha_pre, del1, log=TRUE)-
    dgamma(prop6,del2*alpha_beta0_pre,del2, log=TRUE)-
    dgamma(prop7,del3*alpha_sigma_pre, del3, log=TRUE)-
    dnorm(prop8, mu_beta0_pre, sqrt(var(YY)/2),log=TRUE)-
    dunif(prop9, 0,var(YY), log=TRUE)
  
  if(log(runif(1))>rat){
    prop5 <- alpha_pre
    prop6 <- alpha_beta0_pre
    prop7 <- alpha_sigma_pre
    prop8 <- mu_beta0_pre
    prop9 <- sigma_beta0_pre
    acc1=0
  }
  
  
  
  ###### Proposal distribution for stick-breaking parameters ('theta_k^prime' parameter in the Appendix)
  prop1 <- runif((K-1), pmax(0,w_pre-eps), pmin(1,w_pre+eps)) 
  
  # Stick-breaking construction of the current values
  ww_pre <- NULL
  ww_pre[1] <- ppi.y1[1]
  ww_pre[2:(K-1)] <- sapply(2:(K-1), function(i) w_pre[i] * prod(1 - w_pre[1:(i-1)]))
  ww_pre[K] <- prod(1-w_pre[1:(K-1)])



  # Stick-breaking construction of the proposed values
  pprop1 <- NULL
  pprop1[1] <- prop1[1]
  pprop1[2:(K-1)] <- sapply(2:(K-1), function(i) prop1[i] * prod(1 - prop1[1:(i-1)]))
  pprop1[K] <- prod(1-prop1[1:(K-1)])


  # Gaussian variable for the Copula evaluated at the outcome model with the current parameter values
  h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(Y,X), 1, function(z) 
    sum(ww_pre*pnorm(z[1],beta0_pre+sum(beta_pre*z[2:dimx]),sqrt(sigma_pre)))))))
  # Gaussian variable for the Copula evaluated at the outcome model with the proposed parameter values
  h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(Y,X), 1, function(z) 
    sum(pprop1*pnorm(z[1],beta0_pre+sum(beta_pre*z[2:dimx]),sqrt(sigma_pre)))))))
  
  rat <- logCopula(h,index,h.cur2,R) + loglik(Y,X,pprop1,beta0_pre,beta_pre,sigma_pre) +
    logprior(YY,X,prop1,beta0_pre,beta_pre,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) +
    sum(dunif(w_pre, pmax(0,prop1-eps), pmin(1,prop1+eps), log=TRUE)) -
    logCopula(h,index,h.cur1,R) - loglik(Y,X,ww_pre,beta0_pre,beta_pre,sigma_pre) -
    logprior(YY,X,w_pre,beta0_pre,beta_pre,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) -
    sum(dunif(prop1, pmax(0,w_pre-eps), pmin(1,w_pre+eps), log=TRUE))
  
  if(log(runif(1))>rat){
    prop1 <- w_pre
    acc2 <- 0
    pprop1 <- ww_pre
  }
  
  
  ###### Proposal distribution for the intercept parameter
  prop2 <- pprop2 <- beta0_pre # define variables

  for(W in 1:K){

    # Adaptive Metropolis-Within-Gibbs 
    # The variables are updated one at a time ('Optimal Proposal Distributions and Adaptive MCMC'; Rosenthal (2010))
    pprop2[W] <- rnorm(1,mean=beta0_pre[W], sqrt(exp(2*cov1[W,W]))) 

    # Gaussian variable for the Copula evaluated at the outcome model with the current parameter values
    h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(Y,X), 1, function(z) 
      sum(pprop1*pnorm(z[1],prop2+sum(beta_pre*z[2:dimx]),sqrt(sigma_pre)))))))
    # Gaussian variable for the Copula evaluated at the outcome model with the proposed parameter values
    h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(Y,X), 1, function(z) 
      sum(pprop1*pnorm(z[1],pprop2+sum(beta_pre*z[2:dimx]),sqrt(sigma_pre)))))))
    
    rat <- logCopula(h,index,h.cur2,R) + loglik(Y,X,pprop1,pprop2,beta_pre,sigma_pre) +
      logprior(YY,X,prop1,pprop2,beta_pre,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) +
      dnorm(prop2[W], pprop2[W], sqrt(exp(2*cov1[W,W])),log=TRUE) -
      logCopula(h,index,h.cur1,R) - loglik(Y,X,pprop1,prop2,beta_pre,sigma_pre) -
      logprior(YY,X,prop1,prop2,beta_pre,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) -
      dnorm(pprop2[W], prop2[W], sqrt(exp(2*cov1[W,W])),log=TRUE)
    
    if(log(runif(1))>rat){
      pprop2[W] <- beta0_pre[W]
      acc4[W] <- 0
    }else{
      prop2[W] <- pprop2[W]
    }
  }
  
  
  ###### Proposal distributions for the regression coefficients (except the intercept)
  prop3 <- pprop3 <- beta_pre
  
  for(W in 1:dim.cov){

    # Adaptive Metropolis-Within-Gibbs 
    # The variables are updated one at a time ('Optimal Proposal Distributionsand Adaptive MCMC'; Rosenthal (2010))  
    pprop3[W] <- rnorm(1,mean=beta_pre[W], sqrt(exp(2*cov2[W,W]))) 
    
    # Gaussian variable for the Copula evaluated at the outcome model with the current parameter values
    h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(Y,X), 1, function(z) 
      sum(pprop1*pnorm(z[1],prop2+sum(prop3*z[2:dimx]),sqrt(sigma_pre)))))))
    # Gaussian variable for the Copula evaluated at the outcome model with the proposed parameter values
    h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(Y,X), 1, function(z) 
      sum(pprop1*pnorm(z[1],prop2+sum(pprop3*z[2:dimx]),sqrt(sigma_pre)))))))
    
    rat <- logCopula(h,index,h.cur2,R) + loglik(Y,X,pprop1,prop2,pprop3,sigma_pre) +
      logprior(YY,X,prop1,prop2,pprop3,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) +
      dnorm(prop3[W], pprop3[W], sqrt(exp(2*cov2[W,W])),log=TRUE) -
      logCopula(h,index,h.cur1,R) - loglik(Y,X,pprop1,prop2,prop3,sigma_pre) -
      logprior(YY,X,prop1,prop2,prop3,sigma_pre,alpha=prop5,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop) -
      dnorm(pprop3[W], prop3[W], sqrt(exp(2*cov2[W,W])),log=TRUE)

    if(log(runif(1))>rat){
      pprop3[W] <- beta_pre[W]
      acc6[W] <- 0
    }else{
      prop3[W] <- pprop3[W]
    }
  }
  
  
  ###### Proposal distribution for the variance     
  prop4 <- pprop4 <- sigma_pre
  for(W in 1:K){
    
    pprop4[W] <- runif(1, 0, var(YY))

    # Gaussian variable for the Copula evaluated at the outcome model with the current parameter values
    h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(Y,X), 1, function(z) 
      sum(pprop1*pnorm(z[1],prop2+sum(prop3*z[2:dimx]),sqrt(prop4)))))))
    # Gaussian variable for the Copula evaluated at the outcome model with the proposed parameter values
    h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,apply(cbind(Y,X), 1, function(z)
      sum(pprop1*pnorm(z[1],prop2+sum(prop3*z[2:dimx]),sqrt(pprop4)))))))
    
    rat <- logCopula(h,index,h.cur2,R) + loglik(Y,X,pprop1,prop2,prop3,pprop4) +
      sum(log(dgamma(1/pprop4,prop7, prop7*var(YY)/2)))+log(dunif(prop4[W], 0, var(YY))) -
      logCopula(h,index,h.cur1,R) - loglik(Y,X,pprop1,prop2,prop3,prop4) -
      sum(log(dgamma(1/prop4,prop7, prop7*var(YY)/2)))-log(dunif(pprop4[W], 0, var(YY)))

    if(log(runif(1))>rat){
      pprop4[W] <- sigma_pre[W]
      acc5[W] <- 0
    }else{
      prop4[W] <- pprop4[W]
    }
  }
  
  # return acceptance indicators and new values
  return(c(acc1,prop1,prop2,prop3,prop4,prop5,prop6,prop7,prop8,prop9,acc2,acc3,acc4,acc5,acc6))
}


