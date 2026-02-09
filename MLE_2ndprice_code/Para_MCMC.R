## Wrap up orginal function
MCMC_gamma <- function(ret, F.x){
  if(ret$class != "SecondPriceAuction_PolyaTree"){
    stop("This fucntion works only with data in the class 'SecondPriceAuction_PolyaTree.Rawdata'.
         See e.g., 'Data_Gen_Unobserved_PolyaTree' where we generate such data.")
  } 
  dat <- ret$dat
  Ka <- ret$Ka
  
  if (Ka==1000) {scale.gamma <- 0.02}
  if (Ka==100) {scale.gamma <- 0.05}
  gamma.init <- c(0.3,0.3);
  
  
  result.para <- EstCDF.para(dat, gamma.init, c(1,2), post.den.gamma, pfunc.gamma, dfunc.gamma, scale=scale.gamma)
  
  F.y <- sapply(F.x, function(x) result.para$cdf(x))
  return(list(MCMC_Mat=result.para$para.postsamp,
              fit_CDF=result.para$cdf,fit_PDF=result.para$den,
              F.x = F.x, F.y = F.y))
}

## Wrap up orginal function
MCMC_tNorm <- function(ret, F.x){
  if(ret$class != "SecondPriceAuction_PolyaTree"){
    stop("This fucntion works only with data in the class 'SecondPriceAuction_PolyaTree.Rawdata'.
         See e.g., 'Data_Gen_Unobserved_PolyaTree' where we generate such data.")
  } 
  dat <- ret$dat
  Ka <- ret$Ka
  
  if (Ka==1000) {scale.tnorm  <- 0.1}
  if (Ka==100) {scale.tnorm  <- 0.2}
  tnorm.init <- c(0,10);
  
  
  result.para <- EstCDF.para(dat, tnorm.init, 2, post.den.tnorm, pfunc.tnorm, dfunc.tnorm, scale=scale.tnorm)
  
  F.y <- sapply(F.x, function(x) result.para$cdf(x))
  return(list(MCMC_Mat=result.para$para.postsamp,
              fit_CDF=result.para$cdf,fit_PDF=result.para$den,
              F.x = F.x, F.y = F.y))
}




#################################################################################
## Original codes
#################################################################################


# Estimation of CDF using Parametric method
# -------------------------------------------------------------------------------
# dat: auction data
# init.para: Initial parameter for the parametric form
# pos.ind: +ve indicator
# post.den: posterior density for the parametric form
# 
EstCDF.para <- function(dat, init.para, pos.ind, post.den, pfunc, dfunc, nsim=2000, scale){
  k <- length(init.para)
  result <- list(para.postsamp=matrix(nrow=nsim, ncol=k), cdf=function(x){return(x)}, den=function(x){return(x)})
  para.postsamp <- MCMC(init.para=init.para, pos.ind=pos.ind, post.den=post.den, n=dat[,1], y=dat[,2], nsim=nsim, scale=scale)
  para.postsamp <- para.postsamp[-(1:floor(nsim/2)),]	# burn-in
  
  result$para.postsamp <- para.postsamp
  nmix <- nrow(result$para.postsamp)
  result$cdf <- function(x){ return(sum(pfunc(x, para.postsamp)) / nmix)}
  result$den <- function(x){ return(sum(dfunc(x, para.postsamp)) / nmix)}
  return(result)
}

# Generic Metropolis-Hasting algorithm (for parametric models) 
# ------------------------------------------------------------
MCMC <- function(init.para, pos.ind=NULL, post.den, n, y, nsim=2000, scale){
  acc <- 0
  postsamp <- matrix(nrow = nsim, ncol = length(init.para))
  para <- init.para
  LLold <- post.den(para, n, y)
  for (i in 1:nsim){
    # Draw theta and sigma2 conditional on n's and phi #
    paranew <- rnorm(length(init.para), para, scale)
    if (!is.null(pos.ind)){
      paranew[pos.ind] <- abs(paranew[pos.ind])
    }
    LLnew <- post.den(paranew, n, y)
    r <- min(1, exp(LLnew - LLold))
    if (rbinom(1,1,r)==1){
      LLold <- LLnew
      para <- paranew
      acc <- acc + 1
    }
    postsamp[i,] <- c(para)
    # if (i %% 100 == 0) { cat("Iteration ",i, " done\n");}
  }
  # cat(acc / nsim, " \n")
  return(postsamp)
}

# Functions for Gamma distribution
# ------------------------------------------
priorLL.para.gamma <- function(para){
  return(sum(dnorm(para,0,100, log=TRUE)))
}

LL1.gamma <- function(p1, p2, n1, y1){
  return(log(n1*(n1-1)) + log(1-pgamma(y1,p1,p2)) + dgamma(y1,p1,p2,log=TRUE) + (n1-2)*pgamma(y1,p1,p2,log.p=TRUE))
}

LL.gamma <- function(para, n, y){
  k <- length(y)
  p1 <- para[1]
  p2 <- para[2]
  loglik <- numeric(k)
  
  for (i in 1:k){
    loglik[i] <- LL1.gamma(p1, p2, n[i], y[i])
  }
  return(sum(loglik))
}

post.den.gamma <- function(para, n, y){
  return(priorLL.para.gamma(para) + LL.gamma(para, n, y))
}	

pfunc.gamma <- function(x, para.mat){
  return(pgamma(x, para.mat[,1], para.mat[,2]))
}

dfunc.gamma <- function(x, para.mat){
  return(dgamma(x, para.mat[,1], para.mat[,2]))
}

# Functions for Truncated-Normal distribution 
# ----------------------------------------
priorLL.para.tnorm <- function(para){
  return(sum(dnorm(para,0,100, log=TRUE)))
}

LL1.tnorm <- function(theta, sigma2, n1, y1){
  return(log(n1*(n1-1)) + log(1-ptnorm(y1,mean=theta,sd=sqrt(sigma2),lower=0)) + dtnorm(y1,mean=theta,sd=sqrt(sigma2),lower=0,log=TRUE) + (n1-2)*ptnorm(y1,mean=theta,sd=sqrt(sigma2),lower=0,log.p=TRUE))
}

LL.tnorm <- function(para, n, y){
  k <- length(y)
  theta <- para[1]
  sigma2 <- para[2]
  loglik <- numeric(k)
  
  for (i in 1:k){
    loglik[i] <- LL1.tnorm(theta, sigma2, n[i], y[i])
  }
  return(sum(loglik))
}

post.den.tnorm <- function(para, n, y){
  return(priorLL.para.tnorm(para) + LL.tnorm(para, n, y))
}	

pfunc.tnorm <- function(x, para.mat){
  return(ptnorm(rep(x, nrow(para.mat)), mean=para.mat[,1], sd = sqrt(para.mat[,2]), lower = 0))
}

dfunc.tnorm <- function(x, para.mat){
  return(dtnorm(rep(x, nrow(para.mat)), mean=para.mat[,1], sd = sqrt(para.mat[,2]), lower = 0))
}
