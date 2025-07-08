library(rmutil)   # For generating from pareto distribution

## auxiliary function to generate reserve data and all for the Table
Setting_Gen_Tab1 <- function(Ka,
                              method = c("unif", "pareto", "gamma", "beta",
                                         "piecewise_unif"),
                              para1, para2){
  
  if(method == "unif"){
    method.in <- "unif"
    para1 <- 1
    para2 <- 20
    x.med.vec <- seq(0.1, 20, by= 0.995)
    r_k <- runif(Ka, 0.1, 3)
    reserve.price.Cutoff <- 1
  }else if( method == "piecewise_unif"){
    method.in <- "piecewise_unif"
    para1 <- 2
    para2 <- 4
    x.med.vec <- seq(0, 4, by = 4/20)
    r_k <- runif(Ka, 0.1, 1.5)
    reserve.price.Cutoff <- 1
  }else if (method == "pareto"){
    method.in <- "pareto"
    para1 <- 3
    para2 <- 100
    x.med.vec <- seq(0, 20, by= 1)
    r_k <- runif(Ka, 0.001, 0.1)
    reserve.price.Cutoff <- 0.05
  } else if( method == "gamma"){
    method.in <- "gamma"
    para1 <- 10
    para2 <- 2
    x.med.vec <- seq(0, 12, by= 12/20)
    r_k <- runif(Ka, 0.1, 3)
    reserve.price.Cutoff <- 1
  }else if (method == "beta"){  
    method.in <- "beta"
    para1 <- 2
    para2 <- 2
    x.med.vec <- seq(0, 1, by= 1/20)
    r_k <- runif(Ka, 0.001, 0.2)
    reserve.price.Cutoff <- 0.05
  }else{
    stop("The method is not supported. Please use one of the following methods: 
         'unif', 'pareto', 'gamma', 'beta', 'piecewise_unif'.")
  }
  
  return(list(method.in = method.in,
              para1 = para1,
              para2 = para2,
              x.med.vec = x.med.vec,
              r_k = r_k,
              reserve.price.Cutoff = reserve.price.Cutoff))
}

Setting_Gen_Tab1_F0 <- function(data,
                             method = c("unif", "pareto", "gamma", "beta",
                                        "piecewise_unif"),
                             para1, para2){
  ## True F:
  if(method == "unif"){
    trueF <- punif(data$Z_iT_i[,1], min = para1, max = para2)
  }else if( method == "piecewise_unif"){
    trueF <- sapply(data$Z_iT_i[,1],
                    function(y) {piecewise.unif.cdf(y,
                                                    para1 = para1,
                                                    para2 = para2
                    )
                    }
    )
  }else if (method == "pareto"){
    trueF <- ppareto(data$Z_iT_i[,1], m = para1, s = para2)
  } else if( method == "gamma"){
    trueF <- pgamma(data$Z_iT_i[,1], shape = para1, rate = para2)
  }else if (method == "beta"){  
    trueF <- pbeta(data$Z_iT_i[,1], shape1 = para1, shape2 = para2)
  }else{
    stop("The method is not supported. Please use one of the following methods: 
         'unif', 'pareto', 'gamma', 'beta', 'piecewise_unif'.")
  }
  
  return(list(trueF = trueF))
}

Setting_Gen_Tab1_F0_cdf <- function(method = c("unif", "pareto", "gamma", "beta",
                                           "piecewise_unif"),para1, para2){
  ## True F:
  if(method == "unif"){
    trueF <- function(y) punif(y, min = para1, max = para2)
  }else if( method == "piecewise_unif"){
    trueF <-function(y)  piecewise.unif.cdf(y,
                                            para1 = para1,
                                            para2 = para2)
  }else if (method == "pareto"){
    trueF <- function(y) ppareto(y, m = para1, s = para2)
  } else if( method == "gamma"){
    trueF <- function(y) pgamma(y, shape = para1, rate = para2)
  }else if (method == "beta"){  
    trueF <- function(y) pbeta(y, shape1 = para1, shape2 = para2)
  }else{
    stop("The method is not supported. Please use one of the following methods: 
         'unif', 'pareto', 'gamma', 'beta', 'piecewise_unif'.")
  }
  
  return(list(trueF = trueF))
}

Setting_Gen_Tab1_TV <- function(e1,data,trueF,
                                method = c("unif", "pareto", "gamma", "beta",
                                           "piecewise_unif"),
                                para1, para2){
  if(method == "unif"){
    res <- TotalVarDist(e1 = e1, e2 = Unif(Min = para1, Max = para2), asis.smooth.discretize = "smooth")
  }else if( method == "piecewise_unif"){
    res <- TotalVarDist(e1 = e1,
                     e2 = UnivarMixingDistribution(Unif(Min = 1, Max = para1), Unif(Min = para1+1, Max = para2)),
                     asis.smooth.discretize = "smooth")
  }else if (method == "pareto"){
    e2.pareto <- AbscontDistribution(p = linear_interpolation_cdf(x = data$Z_iT_i[,1], F.x = trueF))
    res <- TotalVarDist(e1 = e1, e2 = e2.pareto, asis.smooth.discretize = "smooth")
  } else if( method == "gamma"){
    res <- TotalVarDist(e1 = e1, e2 = Gammad(shape = para1, scale = 1/para2), asis.smooth.discretize = "smooth")
  }else if (method == "beta"){  
    res <- TotalVarDist(e1 = e1, e2 = Beta(shape1 = para1, shape2 = para2), asis.smooth.discretize = "smooth")
  }else{
    stop("The method is not supported. Please use one of the following methods: 
         'unif', 'pareto', 'gamma', 'beta', 'piecewise_unif'.")
  }
  
  return(list(TV = res))
}


#######################################################################################
## auxiliary to generate reserve data and all for Med1
Setting_Gen_Med1 <- function(Ka,
                              method = c("unif", "pareto", "gamma", "beta"),
                              para1, para2){
  
  if(method == "unif"){
    method.in <- "unif"
    para1 <- 1
    para2 <- 20
    y.seq <- seq(1, 20, by = 0.001) ## should start from 0
    F.0 <- function(x) punif(x, para1, para2) #(SOURAV)F.mle(r_k), where F.mle is the estimate of the MLE in the real data example and is F0 for simulation use case
    r_k <- F.0(runif(Ka, 0.1, 3))
    reserve.price.Cutoff <- F.0(1.1)
  }else if (method == "pareto"){
    method.in <- "pareto"
    para1 <- 3
    para2 <- 100
    y.seq <- seq(0.001, 25, by = 0.05)
    F.0 <- function(x) ppareto(x, para1, para2)
    r_k <- F.0(runif(Ka, 0.001, 0.1)) #(SOURAV)F.mle(r_k), where F.mle is the estimate of the MLE in the real data example and is F0 for simulation use case
    reserve.price.Cutoff <-  F.0(0.05)
  } else if( method == "gamma"){
    method.in <- "gamma"
    para1 <- 10
    para2 <- 2
    y.seq <- seq(0.001, 15, by = 0.05)
    F.0 <- function(x) pgamma(x, shape = para1, rate = para2)
    r_k <- F.0(runif(Ka, 0.1, 3))
    reserve.price.Cutoff <- F.0(1)
  }else if (method == "beta"){  
    method.in <- "beta"
    para1 <- 2
    para2 <- 2
    y.seq <- seq(0, 1, by= .01)
    F.0 <- function(x) pbeta(x, shape1 = para1, shape2 = para2)
    r_k <- F.0(runif(Ka, 0.001, 0.2))
    reserve.price.Cutoff <- F.0(0.05)
  }else{
    stop("The method is not supported. Please use one of the following methods: 
         'unif', 'pareto', 'gamma', 'beta'.")
  }
  
  return(list(method.in = method.in,
              para1 = para1,
              para2 = para2,
              y.seq = y.seq,
              F.0 = F.0,
              r_k = r_k,
              reserve.price.Cutoff = reserve.price.Cutoff))
}


