cost_def <- function(method = c("unif", "pareto", "gamma", "beta",
                                "piecewise_unif"),
                     para1, para2){
  ## set cost to be 0.25th quantile of the distribution
  if(method == "unif"){
    max_cutpt <- para2
    cost <- 5.75
  }else if( method == "piecewise_unif"){
    max_cutpt <- para2
    cost <- 1.5
  }else if (method == "pareto"){
    max_cutpt <- 22
    cost <- 0.856
  }else if( method == "gamma"){
    max_cutpt <- 12.5
    cost <- 3.863
  }else if (method == "beta"){  
    max_cutpt <- 1
    cost <- 0.326
  }else{
    stop("The method is not supported. Please use one of the following methods: 
         'unif', 'pareto', 'gamma', 'beta', 'piecewise_unif'.")
  }
  return(c(cost=cost, max_cutpt=max_cutpt))
}

ProfitMetric <- function(e1,                              
                         method = c("unif", "pareto", "gamma", "beta",
                                    "piecewise_unif"),
                         para1, para2){
  
  ## codes from George and Hui (2012) to compute Revenue and price
  Compute.Opt.Price <- function(cdf, cost=5, max.price=20){
    result <- list(opt.price=numeric(1), profit=numeric(1))
    tempgrid <- seq(cost, max.price, 0.001)
    val <- sapply(tempgrid, Profit.Func, cdf=cdf, cost=cost)
    result$opt.price <- tempgrid[which.max(val)]
    result$profit <- max(val)
    return(result)
  }
  
  Profit.Func <- function(x, cdf, cost){
    return((1 - cdf(x)) * (x - cost))
  }
  
  

  max_cutpt <- cost_def(method, para1, para2)["max_cutpt"]
  cost <- cost_def(method, para1, para2)["cost"]
  
  e1_cdf <- p(e1)
  
  
  Profit_all <- Compute.Opt.Price(e1_cdf, cost, max.price=max_cutpt)
  max_price <- Profit_all$opt.price
  est_Profit <- Profit_all$profit
  Profit <- Profit.Func(max_price, Setting_Gen_Tab1_F0_cdf(method,para1, para2)[[1]], cost)
  
  return(list(max_price = max_price,
              Profit = Profit,
              est_Profit = est_Profit))
}


