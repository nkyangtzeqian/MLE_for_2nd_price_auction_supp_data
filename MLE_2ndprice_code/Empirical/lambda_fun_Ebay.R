source("lambda_fun.R")
rm(list=c("h_list_gen"))

h_list_gen <- function(h_index){
  force(h_index)
  ## h_list: h and H functions for one stage
  ## h_index=1: constant; h_index=2: proportional ; h_index=3: 
  if(h_index==1){ return(list(function(t){1}, function(t){t})) }
  if(h_index==2){ return(list(
    function(t){1+3*t}, 
    function(t){t + 3*t^2/2}
  )) }
}
