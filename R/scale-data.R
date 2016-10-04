gen_scale_data_funcs <- function(x){
  #scale to min = -1, max = 1
  xmax <- max(x) 
  b <- 2/(xmax - min(x))
  c <- 1 - b * xmax
  
  return(
    list(scale   = function(x){ b * x + c },
         unscale = function(x){ (x - c)/b }
    )
  )
  
}