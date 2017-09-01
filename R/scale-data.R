#' Create function to scale and unscale data.
#' 
#' Scale x variables to be between -1 and 1 which improves stability of algorithms working with polynomials
#' @param x vector of data.
#' @export
#' @examples
#' x <- rnorm(100, mean = 3)
#' sc <- gen_scale_data_funs(x)
#' x_scaled <- sc$scale(x)
#' all(abs(sc$unscale(x_scaled) - x) < 1e-10)
#' 
gen_scale_data_funs <- function(x){
  #scale to min = -1, max = 1
  xmax <- max(x) 
  b <- 2/(xmax - min(x))
  c <- 1 - b * xmax
  
  fs <- list(scale   = function(x){ b * x + c },
             unscale = function(x){ (x - c)/b }
  )
  
  attr(fs, "b") <- b
  attr(fs, "c") <- c
  
  return(fs)
  
}
