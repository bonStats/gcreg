#' WIP: Methods for manipulating oracle functions
#' 
#' scale_oracle(): Specify one of data or scaling function
#' 
#' @param f_oracle oracle function
#' @param x data
#' @param f_scale scaling function 
#' @param ... oracle function(s)
#' @return function that returns TRUE if point satisfies oracle
#' @export
#' @examples
#' #To do
#'

combine_oracles <- function(...) { # turn into combine.oracle? test if is oracle
  # to do
  oracles <- list(...)
  
  function(p, ...) {
    
    all(apply(oracles, p = p, ...))
      
    
  }
  
}