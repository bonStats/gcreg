#' WIP: Combine oracle functions to be used in contrained estimation procedure
#' 
#' @param ... oracle functions to be combined
#' @return function that returns TRUE if point satisfies oracle
#' @export
#' @examples
#' #To do
#'

combine_oracles <- function(...) { # turn into combine.oracle? test if is oracle
  
  oracles <- list(...)
  
  function(p, ...) {
    
    all(apply(oracles, p = p, ...))
      
    
  }
  
}