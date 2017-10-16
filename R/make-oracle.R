#' WIP: Methods for manipulating oracle functions
#' 
#' 
#' @param oracle oracle function, returns TRUE if point is in constrained region, FALSE otherwise.
#' @param region region over which constraint is applied.
#' @return function that returns TRUE if point satisfies oracle, FALSE otherwise.
#' @export
#' @examples
#' #To do
#'

make_oracle <- function(oracle, region = c(-Inf,Inf)){
  
  if(!is.function(oracle)) stop("oracle must be a function")
  
  .oracle <- function(x){
    .x <- matrix(as.numeric(x), ncol = 1)
    res <- oracle(.x)
    attr(res, "region") <- region
    return(res)
    }
  
  return(.oracle)
  
}
