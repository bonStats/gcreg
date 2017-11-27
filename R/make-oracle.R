#' Create oracle function for constrained polynomial fitting
#' 
#' For use with \code{\link{cpm}}. See getting started vignette for details.
#' 
#' @param oracle oracle function, returns TRUE if point is in constrained region, FALSE otherwise.
#' @param region region over which constraint is applied.
#' @return function that returns TRUE if point satisfies oracle, FALSE otherwise.
#' @export
#' @examples
#' \dontrun{
#' constraint_oracle <- make_oracle( function(b) { b[1,] > 1 })
#' # In general, variable b should be treated as a matrix.
#' }
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
