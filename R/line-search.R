#' Search for boundary of permissible region between current value and desired value.
#' 
#' @param oracle_fun oracle function to determine membership to constrained space.
#' @param EPS precision of line search algorithm
#' @param maxit maximum number of iterations
#' @param ... extra arguments. Currently ignored.
#' @keywords internal
#' @return Function to conduct line search with specified parameters.


gen_line_search <- function(oracle_fun, EPS = 1e-05, maxit = 200, ...) {
  
  oracle_ <- oracle_fun
  
  function(index, cur, aim){
    
    cur_ <- matrix(cur,ncol = 1)
    aim_ <- matrix(aim,ncol = 1)
    
    i <- 1
    while(i < maxit & abs(cur_[index,] - aim_[index,]) > EPS){
      i <- i + 1
      test_par <- replace(cur_, index, (cur_[index,] + aim_[index,])/2)
      
      if(oracle_(test_par)){
        cur_ <- test_par
      } else {
        aim_ <- test_par
      }
      
    }
    
    attr(cur_,"convg") <- abs(cur_[index,] - aim_[index,]) <= EPS
    attr(cur_,"iter") <- i
    
    return(cur_)
    
  }
}
