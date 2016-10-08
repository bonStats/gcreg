#' Search for boundary of permissible region between current value and desired value.
#' 
#' @param oracle_fun oracle function to determine membership to constrained space.
#' @return Function to conduct line search with specified parameters.
#' @examples
#' mono_ls <- gen_line_search(oracle_fun = is_monotone)
#' mono_ls(cur = 0:5, aim = c(0,1,2,3,4,-5),index = 6)
#' 
#' mono_ls <- gen_line_search(oracle_fun = is_monotone, region = c(0,1))
#' mono_ls(cur = 0:5, aim = c(0,1,2,3,4,-5),index = 6)

gen_line_search <- function(oracle_fun, EPS = 1e-05, max_it = 200, ...) {
  
  function(cur, index, aim){
    
    .cur <- matrix(cur,ncol = 1)
    .aim <- matrix(aim,ncol = 1)
    
    i <- 1
    while(i < max_it & abs(.cur[index,] - .aim[index,]) > EPS){
      i <- i + 1
      test_par <- replace(.cur, index, (.cur[index,] + .aim[index,])/2)
      
      if(oracle_fun(test_par,...)){
        .cur <- test_par
      } else {
        .aim <- test_par
      }
      
    }
    
    attr(.cur,"convg") <- abs(.cur[index,] - .aim[index,]) <= EPS
    attr(.cur,"iter") <- i
    
    return(.cur)
    
  }
}

mono_ls <- gen_line_search(oracle_fun = is_monotone, region = c(0,1))
mono_ls(cur = 0:5, aim = c(0,1,2,3,4,-5),index = 6)

is_monotone(c(0:4,-4.99999),region = c(0,1))
