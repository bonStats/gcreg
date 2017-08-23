#' Create controls attributes which tune the COLS algorithm.
#' 
#' @param method "simple" or "best-step".
#' @param maxit maximum number of iterations.
#' @param tol tolerance dictating when convergence is reached.
#' @param step_start the initial step size taken.
#' @param step_increment how step size increases over iterations.
#' @return a list of control parameters.
#' @export
#' @examples
#' #TO DO


cols_control <- function(method = "best-step", maxit = 200, tol = 1e-05, step_start = 0.7, step_increment = 0.05){
  
  cntrl <- list(
    method = match.arg(method, c("up-walk","down-walk","best-step")), 
    maxit = ifelse(maxit > 0, maxit, 200), 
    tol = ifelse(tol > 0, tol, 1e-05), 
    step_start = ifelse(step_start > 0 & step_start <= 1, step_start, 0.7), 
    step_increment = ifelse(step_increment >= 0 & step_increment <= 1, step_increment, 0.05)
  )
  
  if(maxit <= 0){
    warning("maxit must be greater than 0, reset to ", cntrl$maxit)
  }
  
  if(tol <= 0){
    warning("tol must be greater than 0, reset to ", cntrl$tol)
  }
  
  if(step_start <= 0 | step_start > 1){
    warning("step_start should be in set (0,1], reset to ", cntrl$step_start)
  }
  
  if(step_increment < 0 | step_increment >= 1){
    warning("step_increment should be in set [0,1), reset to ", cntrl$step_increment)
  }
  
  return(cntrl)
  
}
