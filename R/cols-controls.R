#' Create controls attributes which tune the COLS algorithm.
#' 
#' @param method one of "best-step", "once-best-step", "up-walk", "down-walk", "avoid-boundary"
#' @param maxit maximum number of iterations.
#' @param tol tolerance dictating when convergence is reached.
#' @param step_start the initial step size taken.
#' @param step_increment how step size increases over iterations.
#' @param maxit_bounces Maximum number of bounces away from the solution.
#' @param maxit_bumps  Maximum number of bumps to test for each solution.
#' @param tol_bounce tolerance level bounce iterations.
#' @return a list of control parameters.
#' @export
#' @examples
#' cols_control(maxit = 100)
#' 


cols_control <- function(method = "best-step", maxit = 200, tol = 1e-05, step_start = 0.7, step_increment = 0.05, maxit_bounces = 2, maxit_bumps = 5, tol_bounce = tol){
  
  cntrl <- list(
    method = match.arg(method, c("once-best-step","up-walk","down-walk","best-step","avoid-boundary")), 
    maxit = ifelse(maxit > 0, maxit, 200), 
    tol = ifelse(tol > 0, tol, 1e-05), 
    step_start = ifelse(step_start > 0 & step_start <= 1, step_start, 0.7), 
    step_increment = ifelse(step_increment >= 0 & step_increment <= 1, step_increment, 0.05),
    maxit_bounces = ifelse(maxit_bounces > 0, maxit_bounces, 2), 
    maxit_bumps = ifelse(maxit_bumps > 0, maxit_bumps, 5),  
    tol_bounce = ifelse(tol_bounce > 0, tol_bounce, 1e-04)
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
  
  if(maxit_bounces <= 0){
    warning("maxit_bounces must be greater than 0, reset to ", cntrl$maxit_bounces)
  }
  
  if(maxit_bumps <= 0){
    warning("maxit_bumps must be greater than 0, reset to ", cntrl$maxit_bumps)
  }
  
  if(tol_bounce <= 0){
    warning("tol_bounce must be greater than 0, reset to ", cntrl$tol_bounce)
  }
  
  return(cntrl)
  
}
