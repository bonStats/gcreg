#' Bounce away from monotone boundary
#' 
#' This function helps if COLS gets stuck in flat spots of the constraint boundary.
#' It will bounce the starting point towards an inner point of the constrained region.
#' 
#' @param gam gamma iterate (should be on boundary of monotonicity)
#' @param Y Y data
#' @param Xo Orthonormal X design matrix 
#' @param poly_basis Polynomial basis for gamma
#' @param oracle_fun Monotonicity oracle function
#' @param region Region of monotonicity constraint
#' @param basis_cv converter function for polynomial basis
#' @param control control list created by \code{cols_control}
#' @return new gamma estimate or original if no better soluton found
#' @export
#' @examples
#' #To do
#'

bounce_monotone <- function(gam, Y, Xo, poly_basis, oracle_fun, region, basis_cv, control){

  Dpoly_basis <- deriv(poly_basis)
  
  frss <- function(gam) sum((Y - Xo %*% gam)^2)
  
  old_rss <- -Inf
  new_rss <- 0
  z <- gam
  
  bb <- 1
  
  while(abs(new_rss - old_rss) > control$tol_bounce & bb <= control$maxit_bounces){
    
    old_rss <- new_rss <- frss(z)
    
    loc <- linear_dist_to_monotone_boundary(z, region = region, basis_cv = basis_cv)$loc
    bump_dir <- sapply(Dpoly_basis, predict, newdata = loc)
    bump_norm <- bump_dir / sqrt(sum(bump_dir^2))
    max_bump <- bump_norm
    new_start <- z + bump
    
    while(!oracle_fun(new_start)){
      max_bump <- max_bump / 2
      new_start <- z + max_bump
      new_start_neg <- z - max_bump
      if(oracle_fun(new_start_neg)){
        new_start <- new_start_neg
      }
    }
    
    bump <- max_bump / 2
    
    for(kk in 1:control$maxit_bumps){
      
      temp_z <- optim_cols(par = new_start, Y = Y, X = Xo, oracle_fun = oracle_fun, control = control)
      
      rss_temp <- frss(temp_z)
      
      if(rss_temp < old_rss) {
        z <- temp_z
        new_rss <- rss_temp
        break
      }
      
      bump <- bump / (2^(kk))
      new_start <- z + bump
      
    }
    
    bb <- bb + 1
  }
  
  return(z)
  
}
