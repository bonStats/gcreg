#' Create discrete (data-based) orthonormal polynomial basis
#' 
#' @param x data to generat
#' @param deg Highest degree of polynomials.
#' @return polylist polynomials from polynom package defining orthonormal basis
#' @export
#' @examples
#' x <- runif(100, min = -1, max = 1)
#' p_basis <- make_disc_orthonormal_basis(x = x, deg = 9)
#' Xo <- sapply(p_basis, predict, newdata = x)
#' all.equal(diag(length(p_basis)), crossprod(Xo), tolerance = 1e-10)


make_disc_orthonormal_basis <- function(x, deg){
  
  if(max(abs(x)) > 1.00001) warning("x values outside of [-1,1] may lead to numerically unstable results, scale the data.")
  
  # data length
  len <- length(x)
  
  if(len <= deg + 1) stop("Length of x must be greater than degree of polynomial")
  
  # basis scale
  b_scale <- 1/sqrt(len)
  # inverse of sample SD
  inv_s_sd <- (var(x) * (1 - 1/len))^(-1/2)
  
  # initialise discrete orthonormal polynomial list
  p_basis <- vector(mode = "list", length = deg + 1)
  
  # initial polynomials
  p_basis[[1]] <- polynom::polynomial(b_scale)
  p_basis[[2]] <- polynom::polynomial(c(-b_scale*inv_s_sd*mean(x),b_scale*inv_s_sd))
  
  # initialise vector for leading coefficients
  p_lc <- lapply(p_basis, FUN = get_lc)
  # create rest of polys (Polynomial operations done by Depends: polynom)
  for( i in 3:(deg+1) ){
    term1 <- - calc_poly_disc_inner_prod(p1 = p_basis[[i-1]] * polynom::polynomial(c(0,1)),
                                         p2 = p_basis[[i-1]], x = x) * p_basis[[i-1]]
    term2 <- p_basis[[i-1]] * polynom::polynomial(c(0,1))
    term3 <- - p_basis[[i-2]] * p_lc[[i-2]] / p_lc[[i-1]]
    
    p_basis[[i]] <- (term1 + term2 + term3) / p_lc[[i-1]]
    p_lc[[i]] <- get_lc(p_basis[[i]]) / calc_poly_disc_norm(p_basis[[i]], x = x) 
    p_basis[[i]] <- p_basis[[i]] * p_lc[[i]]
  }
  
  class(p_basis) <- "polylist"
  
  return(p_basis) 
}


## helper functions ##

calc_poly_disc_inner_prod <- function(p1,p2,x){
  
  sum(predict(p1,x) * predict(p2,x))
  
}

calc_poly_disc_norm <- function(p,x){
# calculate discrete norm of polynomial
  sqrt(sum(predict(p,x)^2))
  
}

get_lc <- function(p) {
#  get leading coefficient
  p[length(p)]
  
}