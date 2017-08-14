#' Create discrete (data-based) orthonormal polynom::polynomial basis
#' 
#' @param x data to generat
#' @param deg Highest degree of polynom::polynomials
#' @return polylist (list of polynom::polynomials) defining orthonormal basis
#' @examples
#' To do
#' t(X) %*% X = diag()
#' crossprod(X) == diag()

make_disc_orthonormal_basis <- function(x, deg){
  
  # data length
  len <- length(x)
  # basis scale
  b_scale <- 1/sqrt(len)
  # inverse of sample SD
  inv_s_sd <- (var(x) * (len - 1)/len)^(-1/2)
  
  # initialise discrete orthonormal polynomial list
  p_basis <- vector(mode = "list", length = deg + 1)
  
  # initial polynomials
  p_basis[[1]] <- polynom::polynomial(b_scale)
  p_basis[[2]] <- polynom::polynomial(c(-b_scale*inv_s_sd*mean(x),b_scale*inv_s_sd))
  
  # initialise vector for leading coefficients
  p_lc <- lapply(p_basis, FUN = get_lc)
  # create rest of polys (NEED TO CHANGE "*" to operation explicitly from polynom)
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


#helper functions

calc_poly_disc_inner_prod <- function(p1,p2,x){
  
  sum(predict(p1,x) * predict(p2,x))
  
}

calc_poly_disc_norm <- function(p,x){
# calculate discrete norm of 
  sqrt(sum(predict(p,x)^2))
  
}

get_lc <- function(p) {
#  get leading coefficient
  p[length(p)]
  
}



# make_disc_orthonormal_basis <- function(x,deg){
#   # checks
#   # TO DO check deg is integer > 0, deg < x
#   # TO DO check x is numeric vector
#   
#   # data length
#   len <- length(x)
#   # basis scale
#   b_scale <- 1/sqrt(len)
#   # inverse of sample SD
#   inv_s_sd <- (var(x) * (len - 1)/len)^(-1/2)
#   
#   # initialise discrete orthonormal polynomial list
#   p_basis <- vector(mode = "list", length = deg + 1)
#   
#   # initial polynomials
#   p_basis[[1]] <- polynom::polynomial(b_scale)
#   p_basis[[2]] <- polynom::polynomial(c(-b_scale*inv_s_sd*mean(x),b_scale*inv_s_sd))
#   
#   # create rest of polys
#   for( i in 3:(deg+1) ){
#     term1 <- - calc_poly_disc_inner_prod(p1 = p_basis[[i-1]] * polynom::polynomial(c(0,1)),p2 = p_basis[[i-1]], x = x) * p_basis[[i-1]]
#     term2 <- p_basis[[i-1]] * polynom::polynomial(c(0,1))
#     term3 <- - p_basis[[i-2]] * get_lc(p_basis[[i-2]])  / get_lc(p_basis[[i-1]]) 
#     
#     p_basis[[i]] <- (term1 + term2 + term3) / get_lc(p_basis[[i-1]]) 
#     p_basis[[i]] <- p_basis[[i]] / calc_poly_disc_norm(p_basis[[i]], x = x) 
#   }
#   
#   class(p_basis) <- "polylist"
#   
#   return(p_basis) 
# }