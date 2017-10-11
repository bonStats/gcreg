#' Create design matrix from polynomial basis.
#' 
#' @param poly_basis polylist defining the (non-monomial) basis to convert from.
#' @param x data for design matrix
#' @return Design matrix for polynomial regression in \code{ploy_basis}
#' @export

poly_basis_design_matrix <- function(poly_basis, x){
  
  if(!polynom::is.polylist(poly_basis)) warning("poly_basis should be polylist object")
  
  sapply(poly_basis, predict, newdata = x)
  
}
