#' Generate functions to convert regression coefficients between a given polynomial basis (e.g. discrete orthonormal) and standard (monomial) polynomial bases.
#' 
#' @param poly_basis polylist defining the (non-monomial) basis to convert from.
#' @return list of length 2 with functions having names 'to_mono' (conversion to monomial basis) and 'to_ortho' (conversion to orthonormal basis)
#' @export

gen_poly_basis_converters <- function(poly_basis){
  
  if(!polynom::is.polylist(poly_basis)) warning("poly_basis should be polylist object")
  
  dm <- length(poly_basis)
  
  conv_mat <- matrix(0, nrow = dm, ncol = dm)
  
  # add polys by column to matrix
  conv_mat[upper.tri(conv_mat, diag = T)] <- c(unlist(poly_basis), recursive = T)
  
  # convert from gamma (orthonormal or other basis) to beta (monomial)
  to_beta <- function(g) {
    
    if(is.matrix(g)) {
      if(ncol(g) != 1 | nrow(g) != dm) warning("g should have 1 column and ", dm, " rows")
    } else {
      if(length(g) != dm) warning("g should be length", dm)
    }
    
    conv_mat %*% g
    
  }
  
  # convert from beta (monomial) to gamma (orthonormal or other basis)
  to_gamma <- function(b) {
    
    if(is.matrix(b)) {
      if(ncol(b) != 1 | nrow(b) != dm) warning("b should have 1 column and ", dm, " rows")
    } else {
      if(length(b) != dm) warning("b should be length ", dm)
    }
    
    backsolve(r = conv_mat, x = b, k = dm, upper.tri = T)
    
  }
  
  return(list(to_mono = to_beta, to_ortho = to_gamma))
  
}
