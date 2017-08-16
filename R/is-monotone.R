#' Check if polynomial is monotone over desired region.
#' 
#' @param p polynomial
#' @param region numeric is polynomial monotone over this region.
#' @param type character one of "increasing", "decreasing" or NULL
#' @return TRUE if polynomial is monotone within region, FALSE otherwise.
#' @export
#' @examples
#' pol <- polynomial(c(1,4,5,2))
#' is_monotone(pol, region = c(1,5))

is_monotone <- function(p, region = c(-Inf,Inf), type = NULL, EPS = 1e-05){
  
  # Checks
  pl <- polynom::as.polynomial(p)
  if(length(region) != 2 | ! is.numeric(region)){
    r <- dput(region)
    stop("Argument in wrong format. region = c(-1,1) for example.")
  }
  if(EPS <= 0) {
    stop("Argument in wrong format. EPS should be a small positive number.")
  }
  if(is.null(type)){
    usetype <- F
  } else {
    usetype <- T
    type <- match.arg(type, choices = c("increasing","decreasing"))
  }

  a <- min(region)
  b <- max(region)
  
  D_pl <- polynom:::deriv.polynomial(pl)
  
  # roots of derivative of pl
  roots <- polyroot(D_pl)
  roots_re <- Re(roots)
  roots_im <- Im(roots)
  
  # critical roots are those that are real and within the regions (within numerical precision)
  crit_roots <- roots_re[abs(roots_im) < EPS & a + EPS < roots_re & roots_re < b - EPS]
  
  # check multiplicity of critical roots
  crit_roots_mltp <- rowSums( outer(crit_roots, crit_roots, function(x, y) abs(x - y) < EPS) )
  
  # test for derivative roots having multiplicity 2
  
  mono <- all(crit_roots_mltp %% 2 == 0)
  attr(mono,"region") <- region
  
  return(mono)

}