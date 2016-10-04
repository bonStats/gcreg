#' Check if polynomial is monotone over desired region.
#' 
#' @param p polynomial
#' @param region numeric is polynomial monotone over this region.
#' @param which character one of "increasing", "decreasing" or NULL
#' @return TRUE if polynomial is monotone within region, FALSE otherwise.
#' @examples
#' pol <- polynomial(c(1,4,5,2))
#' is_monotone(pol, region = c(1,5))

is_monotone <- function(p, region = c(-Inf,Inf), which = NULL, EPS = 1e-05){
  
  # Checks
  pp <- as.polynomial(p)
  if(length(region) != 2 | ! is.numeric(region)){
    r <- dput(region)
    stop("Argument in wrong format. region = c(-1,1) for example.")
  }
  if(EPS <= 0) {
    stop("Argument in wrong format. EPS should be a small positive number.")
  }
  if(is.null(which)){
    use.which <- F
  } else {
    use.which <- T
    .which <- match.arg(which, choices = c("increasing","decreasing"))
  }

  a <- min(region)
  b <- max(region)
  
  D_pp <- deriv(pp)
  
  # roots of derivative of pp
  roots <- polyroot(D_pp)
  roots_re <- Re(roots)
  roots_im <- Im(roots)
  
  # critical roots are those that are real and within the regions (within numerical precision)
  crit_roots <- roots_re[abs(roots_im) < EPS & a - EPS < roots_re & roots_re < b + EPS]
  # check multiplicity of critical roots
  crit_roots_mltplcty <- rowSums( outer(crit_roots, crit_roots, function(x, y) abs(x - y) < EPS) )
  
  # test for derivative roots
  if(any(crit_roots_mltplcty %% 2 != 0)){
    return(FALSE)
  }
  
  # boundaries and inner points can be checked for monotone increasing or decreasing
  crit_boundary <- region[is.finite(region)]
  pos_inner_pts <- c(mean(region), a + 100, b - 100)
  fin_inner_pts <- pos_inner_pts[is.finite(pos_inner_pts)]
  inner_points <- fin_inner_pts[a < fin_inner_pts & fin_inner_pts < b]
  crit_values <- predict(D_pp,c(crit_boundary,inner_points))

  
  if(all(crit_values > - EPS)){
    mono <- ifelse(use.which, "increasing" == .which, TRUE)
    attr(mono,"region") <- region
    attr(mono,"type") <- "monotone increasing"
    return(mono)
  } else if(all(crit_values < EPS)){
    mono <- ifelse(use.which, "decreasing" == .which, TRUE)
    attr(mono,"region") <- region
    attr(mono,"type") <- "monotone decreasing"
    return(mono)
  } else {
    return(FALSE)
  }

}