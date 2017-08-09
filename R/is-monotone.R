#' Check if polynomial is monotone over desired region.
#' 
#' @param p polynomial
#' @param region numeric is polynomial monotone over this region.
#' @param type character one of "increasing", "decreasing" or NULL
#' @return TRUE if polynomial is monotone within region, FALSE otherwise.
#' @examples
#' pol <- polynomial(c(1,4,5,2))
#' is_monotone(pol, region = c(1,5))

is_monotone <- function(p, region = c(-Inf,Inf), type = NULL, EPS = 1e-05){
  
  # Checks
  pl <- as.polynomial(p)
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
  
  D_pl <- deriv(pl)
  
  # roots of derivative of pl
  roots <- polyroot(D_pl)
  roots_re <- Re(roots)
  roots_im <- Im(roots)
  
  # critical roots are those that are real and within the regions (within numerical precision)
  crit_roots <- roots_re[abs(roots_im) < EPS & a + EPS < roots_re & roots_re < b - EPS]
  
  # check multiplicity of critical roots
  crit_roots_mltplcty <- rowSums( outer(crit_roots, crit_roots, function(x, y) abs(x - y) < EPS) )
  
  # test for derivative roots
  if(any(crit_roots_mltplcty %% 2 != 0)){
    mono <- FALSE
    attr(mono,"region") <- region
    attr(mono,"type") <- "monotone increasing or decreasing"
    return(mono)
  }
  
  # inner points can be checked for monotone increasing or decreasing
  finite_reg <- region[is.finite(region)]
  
  if(length(finite_reg) == 2){
    
    inner_point <- mean(finite_reg)
    
  } else if(length(finite_reg) == 1) {
    
    inner_point <- ifelse(finite_reg == a, a + 1, b - 1)
    
  } else {
    
    inner_point <- 0
    
  }

  crit_value <- predict(D_pl,inner_point)
  
  # should i remove EPS 
  
  if(crit_value > 0){
    mono <- ifelse(usetype, "increasing" == type, TRUE)
    attr(mono,"region") <- region
    attr(mono,"type") <- "monotone increasing"
    return(mono)
  } else if(crit_value < 0){
    mono <- ifelse(usetype, "decreasing" == type, TRUE)
    attr(mono,"region") <- region
    attr(mono,"type") <- "monotone decreasing"
    return(mono)
  } else {
    mono <- TRUE
    attr(mono,"region") <- region
    attr(mono,"type") <- "unsure, check different critical value"
    return(mono)
  }

}