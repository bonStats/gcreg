#' Calculate the linear distance to the monotone boundary
#' 
#' 
#' @param gam gamma iterate (should be on boundary of monotonicity)
#' @param region Region of monotonicity constraint
#' @param basis_cv converter function for polynomial basis
#' @param EPS Precision for root finding
#' @return distance to the monotone boundary (intercept of derivative - in monomial basis)
#' @export
#' @examples
#' #To do

linear_dist_to_monotone_boundary <- function(gam, region, basis_cv, EPS = 1e-06){
  
  # polynomial in monomial basis
  bet <- polynomial(basis_cv$to_mono(gam))
  
  reg <- sort(region)
  reset_reg <- F
  
  #replace infinite regions with large testable number
  if(is.infinite(reg[1])){
    reg[1] <- -100
    reset_reg <- T
  }
  
  if(is.infinite(reg[2])){
    reg[2] <- 100
    reset_reg <- T
  } 
  
  # This effects the optimisation for aim_gamma
  
  #calculate first, second derivative + roots
  Dpl <- deriv(bet)
  DDpl <-deriv(Dpl)
  DDpl_roots <- polyroot(DDpl)
  
  # getting real roots only
  re_roots <- Re(DDpl_roots)
  im_roots <- Im(DDpl_roots)
  
  DDpl_re_roots <- re_roots[abs(im_roots) < EPS]
  DDpl_relevant_roots <- DDpl_re_roots[reg[1] <= DDpl_re_roots & DDpl_re_roots <= reg[2]]
  
  crit_v <- c(reg, DDpl_relevant_roots)
  
  # find minimum
  DDpl_crit_vals <- predict(Dpl, newdata = crit_v)
  which_min_Dpl <- which.min(DDpl_crit_vals)
  location_x <- crit_v[which_min_Dpl]
  
  if(reset_reg & which_min_Dpl %in% c(1,2)) warning("linear distance to the monotone boundary based on constrained region = [", paste(reg, collapse = ","), "] instead of [", paste(region, collapse = ","), "]")
  
  # in monomial basis
  linear_dist <- DDpl_crit_vals[which_min_Dpl]
  
  return(list(dist = linear_dist, loc = location_x))
}