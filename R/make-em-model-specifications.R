#' Specify model for a monotone-constrained polynomial LMM
#' 
#' EXPIREMENTAL: Creates a list with model specifications for polynomial constrained mixed effect model with monotonic mean and (potentially) monotonic random effects.
#' For use with \code{\link{constrained_lmm_em}} and \code{\link{constrained_lmm_mcem}}. See mixed effects vignette for usage details.
#'
#' The interface and default actions of this function are under development and may change without warning. See the mixed effects vignette for usage details.
#' 
#' @param formula Usage: y ~ x, replaced by actual names.
#' @param group_name character column name of group in data.
#' @param data data.frame containing y and x.
#' @param p_degree degree of mean polynomial.
#' @param r_degree degree of polynomial random effects.
#' @param r_constrained should the polynomial random effects be constrained?
#' @param mcontr_region Over what region should monotonicity apply to the mean polynomial?
#' @param rcontr_region Over what region should monotonicity apply to the random effects polynomials? Defaults to (and generally should be) \code{mcontr_region}.
#' @keywords internal
#' @return list speficying model to be passed to fitting methods.


make_em_model_specs <- function(formula, group_name = "grp", data, p_degree, r_degree, r_constrained = F, mcontr_region = c(-1,1), rcontr_region = mcontr_region){
  
  vars <- all.vars(formula)
  if(length(vars) != 2) stop("Specify forumla as: 'y ~ x' and specify degree by 'p_degree'")
  
  old_vars <- c(vars[1], vars[2], group_name)
  
  name_changes <- setNames(c("y","x","grp"), old_vars)
  
  new_data <- data[,old_vars]
  colnames(new_data) <- name_changes[old_vars]
  
  md_specs <-  list()
    
  md_specs$N <- nrow(new_data)
  md_specs$g <- length(unique(new_data$grp))
  md_specs$p = p_degree + 1 
  md_specs$r = r_degree + 1 
  md_specs$constr_r = r_constrained
  
  if(max(mcontr_region) < max(rcontr_region) | min(mcontr_region) > min(rcontr_region)) warning("Defining constrained range of grp to be wider than constrained range of mean not advised")
  
  return(list(specs = md_specs, 
              dat = new_data, 
              name_changes = name_changes, 
              mcontr_region = mcontr_region,
              rcontr_region = rcontr_region,
              group_ids = unique(new_data$grp)
  ))
}
