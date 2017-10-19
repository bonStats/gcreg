#' Create list with model specifications for polynomial constrained mixed effect model with monotonic mean and (potentially) monotonic random effects.
#' 
#' EXPIREMENTAL: See vignette for usage details.
#' 
#' @param formula Usage: y ~ x, replaced by actual names.
#' @param data data.frame containing y and x.
#' @param p_degree degree of mean polynomial.
#' @param r_degree degree of polynomial random effects.
#' @param r_constrained should the polynomial random effects be constrained? Only implemented for r_degree= 0, 1.
#' @param mcontr_region Over what region should monotonicity apply to the mean polynomial?
#' @param rcontr_region Over what region should monotonicity apply to the random effects polynomials?
#' @param group_name character column name of group in data.
#' @keywords internal
#' @return list to be passed to constrained_lmm_em() for fitting model.


make_em_model_specs <- function(formula, data, p_degree, r_degree, r_constrained = F, mcontr_region = c(-1,1), rcontr_region = c(-1,1), group_name = "grp"){
  
  vars <- all.vars(formula)
  if(length(vars) != 2) stop("Specify forumla as: 'y ~ x' and specify degree by 'p_degree'")
  
  name_changes <- setNames(list(vars[1], vars[2], group_name), c("y","x","grp"))
  
  new_data <- data %>% 
    rename_(.dots = name_changes)
  
  md_specs <-  new_data %>%
    summarise(N = n(),
              g = n_distinct(new_data$grp) # avoiding NSE in package.
    ) %>%
    mutate(p = p_degree + 1, 
           r = r_degree + 1, 
           constr_r = r_constrained
    )
  
  if(with(md_specs, r >= 3 & constr_r)){
    warning("Constrained random effects not yet implemented for r_degree >= 2, setting 'constr_r' to FALSE")
    md_specs$constr_r <- F
  }
  
  if(max(mcontr_region) < max(rcontr_region) | min(mcontr_region) > min(rcontr_region)) warning("Defining constrained range of grp to be wider than constrained range of mean not advised")
  
  return(list(specs = md_specs, 
              dat = new_data, 
              name_changes = name_changes, 
              mcontr_region = mcontr_region,
              rcontr_region = rcontr_region,
              group_ids = unique(new_data$grp)
  ))
}