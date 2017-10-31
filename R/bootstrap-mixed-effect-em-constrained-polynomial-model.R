#' Bootstrap for fitting polynomial constrained mixed effect model with monotonic mean and (potentially) monotonic random effects.
#' 
#' EXPIREMENTAL: For use with block-random effects (non-crossed), e.g. subject-specific random effects.
#' This function is under development and may change without warning. See the mixed effects vignette for usage details.
#' 
#' Uses case bootstrapping, where a sample of groups is chosen with replacement and the new estimates are used.
#' 
#' @param model model specification made by gcreg:::make_em_model_specs().
#' @param N 
#' @param maxit maximum number of iterations for each resample.
#' @param tol relative tolerance for convergence of parameter values for each resample.
#' @param start start values for parameters for each resample. Defaults to NULL.
#' @param verbose output current iteration and quasi-log-likelihood value in each resample.
#' @param save_steps should each em step be saved and return in a list?
#' @param check_up should the algorithm check if the quasi-log-likelihood is improving each iteration?
#' @param verbose_boot Output when a bootstrap resample is finished.
#' @param parallel Execute in parallel. This will turn off verbose settings. ONLY IMPLEMENTED FOR LINUX/MAC
#' @param cores Number of cores to use.
#' @keywords internal
#' @return list including estimated mixed effects model.
#' 

case_bootstrap_constrained_lmm_em <- function(model, N, maxit = 25, tol = 1e-02, start, verbose = F, save_steps = F, check_up = F, verbose_boot = verbose, parallel = F, cores = 1){
  
  new_fits <- list()
  
  if(!missing(start)){
    start_list <- list(
      gam = start$gam,
      sig2 = start$sig2,
      omeg = start$omeg
    ) 
  } else {
    start_list <- NULL
  }
  
  make_model_from_sampled_data <- function(model){
    
    .md <- model
    
    sampled_groups <- sample(model$group_ids, replace = T)
    
    not_old_group_col <- which(colnames(model$dat) != "grp")
    
    newdata <- model$dat[0,not_old_group_col]
    
    for(gr in 1:length(sampled_groups)){
      
      newdata <- rbind(
        cbind(model$dat[model$dat$grp == sampled_groups[gr],not_old_group_col], grp = gr),
        newdata
        )
      
    }
    
    .md$dat <- newdata
    .md$group_ids <- 1:length(sampled_groups)
    
    return(.md)
    
  }
  
  if(!parallel){
  
    for(i in 1:N){
      
      new_md <- make_model_from_sampled_data(model)
      
      new_fits[[i]] <- constrained_lmm_em(model = new_md, 
                                          maxit = maxit, 
                                          tol = tol, 
                                          start = start_list, 
                                          verbose = verbose, 
                                          save_steps = save_steps, 
                                          check_up = check_up)
      
      if(verbose_boot) cat("\nBootstrap it =", i, "\n")
      if(!new_fits[[i]]$converged) warning(new_fits[[i]]$status, ". Consider increasing maxit")
      
    }
  
  } else {
    #parallel
    
    gen_boot_clmm_em <- function(model, maxit, tol, start, check_up){
      md <- model
      .maxit <- maxit
      .tol <- tol
      .start <- start
      .check_up <- check_up
      return(
        function(...) {constrained_lmm_em(model = make_model_from_sampled_data(md),
                         maxit = .maxit,
                         tol = .tol,
                         start = .start,
                         verbose = F,
                         save_steps = F,
                         check_up = .check_up)}
        )
    }
    
    boot_clmm_em <- gen_boot_clmm_em(model = model, maxit = maxit, tol = tol, start = start_list, check_up = check_up)
    
    new_fits <- mclapply(1:N, FUN = boot_clmm_em, mc.cores = cores)
    
  }
  
  return(new_fits)
  
}