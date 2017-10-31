#' Fit polynomial constrained mixed effect model with monotonic mean and monotonic random effects.
#' 
#' EXPIREMENTAL: For use with block-random effects (non-crossed), e.g. subject-specific random effects.
#' This function is under development and may change without warning. See the mixed effects vignette for usage details.
#' 
#' @param model model specification made by make_em_model_specs().
#' @param maxit maximum number of iterations.
#' @param tol relative tolerance for convergence of parameter values.
#' @param start start values for parameters. Defaults to NULL.
#' @param verbose output current iteration and quasi-log-likelihood value.
#' @param save_steps should each em step be saved and return in a list?
#' @param check_up should the algorithm check if the quasi-log-likelihood is improving each iteration?
#' @param parallel Execute Monte Carlo expectations in parallel. This will turn off verbose settings for expectations. ONLY IMPLEMENTED FOR LINUX/MAC
#' @param cores Number of cores to use.
#' @param ni_ssize Number of realisations to used for numerical integration.
#' @keywords internal
#' @return list including estimated mixed effects model.

constrained_lmm_mcem <- function(model, maxit = 200, tol, start = NULL, verbose = F, save_steps = F, check_up = F, parallel = F, cores = 1, ni_ssize = 20000){
  
  if(missing(model)) stop("model argument must be specified. Use gcreg:::make_em_model_specs().")
  
  if(with(model$specs, r < 3 | !constr_r)) stop("For unconstrained random effects or r < 3, please use gcreg:::constrained_lmm_em() instead.")
  
  # Algorithm constants, functions
  
  # create orthonormal polynomial basis
  pl_base <- make_disc_orthonormal_basis(x = model$dat$x, deg = model$specs$p - 1)
  # gradient of basis
  Dpl_base <- deriv(pl_base)
  
  # converter matrices between bases for mean and random effects
  mean_pl_cv <- gen_poly_basis_converters(pl_base)
  grp_pl_cv <- gen_poly_basis_converters(pl_base[1:model$specs$r])
  
  # oracle function
  monotone_oracle <- function(b) is_monotone(p = mean_pl_cv$to_mono(b), region = model$mcontr_region)
  
  # Orthonormal X matrix from basis and x - mean curve
  Xo <- poly_basis_design_matrix(poly_basis = pl_base, x = model$dat$x)
  
  # observation matrix
  Y <- model$dat$y
  
  # Orthonormal Z matrices - group curve
  Zo_i <- list() 
  for(gr in 1:model$specs$g){
    Zo_i[[gr]] <- Xo[model$dat$grp == model$group_ids[gr], 1:model$specs$r, drop=F]
  }
  
  Zo <- bdiag(Zo_i)
  
  # Model calculated constants
  Xo_t_Y <- crossprod(Xo, Y) # unconstrained solution for fixed effect regression
  Xo_t_Zo <- crossprod(Xo, Zo)
  Zo_t_Xo <- t(Xo_t_Zo)
  Zo_t_Y <- crossprod(Zo, Y)
  Zo_t_Zo <- crossprod(Zo)
  
  # mv norm (independent) random sample for stochastic evaluation of penalty integral in M-step
  mvnorm_rs <- matrix(rnorm(ni_ssize * model$specs$r), nrow = ni_ssize)
  
  # trace level for optim
  tr_level <- 2
  
  #### Functions ####
  
  #### E-step functions ####
  
  E_step <- function(em_list){
    
    .em <- em_list
    
    ## last iterates' variance matrices
    R <- bdiag(rep(list(.em$sig2),times = model$specs$N))
    H <- make_L_mat(.em$omeg) %>% tcrossprod()
    G <- bdiag(rep(list(H),times = model$specs$g))
    
    ## unconstrained mean/variance matrices of random effects
    M_U_uc <- mean_re_unconstrained(gam = .em$gam, G = G, R = R)
    V_U_uc <- var_re_unconstrained(G = G, R = R)
    
    if(with(model$specs, r >= 3 & constr_r)){ # constrained and r>2
      
      ## constrained mean/variance matrices of random effects
      mv_list <- make_uncon_mean_var_list(re_mean = M_U_uc, re_var = V_U_uc)
      
      if(parallel){
        
        re_iter_list <- mclapply(mv_list, FUN = sample_re_fun, gam_mean = .em$gam, mc.cores = cores)
        
      } else {
        
        re_iter_list <- list()
        
        if(verbose) cat("---- r mc ----\n")
        
        for(i in 1:model$specs$g){
          re_iter_list[[i]] <- sample_re_fun(mean_var_ls = mv_list[[i]], gam_mean = .em$gam)
        }
        
        if(verbose) cat("\n--------------\n")
      }
      
      # constrained mean and variance matrices
      .em$re_mean <- matrix(c(lapply(re_iter_list, colMeans), recursive = T), ncol = 1)
      .em$re_var  <- bdiag(lapply(re_iter_list, var))
      
    } else {
      
      stop("For unconstrained random effects or r < 3, please use gcreg:::constrained_lmm_em() instead.")
      
    } 
    
    return(.em)
    
  }
  
  # helper functions - E step
  
  # make lower diagonal matrix from omega coefficients of lower-diagonal log-Cholesky decomposition
  make_L_mat <- function(omeg){
    L <- matrix(nrow = model$specs$r, ncol = model$specs$r)
    diag(L) <- exp(omeg[1:model$specs$r])
    L[upper.tri(L)] <- omeg[(model$specs$r+1):length(omeg)] # populates along row...
    L[lower.tri(L)] <- 0
    return(L)
  }
  
  mean_re_unconstrained <- function(gam, G, R){
    G_Zo_t <- tcrossprod(G, Zo)
    G_Zo_t %*% solve(Zo %*% G_Zo_t + R, Y - Xo %*% gam)
  }
  
  var_re_unconstrained <- function(G, R){
    G_Zo_t <- tcrossprod(G, Zo)
    G - G_Zo_t %*% solve(Zo %*% G_Zo_t + R, t(G_Zo_t))
  }
  
  # returns in terms of orthonormal basis
  linear_dist_to_monotone_boundary_em <- function(gam, EPS = 1e-06){
    
    # polynomial in monomial basis
    bet <- polynomial(mean_pl_cv$to_mono(gam))
    
    #check if monotone
    is_mono <- is_monotone(bet, region = model$mcontr_region)
    
    if(!is_mono & any(is.infinite(model$mcontr_region))) stop("linear_dist_to_monotone_boundary_em() doesn't deal with infinite support and non-monotone polynomials")
    # This effects the optimisation for aim_gamma
    
    #calculate first, second derivative + roots
    Dpl <- deriv(bet)
    DDpl <-deriv(Dpl)
    DDpl_roots <- polyroot(DDpl)
    
    # getting real roots only
    re_roots <- Re(DDpl_roots)
    im_roots <- Im(DDpl_roots)
    
    DDpl_re_roots <- re_roots[abs(im_roots) < EPS]
    DDpl_relevant_roots <- DDpl_re_roots[min(model$mcontr_region) <= DDpl_re_roots & DDpl_re_roots <= max(model$mcontr_region)]
    
    crit_v <- c(model$mcontr_region[is.finite(model$mcontr_region)], DDpl_relevant_roots)
    
    # find minimum
    DDpl_crit_vals <- predict(Dpl, newdata = crit_v)
    which_min_Dpl <- which.min(DDpl_crit_vals)
    location_x <- crit_v[which_min_Dpl]
    # in monomial basis
    linear_dist <- DDpl_crit_vals[which_min_Dpl]
    
    # # in desired basis
    # linear_dist_o <- grp_pl_cv$to_ortho(c(0,linear_dist))[2]
    
    return(list(dist = linear_dist, loc = location_x))
  } 
  
  rw_sample_mono <- function(re_start, re_mono_oracle, re_mean, re_var, scal, maxit, burnin = 1000, thin = 1){
    
    re_st <- matrix(0, nrow = maxit + burnin, ncol = length(re_start))
    re_st[1,] <- re_start
    
    accepted <- 0 
    
    prop_var <- (0.5)^2 * scal * re_var
    
    for(i in 2:(burnin + maxit)){
      
      re_st[i,] <- mvtnorm::rmvnorm(n = 1, mean = re_st[i-1,], sigma = prop_var)
      
      if(re_mono_oracle(re_st[i,])){
        accept_prop <- exp(mvtnorm::dmvnorm(x = re_st[i,], mean = re_mean, sigma = re_var, log = T) - 
                             mvtnorm::dmvnorm(x = re_st[i-1,], mean = re_mean, sigma = re_var, log = T)
        )
        
        if(runif(1) < accept_prop){ #accept
          accepted <- (accepted + 1) * (i > burnin)
        } else {
          re_st[i,] <- re_st[i-1,] # reject
        }
        
      } else {
        re_st[i,] <- re_st[i-1,] # reject
      }
      
    }
    
    re_st <- re_st[burnin + seq(0, maxit, by = thin),]
    
    attr(re_st,"accept_rate") <- accepted/maxit
    
    return(re_st)
    
  }
  
  make_uncon_mean_var_list <- function(re_mean, re_var){
    
    re_list <- vector(mode = "list", length = model$specs$g)
    
    for(i in 1:model$specs$g){
      
      which_el <- (model$specs$r * (i - 1) + 1):(model$specs$r * i)
      
      re_list[[i]] <- list(
        re_mean = re_mean[which_el,],
        re_var = as.matrix(re_var[which_el,which_el])
      )
      
    }
    
    return(re_list)
  }
  
  
  gen_sample_random_effects <- function(model, min_sample, rjs_maxit, thin_rw = 5, verbose = F){
    
    .model_specs <- model$specs
    .rcontr_region <- model$rcontr_region
    .mean_pl_cv <- mean_pl_cv
    
    function(mean_var_ls, gam_mean){
      
      .re_oralce <- function(gam_re){
        
        subject_gam <- gam_mean + c(gam_re, rep(0, times = with(.model_specs, p - r)))
        
        gcreg::is_monotone(p = .mean_pl_cv$to_mono(subject_gam), region = .rcontr_region)
        
      }
      
      total <- 0
      i <- 0
      
      rs <- matrix(nrow = 0, ncol = .model_specs$r)
      
      while(total < min_sample & i <= rjs_maxit){
        
        new_rs_rjs <- mvtnorm::rmvnorm(n = min_sample * 2,
                                       mean = mean_var_ls$re_mean,
                                       sigma = mean_var_ls$re_var
        )
        
        keep <- apply(new_rs_rjs, MARGIN = 1, .re_oralce)
        
        rs <- rbind(rs, new_rs_rjs[keep,])
        
        total <- total + sum(keep)
        i <- i + 1
        
        if(total < (i * min_sample / rjs_maxit)){ break } 
        
      }
      
      if(verbose) cat("rjs:",total)
      
      if(total < min_sample){
        
        if(nrow(rs) > 0){
          re_start <- rs[nrow(rs),]
        } else {
          re_start <- mean_var_ls$re_mean + c(0,1,rep(0,.model_specs$r - 2))
          while(!.re_oralce(re_start)){
            re_start <- re_start + c(0,1,rep(0,.model_specs$r - 2))
          }
        }
        
        new_rs_rw <- rw_sample_mono(re_start = re_start, 
                                    re_mono_oracle = .re_oralce, 
                                    re_mean = mean_var_ls$re_mean, 
                                    re_var = mean_var_ls$re_var, 
                                    scal = 0.5, 
                                    maxit = min_sample * thin_rw, 
                                    thin = thin_rw)
        
        rs <- rbind(rs, new_rs_rw)
        
      }
      
      if(verbose) cat(" rws:",nrow(rs) - total,";\t")
      
      return(rs)
      
    }
    
  }
  
  sample_re_fun <- gen_sample_random_effects(model = model, 
                                             min_sample = 1000, 
                                             rjs_maxit = 5, 
                                             thin_rw = 5, 
                                             verbose = verbose & !parallel)
  
  ####

  
  #### M-step functions ####
  
  M_step <- function(em_list, step_size = 0.1){ #could remove step size
    
    .em <- em_list
    
    par_st <- gamma_init
    
    if(verbose) cat("*optim aim gamma*\n")
    
    Dql_gamma_wrapper <- gen_calc_Dql_gamma_wrapper(em_list = .em)
    
    aim_gamma <- optim(par = em_list$gam,
                       fn = ql_hood_gamma_wrapper,
                       gr = Dql_gamma_wrapper,
                       em_list = em_list,
                       method = "BFGS",
                       control = list(maxit = 4, trace = tr_level)
    )$par
    
    if(verbose) cat("*optim cols gamma*\n")
    
    .em$gam <- optim_cols(par = em_list$gam, 
                          Y = Y - Zo %*% em_list$re_mean,
                          X = Xo,
                          oracle_fun = monotone_oracle,
                          aim_gamma = aim_gamma
    )
    
    # avoid getting stuck on boundary
    .em <- bounce_monotone_em(em_list = .em, aim_gam = aim_gamma, maxit_bounces = 4, maxit_bumps = 10, tol = 1e-02)
    
    
    # errors
    new_resid <- Y - Xo %*% .em$gam - Zo %*% .em$re_mean
    # residual sum of squares (accounting for mean of random effects)
    new_RSS <- (crossprod(new_resid))[1,1]
    
    # Optimise sigma R
    .em$sig2 <- (trc(Zo_t_Zo %*% .em$re_var) + new_RSS) / model$specs$N
    
    # Optimise omegas
    if(model$specs$r == 1){
      
      .em$omeg <- (trc(.em$re_var) + crossprod(.em$re_mean))[1,1] / model$specs$g
      
    } else { 
      
      if(verbose) cat("*optim omega*\n")
      
      Dql_omega_wrapper <- gen_calc_Dql_omega_wrapper(em_list = .em)
      
      opt_omegas <- optim(par = .em$omeg,
                          fn = ql_hood_omega_wrapper,
                          gr = Dql_omega_wrapper,
                          em_list = .em,
                          method = "L-BFGS-B",
                          lower = with(model$specs, c(rep(-10, times = r), rep(-20, times = (r - 1) * r / 2))),
                          upper = with(model$specs, c(rep( 10, times = r), rep( 20, times = (r - 1) * r / 2))),
                          control = list(maxit = 6, trace = tr_level)
                          # lower/upper bounds stop numerical difficulty with solve(L). Should generically decide on these based on problem size.
      )
      
    }
    
    .em$omeg <- opt_omegas$par
    
    return(.em)
  }
  
  
  # helper functions - M step
  
  Dpenalty_gamma <- function(gam_param, omeg, boundary_rs, total_rs){
    
    L <- make_L_mat(omeg)
    U <- boundary_rs %*% t(L) # mv normal random sample with mean = 0 and var = H = L %*% t(L)
    
    re_oralce <- function(gam_re, gam){
      
      subject_gam <- gam + c(gam_re, rep(0, times = with(model$specs, p - r)))
      
      gcreg::is_monotone(p = mean_pl_cv$to_mono(subject_gam), region = model$rcontr_region)
      
    }
    
    penalty_fun <- function(gam){
      in_region <- sum(apply(U, MARGIN = 1, re_oralce, gam = gam))
      in_region / total_rs
    } 
    
    return(
      numDeriv::grad(func = penalty_fun, x = gam_param, method = "simple", method.args=list(eps=1e-02))
    )
    
  }
  
  # quasi-likelihood derivative wrt gamma
  gen_calc_Dql_gamma <- function(em_list_init){

    # find random sample that are on boundary to speed up computation.
    L_init <- make_L_mat(em_list_init$omeg)
    U <- mvnorm_rs %*% t(L_init)
    subject_gams <- apply(U, MARGIN = 1, FUN = function(re){em_list_init$gam + c(re, rep(0, times = with(model$specs, p - r)))})
    
    U_lin_dist_to_boundary <- apply(subject_gams, MARGIN = 2, 
                                    FUN = function(re_gam){ linear_dist_to_monotone_boundary_em(re_gam)$dist} )
    
    sample_on_boundary <- order(U_lin_dist_to_boundary)[1:floor(nrow(U) * 0.01)] # use 1% closest to boundary.
    
    penalty_val <- exp(log_normalising_penalty(L = L_init, gam_mean = em_list_init$gam, all_integrals = F))
    
    function(em_list){
    
      Drss <- with(em_list, - 2 * (Xo_t_Y - gam - Xo_t_Zo %*% re_mean) / sig2)
  
      if(with(model$specs, r > 2 & constr_r)){
        
        Dpen <- Dpenalty_gamma(gam_param = em_list$gam, 
                               omeg = em_list$omeg,
                               boundary_rs = mvnorm_rs[sample_on_boundary,], 
                               total_rs = nrow(mvnorm_rs)
        )
        
        Dlog_pen <- model$specs$g * Dpen / penalty_val
        
      }
  
      return(matrix(Drss + 2 * Dlog_pen))
      
    }

  }
  
  # calc_Dc_beta <- function(gam, min_x){ # using evelope theorem
  #   
  #   sapply(Dpl_base, predict, newdata = min_x)
  #   
  # }
  
  gen_calc_Dql_gamma_wrapper <- function(em_list){
    
    calc_Dql_gamma <- gen_calc_Dql_gamma(em_list)
    .em <- em_list
    
    function(gam, ...){
      .em$gam <- gam
      calc_Dql_gamma(em_list = .em)
    }
    
  }
  
  trc <- function(x){
    sum(diag(x))
  }
  
  
  log_normalising_penalty <- function(L, gam_mean, all_integrals = T){ # eta(beta) in paper
    #mvnorm_rs = matrix(rnorm(n * ncol(sigma)), nrow = n), fixed to avoid computational costs.
    U <- mvnorm_rs %*% t(L) # mv normal random sample with mean = 0 and var = H = L %*% t(L)
    
    re_oralce <- function(gam_re){
      
      subject_gam <- gam_mean + c(gam_re, rep(0, times = with(model$specs, p - r)))
      
      gcreg::is_monotone(p = mean_pl_cv$to_mono(subject_gam), region = model$rcontr_region)
      
    }
    
    in_region <- sum(apply(U, MARGIN = 1, re_oralce))
    
    if(in_region == 0) return(1e08)
    
    mult <- ifelse(all_integrals, model$specs$g, 1)
    
    
    return(
      mult * ( log(in_region) - log(nrow(U)) )
    )
    
  }
  
  ql_hood <- function(em_list){
    
    # residual sums of squares for gamma, after taking expectation
    RSS <- with(em_list, (crossprod(Y - Xo %*% gam - Zo %*% re_mean)[1,1] + trc(Zo %*% re_var %*% t(Zo))) / sig2)
    
    log_det_R <- model$specs$N * log(em_list$sig2)
    
    L <- make_L_mat(em_list$omeg)
    
    if(class(try( L_inv <- solve(L), silent=T)) != "matrix") return(1e08)
    
    G_L_inv <- bdiag(rep(list(L_inv), times = model$specs$g))
    
    G_L_inv_re_mean <- G_L_inv %*% em_list$re_mean
    
    # sum of squares for random effects, after taking expectation
    SSRE <- (crossprod(G_L_inv_re_mean) + trc(crossprod(G_L_inv) %*% em_list$re_var))[1,1]
    log_det_L <- model$specs$g * log(det(L)) # or sum(log(diag(L))) # multiply by 2 before returning below
    
    penalty <- log_normalising_penalty(L = L, gam_mean = em_list$gam) # implicitly uses mvnorm_rs
    
    RSS + SSRE + log_det_R + 2 * (log_det_L + penalty) # deviance scale
  }
  
  ql_hood_omega_wrapper <- function(omeg, em_list){
    
    .em <- em_list
    .em$omeg <- omeg
    
    ql_hood(em_list = .em)
    
  }
  
  ql_hood_gamma_wrapper <- function(gam, em_list){
    
    .em <- em_list
    .em$gam <- gam
    
    ql_hood(em_list = .em)
    
  }
  
  Dpenalty_omega <- function(omeg, gam_mean, boundary_rs, total_rs){
    
    re_oralce <- function(gam_re){
      
      subject_gam <- gam_mean + c(gam_re, rep(0, times = with(model$specs, p - r)))
      
      gcreg::is_monotone(p = mean_pl_cv$to_mono(subject_gam), region = model$rcontr_region)
      
    }
    
    penalty_fun <- function(omeg){
      L <- make_L_mat(omeg)
      U <- boundary_rs %*% t(L) # mv normal random sample with mean = 0 and var = H = L %*% t(L)
      in_region <- sum(apply(U, MARGIN = 1, re_oralce))
      in_region / total_rs
    } 
   
    return(
      numDeriv::grad(func = penalty_fun, x = omeg, method = "simple", method.args=list(eps=1e-02))
    )
    
  }
  
  gen_calc_Dql_omega <- function(em_list_init){

    # find random sample that are on boundary to speed up computation.
    L_init <- make_L_mat(em_list_init$omeg)
    U <- mvnorm_rs %*% t(L_init)
    subject_gams <- apply(U, MARGIN = 1, FUN = function(re){em_list_init$gam + c(re, rep(0, times = with(model$specs, p - r)))})
  
    U_lin_dist_to_boundary <- apply(subject_gams, MARGIN = 2, 
                                    FUN = function(re_gam){ linear_dist_to_monotone_boundary_em(re_gam)$dist} )
    
    sample_on_boundary <- order(U_lin_dist_to_boundary)[1:floor(nrow(U) * 0.01)] # use 1% closest to boundary.
    
    penalty_val <- exp(log_normalising_penalty(L = L_init, gam_mean = em_list_init$gam, all_integrals = F))
    
    function(em_list){
      L <- make_L_mat(em_list$omeg)
  
      if(class(try(solve(L),silent=T)) != "matrix") return(rep(1e08, times = length(em_list$omeg)))
  
      L_inv <- solve(L)
  
      G_L_inv <- bdiag(rep(list(L_inv), times = model$specs$g))
  
      H_inv <- crossprod(L_inv)
  
      G_inv <- bdiag(rep(list(H_inv), times = model$specs$g))
  
      DL_omega_diag_elements <- list()
  
      for(i in 1:model$specs$r){
        DL_omega_diag_elements[[i]] <- matrix(0, ncol = model$specs$r, nrow = model$specs$r)
        DL_omega_diag_elements[[i]][i,i] <- L[i,i]
      }
  
      DL_omega_upper_elements <- list()
  
      for(i in 1:(model$specs$r-1)){
        for(j in (i+1):model$specs$r){
          dL <- matrix(0, nrow = model$specs$r, ncol = model$specs$r)
          dL[i,j] <- 1
          DL_omega_upper_elements[[length(DL_omega_upper_elements)+1]] <- dL
        }
      }
  
      DG_L_omega <- c(DL_omega_diag_elements, DL_omega_upper_elements) %>% lapply(FUN = function(dL){
        bdiag(rep(list(dL),times = model$specs$g))})
  
      Dql_omega_d <- sapply(DG_L_omega, FUN = function(dGL){
        var_temp <- G_inv %*% dGL %*% G_L_inv
        dlomega <- 2 * (trc(G_L_inv %*% dGL) - trc(var_temp %*% em_list$re_var) - t(em_list$re_mean) %*% var_temp %*% em_list$re_mean)
        dlomega[1,1]
      })
  
      # penalty term derivative
      if(with(model$specs, r > 2 & constr_r)){
        
       Dpen <- Dpenalty_omega(omeg = em_list$omeg, 
                       gam_mean = em_list$gam,
                       boundary_rs = mvnorm_rs[sample_on_boundary,], 
                       total_rs = nrow(mvnorm_rs)
                         )
        Dlog_pen <- model$specs$g * Dpen / penalty_val
  
      }
      
      return(matrix(Dql_omega_d + 2 * Dlog_pen, ncol = 1))
    }

  }
  
  gen_calc_Dql_omega_wrapper <- function(em_list){

    calc_Dql_omega <- gen_calc_Dql_omega(em_list)
    .em <- em_list
    
    function(omeg, ...){
      .em$omeg <- omeg
      calc_Dql_omega(em_list = .em)
    }

  }
  
  bounce_monotone_em <- function(em_list, aim_gam, maxit_bounces = 2, maxit_bumps = 5, tol = 1e-02){
    
    Dpl_base <- deriv(pl_base)
    
    old_ql_hood <- -Inf
    new_ql_hood <- 0
    cur_em <- em_list
    
    bb <- 1
    
    #while(abs(new_ql_hood - old_ql_hood) > tol & bb <= maxit_bounces){
    while(new_ql_hood > old_ql_hood & bb <= maxit_bounces){
      
      old_ql_hood <- new_ql_hood <- ql_hood(cur_em)
      
      loc <- linear_dist_to_monotone_boundary_em(gam = cur_em$gam)$loc
      bump_dir <- sapply(Dpl_base, predict, newdata = loc)
      bump_norm <- bump_dir / sqrt(sum(bump_dir^2))
      max_bump <- bump_norm
      new_start <- cur_em$gam + max_bump
      
      while(!monotone_oracle(new_start)){
        max_bump <- max_bump / 2
        new_start <- cur_em$gam + max_bump
        new_start_neg <- cur_em$gam - max_bump
        if(monotone_oracle(new_start_neg) & !monotone_oracle(new_start)){
          new_start <- new_start_neg
          max_bump <- - max_bump
        }
      }
      
      bump <- max_bump / 2
      
      for(kk in 1:maxit_bumps){
        
        temp_em <- cur_em
        temp_em$gam <- optim_cols(par = new_start, aim_gamma = aim_gam, Y = Y, X = Xo, oracle_fun = monotone_oracle)
        
        temp_ql_hood <- ql_hood(temp_em)
        
        if(temp_ql_hood < old_ql_hood) {
          cur_em <-  temp_em
          new_ql_hood <- temp_ql_hood
          break
        }
        
        bump <- bump / (2^(kk))
        new_start <- cur_em$gam + bump
        
      }
      
      bb <- bb + 1
    }
    #cat(cur_em$gam,"\n")
    return(cur_em)
    
  }
  
  ####
  #return(mget(x = ls()))
  # Algorithm storage, initialisation
  
  if(save_steps) em <- list()
  
  # Some choice for initial par in monomial basis
  if(is.null(start$beta)){
    fe_fit <- cpm(formula = y ~ x, data = model$dat, degree = model$specs$p - 1, constraint = "monotone", c_region = model$mcontr_region)
    beta_init <- coef(fe_fit)
  } else {
    beta_init <- start$beta
  }
  
  if(!is_monotone(beta_init, region = model$mcontr_region)) stop("Initial parameter chosen is not monotone over required region")
  
  gamma_init <- mean_pl_cv$to_ortho(beta_init)
  
  # gam = estimate mean polynomial in orthonormal basis
  # sig2 = variance of the error
  # omeg = paramters for log-Cholesky decomposition of the random effects var-covar matrix
  
  if(is.null(start$sig2) | length(start$sig2) != 1){
    sig2_init <- 0.5
  } else {
    sig2_init <- start$sig2
  }
  
  omeg_length <- with(model$specs, r*(r+1)/2)
  
  if(is.null(start$omeg) | length(start$omeg) != omeg_length){
    omeg_init <- rep(0.5, times = omeg_length)
  } else {
    omeg_init <- start$omeg
  }
  
  prev_em <- list(gam = gamma_init,
                  sig2 = sig2_init,
                  omeg = omeg_init
  )
  
  if(save_steps) {
    em[[1]] <- prev_em 
  } else {
    em <- F
  }
  
  # EM Algorithm Components
  
  param_convg_check <- function(em1, em2, tol){
    nms <- c("gam","sig2","omeg")
    
    for(n in nms){
      if(!all(abs((em2[[n]] - em1[[n]])/(em1[[n]] + 0.01)) < tol)) return(F)
    }
    
    return(T)
  }
  
  # Run EM Algorithm
  
  # start vals
  convg <- F 
  convg_status <- paste("did not converge in",maxit, "iterations")
  which_min <- 0
  i <- 1
  
  ql_hood_vals <- numeric(length = maxit)
  first_e_step <- prev_em %>% E_step()
  ql_hood_vals[1] <- first_e_step %>% ql_hood()
  
  if(verbose) cat("it:", 1, "ql:", ql_hood_vals[1], "\n")
  
  while(!convg & i < maxit){
    
    if(i == 1){
      new_em <- first_e_step %>% M_step()
    } else {
      new_em <- prev_em %>% E_step() %>% M_step()
    }
    
    
    ql_hood_vals[i+1] <- new_em %>% ql_hood()
    
    
    # check monotonicity (keeping for testing)
    if(!monotone_oracle(new_em$gam)) stop("Gamma not monotone: it = ", i + 1)
    
    # check convergence
    
    if(param_convg_check(prev_em, new_em, tol = tol)){
      
      convg <- T
      convg_status <- paste("Parameters converged with tol:", tol)
      ql_hood_vals <- ql_hood_vals[1:(i + 1)]
      which_min <- i + 1
      res <- new_em
      
    } else if(check_up & i > 1){
      
      #ql_hood_check <- prev_em %>% E_step() %>% ql_hood()
      ql_hood_check <- ql_hood_vals[i]
      
      if(ql_hood_vals[i+1] > ql_hood_check){ # on divergence scale so should be decreasing
        
        #if(verbose) cat("quasi-log-likelihood did not decrease \n")
        convg <- T
        convg_status <- paste("quasi-log-likelihood did not decrease")
        ql_hood_vals <- ql_hood_vals[1:(i+1)]
        
        best_models <- list(prev2_em, prev_em, new_em) %>% 
          lapply(function(e){e %>% E_step()})
        
        which_min <- best_models %>% sapply(ql_hood) %>% which.min()
        
        res <- best_models[[which_min]]
        
      }
    } 
    
    if(save_steps) em[[i + 1]] <- new_em
    
    i <- i + 1
    
    prev_em <- new_em
    prev2_em <- prev_em
    
    
    
    if(verbose) cat("it:", i, "ql:", ql_hood_vals[i], "\n")
    
  } 
  
  if(maxit <= i){
    res <- new_em %>% E_step()
  }
  
  beta_mean <- mean_pl_cv$to_mono(res$gam)
  
  beta_grp <- lapply(seq(0, nrow(res$re_mean)-1, by = model$specs$r), FUN = 
                       function(i){
                         beta_mean + 
                           c(grp_pl_cv$to_mono(res$re_mean[i + 1:model$specs$r,]), 
                             rep(0, times = model$specs$p - model$specs$r)
                           )
                       }
  )
  
  names(beta_grp) <-  model$group_ids
  
  return(list(beta_mean = beta_mean,
              beta_grp = beta_grp,
              H_mat = tcrossprod(make_L_mat(res$omeg)),
              sig2 = res$sig2,
              res = res,
              em_list = em,
              ql_vals = ql_hood_vals,
              which_min = which_min,
              converged = convg,
              status = convg_status,
              funcs = list(
                ql_hood = ql_hood,
                mean_pl_cv = mean_pl_cv,
                grp_pl_cv = grp_pl_cv,
                monotone_oracle = monotone_oracle,
                E_step = E_step,
                M_step = M_step
              ),
              model = model
  )
  )
  
}
