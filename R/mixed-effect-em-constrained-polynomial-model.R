#' Fit polynomial constrained mixed effect model with monotonic mean and (potentially) monotonic random effects.
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
#' @keywords internal
#' @return list including estimated mixed effects model.

constrained_lmm_em <- function(model, maxit = 200, tol, start = NULL, verbose = F, save_steps = F, check_up = F){
  
  if(missing(model)) stop("model argument must be specified. Use gcreg:::make_em_model_specs().")
  
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
    
    if(with(model$specs, r == 1 | !constr_r)){
      
      .em$re_mean <- M_U_uc
      .em$re_var <- V_U_uc
      
    } else if(model$specs$r == 2){ # constrained and r==2
      
      ## constrained mean/variance matrices of random effects
      # Also refered to as -c(beta). Calculates linear distance from monotonic boundary for r = 2 to be used in truncation of R.E. 
      random_slope_trunc <- - linear_dist_to_monotone_boundary_em(gam = .em$gam)$dist
      # truncation vector for random intercept and slope
      lower_trunc <- rep(c(-Inf, random_slope_trunc), times = model$specs$g)
      # mean and variance from mtmvnorm (wrapper function)
      trunc_MV <- mean_var_re_constrained2(mean_uncon = as.vector(M_U_uc), var_uncon = V_U_uc, lower_trunc = lower_trunc, em_list = .em)
      
      # constrained mean and variance matrices
      .em$re_mean <- trunc_MV$tmean
      .em$re_var  <- trunc_MV$tvar
      
    } else if(model$specs$r == 3){
      stop("r = 3 and constrained random effects not yet implemented")
    } else {
      stop("r > 3 and constrained random effects not implemented")
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
  
  mean_var_re_constrained2 <- function(mean_uncon, var_uncon, lower_trunc, em_list){
    
    # store components of bdaig matrix
    trunc_moments_list <- list()
    
    # find any truncated moments that calculated badly.
    bad_r <- NULL
    
    for(i in 1:model$specs$g){
      elems <- (2 * i - 1):(2 * i) # elements to access
      trunc_moments <- mtmvnorm(mean = mean_uncon[elems], 
                                sigma = var_uncon[elems,elems], 
                                lower = lower_trunc[elems],
                                doComputeVariance = T,
                                pmvnorm.algorithm = TVPACK() # can try different algorithms
      )
      
      # check if returned NaN, use previous mean/var instead
      if(any(is.nan(trunc_moments$tmean)) | any(is.nan(trunc_moments$tvar))){
        if(is.null(em_list$re_mean) &  !is.null(em_list$re_var)){
          trunc_moments_list[[i]] <- with(em_list, list(tmean = re_mean[elems], tvar = re_var[elems, elems]))
        } else {
          trunc_moments_list[[i]] <- list(tmean = mean_uncon[elems], tvar = var_uncon[elems, elems])
        }
        bad_r <- c(bad_r,i)
      } else {
        trunc_moments_list[[i]] <- trunc_moments
      }
    }
    
    if(!is.null(bad_r)){
      warning("moments of random effects for subject",paste(bad_r, collapse = ", "),
              " unable to be truncated, using previous moments this iteration instead",
              immediate. = T)
    }
    
    return(
      list(tmean = lapply(trunc_moments_list, function(x){x$tmean}) %>% c(recursive = T) %>% matrix(ncol = 1),
           tvar =  lapply(trunc_moments_list, function(x){x$tvar}) %>% bdiag()
      )
    )
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
  
  ####
  
  #### M-step functions ####
  
  M_step <- function(em_list, step_size = 0.1){ #remove step size ???
    
    .em <- em_list
    
    par_st <- gamma_init # if on boundary need to do more work.????
    
    aim_gamma <- optim(par = em_list$gam,
                       fn = ql_hood_gamma_wrapper,
                       gr = calc_Dql_gamma_wrapper, 
                       em_list = em_list,
                       method = "BFGS",
                       control = list(maxit = 10)
    )$par
    
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
      
      opt_omegas <- optim(par = .em$omeg,
                          fn = ql_hood_omega_wrapper,
                          gr = calc_Dql_omega_wrapper,
                          em_list = .em,
                          method = "L-BFGS-B",
                          lower = with(model$specs, c(rep(-10, times = r), rep(-20, times = (r - 1) * r / 2))),
                          upper = with(model$specs, c(rep( 10, times = r), rep( 20, times = (r - 1) * r / 2)))
                          # lower/upper bounds stop numerical difficulty with solve(L). Should generically decide on these based on problem size.
      )

    }
    
    .em$omeg <- opt_omegas$par

    return(.em)
  }
  
  
  # helper functions - M step
  
  # quasi-likelihood derivative wrt gamma
  calc_Dql_gamma <- function(em_list){
    
    Drss <- with(em_list, - 2 * (Xo_t_Y - gam - Xo_t_Zo %*% re_mean) / sig2)
    
    if(with(model$specs, r == 2 & constr_r)){
      c_beta <- linear_dist_to_monotone_boundary_em(gam = em_list$gam) # boundary for slope RE
      Dc_beta <- calc_Dc_beta(gam = em_list$gam, min_x = c_beta$loc)
      var_random_slope <- exp(em_list$omeg[2])
      z_scr <- - c_beta$dist[1] / var_random_slope
      Dpenalty <- 2 * model$specs$g * dnorm(z_scr) * Dc_beta / (var_random_slope * (1 - pnorm(z_scr)))
    } else {
      c_beta <- 0
      Dc_beta <- 0
      Dpenalty  <- 0
      if(with(model$specs, r > 2 & constr_r)) warning("Derivative of quasi-likelihood not defined for constrained random effects and r > 2 ")
    }
    
    return(list(drv = matrix(Drss + Dpenalty),
                Drss = matrix(Drss),
                Dpenalty = matrix(Dpenalty),
                c_beta = c_beta,
                Dc_beta = Dc_beta
    )
    )
    
  }
  
  calc_Dc_beta <- function(gam, min_x){ # using evelope theorem
    
    sapply(Dpl_base, predict, newdata = min_x)
    
  }
  
  calc_Dql_gamma_wrapper <- function(gam, em_list){
    
    .em <- em_list
    .em$gam <- gam
    
    calc_Dql_gamma(em_list = .em)$drv
    
  }
  
  trc <- function(x){
    sum(diag(x))
  }
  
  ql_hood <- function(em_list){
    
    # residual sums of squares for gamma, after taking expectation
    RSS <- with(em_list, (crossprod(Y - Xo %*% gam - Zo %*% re_mean)[1,1] + trc(Zo %*% re_var %*% t(Zo))) / sig2)
    
    log_det_R <- model$specs$N * log(em_list$sig2)
    
    L <- make_L_mat(em_list$omeg)
    
    if(class(try(solve(L),silent=T)) != "matrix") return(1e08)
    
    L_inv <- solve(L) 
    
    G_L_inv <- bdiag(rep(list(L_inv), times = model$specs$g))
    
    G_L_inv_re_mean <- G_L_inv %*% em_list$re_mean
    
    
    # sum of squares for random effects, after taking expectation
    SSRE <- (crossprod(G_L_inv_re_mean) + trc(crossprod(G_L_inv) %*% em_list$re_var))[1,1]
    log_det_L <- model$specs$g * log(det(L))
    
    # only works for r = 2
    if(with(model$specs, r == 2 & constr_r)){
      c_beta <- linear_dist_to_monotone_boundary_em(gam = em_list$gam) #boundary for re slope
      penalty <- model$specs$g * log(1 - pnorm(- c_beta$dist / exp(em_list$omeg[2]))) # penalty term # lower.tail = FALSE if monotonically decreasing
    } else {
      penalty <- 0
    }
    
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
  
  calc_Dql_omega <- function(em_list){
    
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
    if(with(model$specs, r == 2 & constr_r)){
      c_beta <- linear_dist_to_monotone_boundary_em(gam = em_list$gam)  # boundary for slope RE
      var_random_slope <- exp(em_list$omeg[2])
      z_scr <- - c_beta$dist[1] / var_random_slope
      Dpenalty <- - 2 * model$specs$g * dnorm(z_scr) * c_beta$dist[1] / ((var_random_slope) * (1 - pnorm(z_scr)))
      Dql_omega_d[2] <- Dql_omega_d[2] + Dpenalty
      
    } 
    
    return(matrix(Dql_omega_d, ncol = 1))
    
  }
  
  calc_Dql_omega_wrapper <- function(omeg, em_list){
    
    .em <- em_list
    .em$omeg <- omeg
    
    calc_Dql_omega(em_list = .em)
    
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
  ql_hood_vals[1] <- prev_em %>% E_step() %>% ql_hood()
  
  if(verbose) cat("it:", 1, "ql:", ql_hood_vals[1], "\n")
  
  while(!convg & i < maxit){
    
    
    new_em <- prev_em %>% E_step() %>% M_step()
    
    
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
      
      ql_hood_check <- prev_em %>% E_step() %>% ql_hood()
      #ql_hood_check <- ql_hood_vals[i]
      
      if(ql_hood_vals[i+1] > ql_hood_check){ # on divergence scale so should be decreasing
        
        if(verbose) cat("quasi-log-likelihood did not decrease \n")
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
              sig = res$sig2,
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
