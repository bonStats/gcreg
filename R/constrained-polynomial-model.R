#' Determine the least squares estimates of a constrained polynomial regression
#' 
#' Estimates univariate polynomials with parameters from a closed convex set. 
#' Some examples are shape constraints (e.g. montonicity, convexity) or simply linear parameter constraints
#' such as \eqn{\beta_{1} >= 0}{\beta [1] >= 0}.
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which the function is called.
#' @param subset NOT YET IMPLEMENTED an optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights NOT YET IMPLEMENTED
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param degree degree of polynomial to be fit
#' @param constraint an optional character string of "monotone" ("convex" coming soon)
#' @param oracle an optional function of class "oracle", returning TRUE when a given point is inside the constrained set and FALSE otherwise.
#' @param start intial value of COLS optimisation.
#' @param c_region the applicable region for the constraint, default is (-Inf,Inf), the real line.
#' @param ... arguments to be passed to control_cols()
#' @return Constrained regression model with class \code{creg}
#' @export

cpm <- function(formula, data, subset, weights, na.action,
                degree, constraint = NULL, oracle = NULL, start, c_region = c(-Inf,Inf), ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(c("formula", "data", "subset", "weights", "na.action"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  if (is.empty.model(mt)) 
    stop("You should specify a regressor variable.")
  NoIntercept <- attr(mt, "intercept") == 0
  attr(mt, "intercept") <- 0
  x <- model.matrix(mt, mf)
  if (NCOL(x) != 1) 
    stop("Regressor variable should be univariate.")
  
  if (missing(degree)) 
    stop("'degree' should be specified") # to be updated
  if (trunc(degree) != degree || degree <= 0 || degree >= length(x)) 
    stop("'degree' should be a positive integer less than length(x).")
  
  # scale_data
  sc_fun_x <- gen_scale_data_funs(x)
  x_sc <- sc_fun_x$scale(x)
  
  sc_fun_y <- gen_scale_data_funs(y)
  y_sc <- sc_fun_y$scale(y)
  
  # orthonormal polynomial basis
  poly_basis <- make_disc_orthonormal_basis(x_sc, deg = degree)
  
  # orthonormal design matrix
  Xo <- sapply(poly_basis, predict, newdata = x_sc) # polynom:::predict.polynomial
  
  # basis converter
  cv <- gen_poly_basis_converters(poly_basis)
  
  # removes orthonormality and data transformation 
  untransform_gam <- function(gam){
    
    # un-orthonormalise (discrete polynomial orthogonalisation)
    sc_xy_beta_par <- as.numeric(cv$to_mono(gam))
    
    len <- length(sc_xy_beta_par)
    
    # unscale y
    sc_x_beta_par <- replace(sc_xy_beta_par, 1, sc_xy_beta_par[1] - attr(sc_fun_y,"c")) / attr(sc_fun_y,"b")
    
    # un-scale (range of x data was set to between -1 and 1).
    b_unsc <- attr(sc_fun_x,"b") ^ (0:degree)
    sc_x_origin_pl <- polynomial(sc_x_beta_par * b_unsc)
    
    beta_pl <- change.origin(sc_x_origin_pl, o = attr(sc_fun_x,"c")/attr(sc_fun_x,"b"))
    
    beta_num <- as.numeric(beta_pl)
    
    beta_par <- rep(0, times = len)
    beta_par[1:length(beta_num)] <- beta_num
    
    attr(beta_par, "beta_pl") <- beta_pl
    
    return(beta_par)

  }
  
  transform_beta <- function(beta){
    
    beta_par <- as.numeric(beta)
    
    len <- length(beta_par)
    
    #scale x (range of x data was set to between -1 and 1)
    b_unsc <- attr(sc_fun_x,"b") ^ (0:degree)
    sc_x_beta_par <- beta_par / b_unsc
    
    sc_x_beta <- as.numeric(
      change.origin(polynomial(sc_x_beta_par), o = - attr(sc_fun_x,"c"))
    )
    
    sc_x_beta_par <- rep(0, times = len)
    sc_x_beta_par[1:length(sc_x_beta)] <- sc_x_beta
    
    # scale y
    sc_xy_beta_par <- as.numeric(sc_x_beta_par) * attr(sc_fun_y,"b")
    sc_xy_beta_par <- replace(sc_xy_beta_par, 1, sc_xy_beta_par[1] + attr(sc_fun_y,"c"))
    
    # un-orthonormalise (discrete polynomial orthogonalisation)
    gam <- as.numeric(cv$to_ortho(sc_xy_beta_par))
    
    return(gam)
    
  }
  
  if (!missing(constraint)){ 
    if(!missing(oracle)) warning("Only specify one of 'constraint' and 'oracle'. Using 'constraint', 'oracle' argument ignored.")
    
    constraint_name <- match.arg(constraint, c("monotone","convex"))
   
    if(constraint_name == "convex"){
     stop("Convexity not yet implemented")
    } else {
      # convert boundary to scale data
      sc_oracle_region <- sc_fun_x$scale(c_region)
      # define oracle on scaled boundaries
        # use cv rather than untransform_gam since monotonicity unaffected by data transformations
      sc_oracle <- function(p) is_monotone(cv$to_mono(p), region = sc_oracle_region)
     if(missing(start)){
       # assumes monotone increasing
       if(degree+1 %% 2 == 0) {
         init_par <- cv$to_ortho(
           rep(c(0.1,1), times = (degree+1)/2)
         )
       } else {
         init_par <- cv$to_ortho(
          c(rep(c(0.1,1), times = degree/2),0.1)
         )
       }
     } else {
       init_par <- cv$to_ortho(start)
     }
      init_par_in <- sc_oracle(init_par)
    }
   
  } else {
    
    # use oracle argument
    if(missing(start)) stop("Please provide initial parameter when using the oracle = argument.")
    init_par_in <- oracle(start)
    oracle_region <- attr(init_par_in, "region")
    init_par <- transform_beta(start)
    
    if(is.null(oracle_region)) {
      warning("No boundaries specified in attributes of return of oracle function, assuming constraining region is: (", paste(c_region,collapse = ","),")")
      oracle_region <- c_region
    } 
    
    if(!is.null(oracle_region) & !missing(c_region)){
      warning("Oracle function specifies region and c_region given. Defaulting to oracle's region")
    }
    
    region_arg_exists <- !is.null(formals(oracle)$region)
    
    if(region_arg_exists & !missing(c_region)){
      sc_oracle <- function(p) {oracle(untransform_gam(p), region = c_region)}
    } else {
      sc_oracle <- function(p) {oracle(untransform_gam(p))}
    }
    
  }
  
  # check starting value is within boundary
  if(!init_par_in) stop(paste("Initial value outside constrained space, check starting value and degree =", degree, "is appropriate."))
  
  if(length(init_par) != degree + 1) stop("length of start/init_par is incorrect for specified degree")
  
  ctrl_list <- cols_control(method = "best-step", ...)
  
  # need to scale init param also
  z <- optim_cols(par = init_par, Y = y_sc, X = Xo, oracle_fun = sc_oracle, control = ctrl_list)
  
  #bounces to deal with flat spots if on monotone boundary (others to be implemented)
  if(!missing(constraint)){
    if(constraint == "monotone" & ctrl_list$maxit_bounces > 0){
      boundary_distance <- linear_dist_to_monotone_boundary(gam = z, region = c_region, basis_cv = cv)$dist
      
      if(abs(boundary_distance) < 1e-04){
        z <- bounce_monotone(gam = z, Y = y_sc, Xo = Xo, poly_basis = poly_basis, oracle_fun = sc_oracle, region = c_region, basis_cv = cv, control = ctrl_list)
      }
    }
  }
 
  beta_par <- untransform_gam(gam = z)
  beta_pl <- attr(beta_par, "beta_pl")
  attributes(beta_par) <- NULL
  
  names(beta_par) <- c("Intercept", paste0("x^",1:degree))
  
  z <- list(beta = beta_par)
  
  z$fitted.values <- predict(beta_pl, newdata = x)
  attributes(z$fitted.values) <- attributes(y)
  
  z$residuals <- y - z$fitted.values
  
  z$RSS <- sum(z$residuals^2)
  
  z$controls <- ctrl_list
  
  z$na.action <- attr(mf, "na.action")
  z$call <- cl
  z$terms <- mt
  # if (model)
  #   z$model <- mf
  # if (ret.x)
  #   z$x <- x
  # if (ret.y)
  #   z$y <- y
  structure(z, class = "creg")
  
  
}


