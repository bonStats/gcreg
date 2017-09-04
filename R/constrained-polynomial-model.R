#' WIP: Determine the least squares estimates of a constrained polynomial regression
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which the function is called.
#' @param subset NOT IMPLEMENTED an optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights NOT IMPLEMENTED
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param degree degree of polynomial to be fit
#' @param constraint an optional character string of "monotone" or "convex"
#' @param oracle an optional function of class "oracle", returning TRUE when a given point is inside the constrained set and FALSE otherwise.
#' @param start intial value of COLS optimisation.
#' @param c_region the applicable region for the constraint, default is (-Inf,Inf), the real line.
#' @param ... arguments to be passed to control_cols()
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
  
  if (!missing(constraint)){ 
    if(!missing(oracle)) warning("Only specify one of 'constraint' and 'oracle'. Using 'constraint', 'oracle' argument ignored.")
    
    constraint_name <- match.arg(constraint, c("monotone","convex"))
   
    if(constraint_name == "convex"){
     stop("Convexity not yet implemented")
    } else {
      # convert boundary to scale data
      sc_oracle_region <- sc_fun_x$scale(c_region)
      # define oracle on scaled boundaries
      sc_oracle <- function(p) is_monotone(cv$to_mono(p), region = sc_oracle_region)
     if(missing(start)){
       # assumes monotone increasing
       init_par <- rep(c(0.1,1), times = floor((degree+1)/2))
     } else {
       init_par <- start
     }
      init_par_in <- sc_oracle(cv$to_ortho(init_par))
    }
   
  } else {
    
    # use oracle argument
    if(missing(start)) stop("Please provide initial parameter when using the oracle = argument.")
    init_par_in <- oracle(start)
    oracle_region <- attr(init_par_in, "region")
    init_par <- start
    if(is.null(oracle_region)) {
      warning("No boundaries specified in attributes of return of oracle function, assuming constraining region is: (", paste(c_region,collapse = ","),")")
      oracle_region <- c_region
    } 
    
    region_arg_exists <- !is.null(formals(oracle)$region)
    
    if(region_arg_exists){
      sc_oracle_region <- sc_fun_x$scale(c_region)
      sc_oracle <- function(p) {oracle(cv$to_mono(p), region = sc_oracle_region)}
    } else {
      if(any(is.finite(c_region))) warning("Applicable region for oracle function not scaled for finite boundary. Please provide oracle function with argument 'region'")
      sc_oracle <- function(p) {oracle(cv$to_mono(p))}
    }
    
  }
  
  # check starting value is within boundary
  if(!init_par_in) stop(paste("Initial value outside constrained space, check starting value and degree =", degree, "is appropriate."))
  
  if(length(init_par) != degree + 1) stop("length of start/init_par is incorrect for specified degree")
  
  # need to scale init param also
  z <- optim_cols(par = cv$to_ortho(init_par), Y = y_sc, X = Xo, oracle_fun = sc_oracle, control = cols_control(method = "best-step", ...))
  
  # un-orthonormalise (discrete polynomial orthogonalisation)
  sc_xy_beta_par <- as.numeric(cv$to_mono(z))
  
  # unscale y
  sc_x_beta_par <- replace(sc_xy_beta_par, 1, sc_xy_beta_par[1] - attr(sc_fun_y,"c")) / attr(sc_fun_y,"b")
  
  # un-scale (range of x data was set to between -1 and 1).
  b_unsc <- attr(sc_fun_x,"b") ^ (0:degree)
  sc_x_origin_pl <- polynomial(sc_x_beta_par  * b_unsc)
  
  beta_pl <- change.origin(sc_x_origin_pl, o = attr(sc_fun_x,"c")/attr(sc_fun_x,"b"))
  beta_par <- as.numeric(beta_pl)
  
  names(beta_par) <- c("Intercept", paste0("x^",1:degree))
  
  y_predict <- predict(beta_pl, newdata = x)
  
  RSS <- sum((y - y_predict)^2)
  
  
  
  return(
    list(beta = beta_par, RSS = RSS, 
         scaled = list(beta = sc_xy_beta_par, oracle = sc_oracle, x = x_sc, y = y_sc)
              )
    )
  
  # z$na.action <- attr(mf, "na.action")
  # z$call <- cl
  # z$terms <- mt
  # if (model) 
  #   z$model <- mf
  # if (ret.x) 
  #   z$x <- x
  # if (ret.y) 
  #   z$y <- y
  # structure(z, class = "polreg")
  
  
}


