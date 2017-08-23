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
#' @param start intial value of COLS optimisation
#' @param ... arguments to be passed to control_cols()
#' @export

cpm <- function(formula, data, subset, weights, na.action,
                degree, constraint = NULL, oracle = NULL, start, ...) {
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
  
  if (!missing(constraint)){
    if(!missing(oracle)) warning("Only specify one of 'constraint' and 'oracle'. Using 'constraint', 'oracle' argument ignored.")
   
    constraint_name <- match.arg(constraint, c("monotone","convex"))
   
    if(constraint_name == "convex"){
     stop("Convexity not yet implemented")
    } else {
     oracle <- is_monotone
     if(missing(start)){
       init_par <- rep(c(0.1,1), times = floor((degree+1)/2))
     } else {
       init_par <- start
     }
    }
   
  } else {
    # use oracle
    if(missing(init_par)) stop("Please provide initial parameter when using the oracle = argument.")
  }
  
  # check starting value is within boundary
  if(!oracle(init_par)) stop(paste("Initial value outside constrained space, check starting value and degree =", degree, "is appropriate."))
  if(length(init_par) != degree + 1) stop("length of start/init_par is incorrect for specified degree")
  
  poly_basis <- make_disc_orthonormal_basis(x, deg = degree)
  
  # orthonormal design matrix
  Xo <- sapply(poly_basis, predict, newdata = x) # polynom:::predict.polynomial
  
  # basis converter
  cv <- gen_poly_basis_converters(poly_basis)
  
  z <- optim_cols(par = cv$to_ortho(init_par), Y = y, X = Xo, oracle_fun = oracle, control = cols_control(method = "best-step", ...))
  
  beta_par <- as.numeric(cv$to_mono(z))
  names(beta_par) <- c("Intercept", paste0("x^",1:degree))
  
  RSS <- sum((y - Xo%*%z)^2)
  
  return(list(beta = beta_par, RSS = RSS))
  
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


