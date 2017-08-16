#' WIP: Determine the least squares estimates of a constrained polynomial regression
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which the function is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights NOT IMPLEMENTED
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param degree degree of polynomial to be fit
#' @param constraint an optional character string of "monotone" or "convex"
#' @param oracle an optional function of class "oracle", returning TRUE when a given point is inside the constrained set and FALSE otherwise.
#' @param ... arguments to be passed to control_cols()
#' @export

cpm <- function(formula, data, subset, weights, na.action,
                degree, constraint = NULL, oracle = NULL, ...) {
  
  m <- match(c("formula", "data", "subset", "weights", "na.action"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  
  if (is.empty.model(mt)) 
    stop("You should specify a regressor variable.")
  
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
     oracle <- is_monotone()
    }
   
  }
  
  poly_basis <- make_disc_orthonormal_basis(x, deg = degree)
  
  Xo <- sapply(poly_basis, predict, newdata = x) # polynom:::predict.polynomial
  
  # grab or set initial value par
  
  
  z <- optim_COLS(par, Y = y, X = Xo, oracle_fun = oracle, control = cols_control(...))
  
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


