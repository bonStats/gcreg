#' Methods for class "creg"
#' 
#' @param object The constrained model
#' @param ... Additional arguments to be passed to function
#' @export
#' @method coef creg
coef.creg <- function(object, ...){
  
  object$beta
  
}

#' Regression coefficients
#' 
#' @param x The constrained model
#' @param ... Additional arguments to be passed to function
#' @export
#' @method formula creg
formula.creg <- function(x, ...){
  
  formula(x$terms)
  
}

#' Fitted values
#' 
#' @param object The constrained model
#' @param ... Additional arguments to be passed to function
#' @export
#' @method fitted creg
fitted.creg <- function(object, ...){
  
  napredict(object$na.action, object$fitted.value)
  
}

#' Regression residuals
#' 
#' @param object The constrained model
#' @param ... Additional arguments to be passed to function
#' @export
#' @method residuals creg
residuals.creg <- function(object, ...){
  
  naresid(object$na.action, object$residuals)
  
}

#' Predicted (or fitted) values
#' 
#' @param object The constrained model
#' @param ... Additional arguments to be passed to function
#' @param newdata Data to be provided to predict()
#' @export
#' @method predict creg
predict.creg <- function(object, newdata, ...){
  
  if(missing(newdata)){ return(fitted(object)) }
  
  newdata <- as.data.frame(newdata)
  
  pred <- predict(polynomial(coef(object)), newdata = newdata)
  
  pred
  
}

#' Print regression model
#' @param x The constrained model
#' @param ... Additional arguments to be passed to function
#' @param digits Minimal number of significant digits for print()
#' @export
#' @method print creg
print.creg <- function(x, digits = max(3, getOption("digits") - 3), ...){
  
  cat("\nConstrained regression model\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cff <- coef(x, ...)
  if (length(cff)) {
    cat("Coefficients:\n")
    print.default(format(cff, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
  
}