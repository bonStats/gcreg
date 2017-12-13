#' gcreg: A package for fitting constrained regression models
#' 
#' Currently under development. See \url{https://github.com/bonStats/gcreg} for updates.
#' 
#' Function \code{\link{cpm}} fits constrained polynomial models with fixed effects.
#' 
#' Hidden functions \code{\link{constrained_lmm_em}} and \code{\link{constrained_lmm_mcem}} fit constrained polynomial models with mixed effects.
#' 
#' Hidden function \code{\link{case_bootstrap_constrained_lmm_em}} implements case bootstrapping for standard errors for constrained mixed effects models.
#' 
#' Access hidden functions with \code{gcreg:::} prefix, e.g. \code{?gcreg:::constrained_lmm_em}.
#'
#' @docType package
#' @name gcreg
#' @import polynom
#' @import Matrix
#' @import dplyr
#' @importFrom mvtnorm rmvnorm dmvnorm TVPACK
#' @importFrom tmvtnorm mtmvnorm
#' @importFrom numDeriv grad
#' @importFrom parallel mclapply
#' @importFrom stats is.empty.model model.matrix model.weights model.response predict var deriv formula napredict naresid coef fitted pnorm dnorm setNames optim rnorm runif
NULL