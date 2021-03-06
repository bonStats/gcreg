% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixed-effect-em-constrained-polynomial-model.R
\name{constrained_lmm_em}
\alias{constrained_lmm_em}
\title{Fit monotone-constrained polynomial LMM with an EM algorithm}
\usage{
constrained_lmm_em(model, maxit = 200, tol, start = NULL, verbose = F,
  save_steps = F, check_up = F)
}
\arguments{
\item{model}{model specification made by \code{\link{make_em_model_specs}}.}

\item{maxit}{maximum number of iterations.}

\item{tol}{relative tolerance for convergence of parameter values.}

\item{start}{start values for parameters. Defaults to NULL.}

\item{verbose}{output current iteration and quasi-log-likelihood value.}

\item{save_steps}{should each em step be saved and return in a list?}

\item{check_up}{should the algorithm check if the quasi-log-likelihood is improving each iteration?}
}
\value{
list including estimated mixed effects model.
}
\description{
EXPIREMENTAL: Fit polynomial constrained mixed effect model with monotonic mean and (potentially) monotonic random effects (\eqn{r <= 2} for monotone random effects). 
See \code{\link{constrained_lmm_mcem}} for (\eqn{r > 2}).
For use with block-random effects (non-crossed), e.g. subject-specific random effects.
}
\details{
The interface and default actions of this function are under development and may change without warning. See the mixed effects vignette for usage details.
}
\keyword{internal}
