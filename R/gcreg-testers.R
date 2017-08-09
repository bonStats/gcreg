#### Test Data ###
library(polynom)
library(MonoPoly)
library(dplyr)

source("R/discrete-orthonormal-polynomial-basis-generator.R")
source("R/line-search.R")
source("R/is-monotone.R")
source("R/polynomial-basis-converter.R")
source("R/constrained-orthonormal-least-squares.R")

N <- 1000
q <- 5
x <- runif(N,min = -1,max = 1)
poly_orth <- make_disc_orthonormal_basis(x, deg = q)
X <- sapply(poly_orth, predict, newdata = x)
all(abs(t(X) %*% X - diag(q+1)) <= 1e-12)
d_true_beta_orth <- runif(q, min = -3, max = 3)
true_beta_orth <- replace(as.numeric(integral(polynomial(d_true_beta_orth))),1,1)
Y <-  X %*% true_beta_orth + rnorm(N, sd = 0.025)


#save(q,x,poly_orth,X,true_beta_orth,Y,file = "../testing_workspaces/gcreg_early_testing_20170228.RData")
#load("../testing_workspaces/gcreg_early_testing_20170228.RData")

# converters
cv <- gen_poly_basis_converters(poly_orth)

xc <- seq(-1,1,length = 200)

plot(x = xc, y = predict(polynomial(cv$to_monomial(true_beta_orth)), newdata = xc), type = "l")
points(x = x, y = Y)

TYPE <- "increasing"

# identity constaint
oracle_id <- function(...) {T}

opCOLS <- optim_COLS(par = rep(0,times = q+1), Y = Y, X = X, oracle_fun = oracle_id)
opCOLS - coef(lm(Y~-1+X))

# parameter inequality constaint
oracle_pic <- function(par) {par[2,] <= 0}

opCOLS <- optim_COLS(par = rep(-1,times = q+1), Y = Y, X = X, oracle_fun = oracle_pic)
op_BFGS <- optim(par = rep(-1,times = q+1), fn = function(beta) crossprod(Y - X %*% beta), upper = c(Inf,0,rep(Inf, times = q-1)), method = "L-BFGS-B")

opCOLS - op_BFGS$par

# monotinicity  constaint
oracle_mon <- function(par) is_monotone(p = cv$to_monomial(par), region = c(-1,1), type = TYPE, EPS = 1e-05)

oracle_mon2 <- function(par) ismonotone(object = cv$to_monomial(par), a = -1, b = 1, type = TYPE)

mult <- 2 * (TYPE == "increasing" ) - 1

starter_par_basis <- cv$to_basis(c(mean(Y), mult * mean(Y)*100,rep(0,times = q-1)))

starter_par_mono <- coef(lm(Y~x))

starter_par_basis1 <- cv$to_basis(c(starter_par_mono,starter_par_mono[2] * 2:q));oracle_mon(starter_par_basis1)
starter_par_basis2 <- cv$to_basis(c(starter_par_mono,starter_par_mono[2] * 2^(1:(q-1))));oracle_mon(starter_par_basis2)
starter_par_basis3 <- cv$to_basis(c(starter_par_mono,rep(0,times = q -1)));oracle_mon(starter_par_basis3)

starter_new <- cv$to_basis(integral(deriv(polynomial(cv$to_monomial(crossprod(X,Y))+1))+10) + 2)


d_crit <- function(pl,region, EPS = 1e-05){
  a <- min(region)
  b <- max(region)
  
  D_pl <- deriv(polynomial(pl))
  DD_pl <- deriv(D_pl)
  roots <- polyroot(DD_pl)
  roots_re <- Re(roots)
  roots_im <- Im(roots)
  
  # critical roots are those that are real and within the regions (within numerical precision)
  crit_roots <- roots_re[abs(roots_im) < EPS & a + EPS < roots_re & roots_re < b - EPS]
  
  predict(D_pl, newdata = c(crit_roots,region[is.finite(region)]))
  
}

d_crit(cv$to_monomial(starter_new), region = region)

#if(!oracle_mon(starter_par_basis)) starter_par_basis <- cv$to_basis(c(0,1,rep(0,times=q-1)))

ctr_simple <- cols_control(method = "simple", tol = 1e-15, step_start = 0.3, step_increment = 0.05)
ctr_bs <- cols_control(method = "best-step", tol = 1e-15, step_start = 0.3, step_increment = 0.05, maxit = 1000)

opCOLS1  <- optim_COLS(par = starter_par_basis, Y = Y, X = X, oracle_fun = oracle_mon2, control = ctr_simple)
opCOLS2  <- optim_COLS(par = starter_par_basis, Y = Y, X = X, oracle_fun = oracle_mon2, control = ctr_bs)
opCOLS3  <- optim_COLS(par = as.numeric(starter_new), Y = Y, X = X, oracle_fun = oracle_mon, control = ctr_bs)
crossprod(Y - X %*% opCOLS3)


oracle_mon(opCOLS3)

cv$to_monomial(opCOLS3) - cv$to_monomial(opCOLS2)

plot(x = xc, y = predict(polynomial(cv$to_monomial(opCOLS3)), newdata = xc), type = "l")
points(x = x, y = Y)
lines(x = xc, y = predict(polynomial(cv$to_monomial(starter_par_basis)), newdata = xc))

system.time(opMP <- monpol(Y ~ x, degree = q, a = -1, b = 1, monotone = TYPE))

oracle_mon(cv$to_basis(coef(opMP)))
oracle_mon2(cv$to_basis(coef(opMP)))
ismonotone(cv$to_monomial(opCOLS3), a = -1, b = 1)
ismonotone(opMP, a = -1, b = 1)

crossprod(Y - X %*% opCOLS1)
crossprod(Y - X %*% opCOLS2)
crossprod(Y - X %*% opCOLS3)

crossprod(Y - X %*% cv$to_basis(coef(opMP)))
crossprod(Y - X %*% crossprod(X,Y))


is_monotone(cv$to_monomial(opCOLS3))
plot(xc, predict(polynomial(tes), newdata = xc) )
is_monotone(tes)

opCOLS1 - cv$to_basis(coef(opMP))
(opCOLS2 - cv$to_basis(coef(opMP)))/cv$to_basis(coef(opMP))

# need to check monotonicity checking functions for 1D change in opCOLS2. See what movements is_monotone allows versus ismonotone
# read optim code
# check what pot_moves is saying

# cross check is monotone function
ismonotone(cv$to_monomial(opCOLS1), a = -1, b = 1)
ismonotone(cv$to_monomial(opCOLS2), a = -1, b = 1)
is_monotone(coef(opMP), region = c(-1,1), EPS = 1e-15)


opMP$RSS
sum(residuals(opMP)^2)
plot(x = xc, y = predict(polynomial(coef(opMP)), newdata = xc), type = "l")
lines(x = xc, y = predict(polynomial(cv$to_monomial(opCOLS3)), newdata = xc), col = "red")
####

optim_COLS(par = cv$to_basis(coef(opMP)), Y = Y, X = X, oracle_fun = oracle_mon, control = ctr_bs) - cv$to_basis(coef(opMP))
