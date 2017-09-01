context("testing functions in R/constrained-orthonormal-least-squares.R script")

test_that("optim_cols() converges to same solution as lm() on unconstrained problem",{
  
  n_test <-c(20,400,8000)
  p_test <- seq(from = 4, to = 13, by = 3)
  
  for(n in n_test){
    
    for(p in p_test){
      
      x <- runif(n, min = -1, max = 1)
      
      b <- rnorm(p)
      
      p_basis <- make_disc_orthonormal_basis(x = x, deg = p - 1)
      
      Xo <- sapply(p_basis, predict, newdata = x)
      
      t_mean <-  Xo%*%b
      
      y <- t_mean + rnorm(n,sd = diff(range(t_mean))/4)
      
      cv_y <- gen_scale_data_funs(y)
      
      Y <- cv_y$scale(y)
      
      sim_dat <- data.frame(Y,Xo)
      
      lm_coef <- as.numeric(coef(lm(Y~.-1, data = sim_dat)))
      
      
      cols_coef <- optim_cols(par = rep(0,times = p), Y = Y, X = Xo, 
                              oracle_fun = function(x) {T}
      )[,1]
      
      expect_equal(cols_coef, lm_coef, tolerance = 1e-05, info = paste("Different coefficients from lm and optim_cols, with p =", p, "and n =", n))
      
      
    }
    
  }
  
  
})


test_that("gcreg::optim_cols() vs MonoPoly::monpol() on general constrained monotonic  problem",{
  
  library(MonoPoly)
  
  n_test <-c(200,400,8000)
  poly_test <- list(
    c(0.5, 0.25, 0, - 0.5, 0, 0.4),
    c(-2, 1, 0.475, -3.167, -0.872, 6.255, 0.506, -3.182)
  )
  
  for(n in n_test){
    
    for(b in poly_test){
      
      compact_lims <- c(runif(1, min = -1, max = -0.2), runif(1, min = 0.2, max = 1))
      
      clim <- sample(list(c(-Inf,Inf),compact_lims),size = 1)[[1]]
      
      x <- runif(n, min = -1, max = 1)
      
      p_basis <- make_disc_orthonormal_basis(x = x, deg = length(b) - 1)
      
      cv <- gen_poly_basis_converters(poly_basis = p_basis)
      
      b_o <- cv$to_ortho(b)
      
      Xo <- sapply(p_basis, predict, newdata = x)
      
      t_mean <-  Xo%*%b_o 
      
      y <- t_mean + rnorm(n,sd = diff(range(t_mean))/10)
      
      cv_y <- gen_scale_data_funs(y)
      
      Y <- cv_y$scale(y)
      
      monpol_coef <- cv$to_ortho(
        coef(
          monpol(formula = Y ~ x, degree = length(b) - 1, a = clim[1], b = clim[2])
        )
      )
      

      cols_coef <- optim_cols(par = cv$to_ortho(rep(c(0.1,1),times = length(b)/2)), Y = Y, X = Xo,
                              function(p) is_monotone(p = cv$to_mono(p), region = clim),
                              control = cols_control(tol = 1e-6, method = "best-step",step_start = 0.7,step_increment = 0.05)
      )[,1]
      
      cols_RSS_better <- sum((Y - Xo%*%cols_coef)^2) < sum((Y - Xo%*%monpol_coef)^2)
      
      RSS_close <- abs(sum((Y - Xo%*%cols_coef)^2) / sum((Y - Xo%*%monpol_coef)^2) - 1) < 1e-02
      
      cols_monpol_equal <- all.equal(cols_coef, monpol_coef, tolerance = 1e-02)
      
      expect(cols_RSS_better | RSS_close | is.logical(cols_monpol_equal), message = paste("monpol out performed cols by more than 1% and different coefficients were found, with p =", length(b), "and n =", n))
    
      expect_true(is_monotone(cv$to_mono(cols_coef), region = clim))  
    
    }
    
  }
  
})

test_that("optim_cols() converges to same solution as monpol() on gap constrained monotonic  problem",{
  
  library(MonoPoly)
  
  n_test <-c(20,400,8000)
  poly_test <- list(
    c(0.5, 0.25, 0, - 0.5, 0, 0.4),
    c(-2, 1, 0.475, -3.167, -0.872, 6.255, 0.506, -3.182)
  )
  
  for(n in n_test){
    
    for(b in poly_test){
      
      compact_lims <- c(runif(1, min = -1, max = -0.2), runif(1, min = 0.2, max = 1))
      
      clim <- sample(list(c(-Inf,Inf),compact_lims),size = 1)[[1]]
      
      x <- c(runif(floor(n/2), min = -1, max = 0.3),runif(ceiling(n/2), min = 0.7, max = 1))
      
      p_basis <- make_disc_orthonormal_basis(x = x, deg = length(b) - 1)
      
      cv <- gen_poly_basis_converters(poly_basis = p_basis)
      
      b_o <- cv$to_ortho(b)
      
      Xo <- sapply(p_basis, predict, newdata = x)
      
      t_mean <-  Xo%*%b_o 
      
      y <- t_mean + rnorm(n,sd = diff(range(t_mean))/10)
      
      cv_y <- gen_scale_data_funs(y)
      
      Y <- cv_y$scale(y)
      
      monpol_coef <- cv$to_ortho(
        coef(
          monpol(formula = Y ~ x, degree = length(b) - 1, a = clim[1], b = clim[2])
        )
      )
      
      cols_coef <- optim_cols(par = cv$to_ortho(rep(c(0.1,1),times = length(b)/2)), Y = Y, X = Xo, 
                              function(p) is_monotone(p, region = clim),
                              control = cols_control(tol = 1e-6, method = "best-step")
      )[,1] 
      
      cols_RSS_better <- sum((Y - Xo%*%cols_coef)^2) < sum((Y - Xo%*%monpol_coef)^2)
      
      RSS_close <- abs(sum((Y - Xo%*%cols_coef)^2) / sum((Y - Xo%*%monpol_coef)^2) - 1) < 1e-02
      
      cols_monpol_equal <- all.equal(cols_coef, monpol_coef, tolerance = 1e-02)
      
      expect(cols_RSS_better | RSS_close | is.logical(cols_monpol_equal), message = paste("monpol out performed cols by more than 1% and different coefficients were found, with p =", length(b), "and n =", n))
      
    }
    
  }
  
})


test_that("optim_cols() converges to same solution as monpol() on gap + outlier constrained monotonic  problem",{
  
  library(MonoPoly)
  
  n_test <-c(20,400,8000)
  poly_test <- list(
    c(0.5, 0.25, 0, - 0.5, 0, 0.4),
    c(-2, 1, 0.475, -3.167, -0.872, 6.255, 0.506, -3.182)
  )
  
  for(n in n_test){
    
    for(b in poly_test){
      
      compact_lims <- c(runif(1, min = -1, max = -0.2), runif(1, min = 0.2, max = 1))
      
      clim <- sample(list(c(-Inf,Inf),compact_lims),size = 1)[[1]]
      
      x <- c(runif(floor(n/2), min = -1, max = 0.3),runif(ceiling(n/2), min = 0.7, max = 1))
      
      p_basis <- make_disc_orthonormal_basis(x = x, deg = length(b) - 1)
      
      cv <- gen_poly_basis_converters(poly_basis = p_basis)
      
      b_o <- cv$to_ortho(b)
      
      Xo <- sapply(p_basis, predict, newdata = x)
      
      t_mean <-  Xo%*%b_o 
      
      y <- t_mean + rnorm(n,sd = diff(range(t_mean))/10)
      
      cv_y <- gen_scale_data_funs(y)
      
      Y <- cv_y$scale(y)
      
      outliers_x <- c(0.5,0.7)
      outliers_Xo <- sapply(p_basis, predict, newdata = outliers_x) 
      outliers_y <- outliers_Xo %*% b_o - runif(2, min = 0, max = 0.5)
      
      x <- c(x,outliers_x)
      p_basis <- make_disc_orthonormal_basis(x = x, deg = length(b) - 1)
      Xo <- sapply(p_basis, predict, newdata = x)
      Y <- rbind(Y,outliers_y)
      
      monpol_coef <- cv$to_ortho(
        coef(
          monpol(formula = Y ~ x, degree = length(b) - 1, a = clim[1], b = clim[2])
        )
      )
      
      cols_coef <- optim_cols(par = cv$to_ortho(rep(c(0.1,1),times = length(b)/2)), Y = Y, X = Xo, 
                              function(p) is_monotone(p, region = clim),
                              control = cols_control(tol = 1e-6, method = "best-step")
      )[,1] 
      
      cols_RSS_better <- sum((Y - Xo%*%cols_coef)^2) < sum((Y - Xo%*%monpol_coef)^2)
      
      RSS_close <- abs(sum((Y - Xo%*%cols_coef)^2) / sum((Y - Xo%*%monpol_coef)^2) - 1) < 1e-02
      
      cols_monpol_equal <- all.equal(cols_coef, monpol_coef, tolerance = 1e-02)
      
      expect(cols_RSS_better | RSS_close | is.logical(cols_monpol_equal), message = paste("monpol out performed cols by more than 1% and different coefficients were found, with p =", length(b), "and n =", n))
      
    }
    
  }
  
})
