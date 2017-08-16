context("testing functions in R/is-monotone.R script")

test_that("is_monotone() true positive results over entire real line", {
  
  p <- polylist()
  p[[1]] <- polynomial(1:6)
  p[[2]] <- polynomial(-(1:6))
  p[[3]] <- polynomial(c(-1,1))
  p[[4]] <- polynomial(c(-1,-1))
  p[[5]] <- polynomial(c(-1,0,0,0,0,0.001))
  p[[6]] <- polynomial(c(-1,0,0,0,0,-0.001))
  
  for(i in 1:length(p)){
    
    expect_true(is_monotone(p[[i]]), info = paste("The polynomial:", p[[i]],"should be monotonic over the real line."))
    
  }

  for(i in 1:10){
    
    # only for odd degree polynomials, increasing
    deg <- 2 * i + 1
    coefs <- runif(deg,min = -1, max = 1)
    coefs <- c(coefs, min(abs(max(coefs) + 0.1), 1))
    
    pl <- polynomial(coefs)
    D_pl <- deriv(pl)
    DD_pl <- deriv(D_pl)
    roots_DD_pl <- polyroot(DD_pl)
    real_roots_DD_pl <- Re(roots_DD_pl)[abs(Im(roots_DD_pl)) < 1e-10]
    
    if(all(abs(Im(polyroot(D_pl))) > 1e-10)){ # if derivative has no real roots
      
      expect_true(is_monotone(pl), info = paste("The polynomial:", pl,"should be monotonic over the real line."))
      
    } else {
      
      # adjust random polynomial so that it is monotone
      min_D_pl <- min(predict(D_pl, newdata = real_roots_DD_pl))
      new_pl <- pl + polynomial(c(0,-min_D_pl*1.01))
      
      expect_true(is_monotone(new_pl), info = paste("The polynomial:", new_pl,"should be monotonic over the real line."))
      
      
    }
    
  }  
  
})


test_that("is_monotone() true negative results over entire real line", {
  
  p <- polylist()
  p[[1]] <- polynomial(c(1:5,-6))
  p[[2]] <- polynomial(-c(1:5,-6))
  p[[3]] <- polynomial(c(-1,1,-1))
  p[[4]] <- polynomial(c(-1,-1,1))
  p[[5]] <- polynomial(c(-1,0,0,0,1,0.001))
  p[[6]] <- polynomial(c(-1,0,0,0,1,-0.001))
  p[[7]] <- polynomial(runif(7, min = -2, max = 2))
  p[[8]] <- polynomial(runif(9, min = -2, max = 2))
  p[[9]] <- polynomial(runif(11, min = -2, max = 2))
  
  for(i in 1:length(p)){
    
    expect_false(is_monotone(p[[i]]), info = paste("The polynomial:", p[[i]],"should not be monotonic over the real line."))
    
  }
  
})


test_that("is_monotone() true positive results over compact or semi-compact space", {
  
    # using isotonic parametrisation from Murray, K. (2015). [Thesis] "Improved monotone polynomial fitting with applications and variable selection"
    
    prec <- 0.01 # precision
    
    for(i in 1:10){
      r <- list()
      r[[1]] <- sort(sample(-10:10,size = 2,replace = F))
      r[[2]] <- c(r[[1]][1],Inf)
      r[[3]] <- c(-Inf,r[[1]][2])
      
      pr1 <- polynomial(c(-r[[1]][1] + prec, 1)) * polynomial(c(r[[1]][2] + prec, -1))
      pr2 <- polynomial(c(-r[[1]][1] + prec, 1))
      pr3 <- polynomial(c(r[[1]][2] + prec, -1))
      p1 <- polynomial(runif(i, min = -2, max = 2))
      p2 <- polynomial(runif(i, min = -2, max = 2))
      scl <- runif(1, min = -2, max = 2)
      interc <- runif(1, min = -2, max = 2)
      
      pl <- polylist()
      pl[[1]] <- scl * integral(p1 * p1 + pr1 * p2 * p2) + interc # monotonic increasing over r1
      pl[[2]] <- scl * integral(p1 * p1 + pr2 * p2 * p2) + interc # monotonic increasing over r2
      pl[[3]] <- - scl * integral(p1 * p1 + pr3 * p2 * p2) + interc # monotonic decreasing over r3
      
      for(j in 1:3){
        
        expect_true(is_monotone(pl[[j]], region = r[[j]]), info = paste("The polynomial:", pl[[j]],"should be monotonic over", paste(sort(r[[j]]),collapse = ", ")))
        
      }
  
    }
  
  })


test_that("is_monotone() true negative results over compact or semi-compact space", {
  
  p <- polylist()
  p[[1]] <- integral(poly.from.roots(c(-1,-1,1,1,5)))
  p[[2]] <- integral(poly.from.roots(c(-1,-1,1,1,2,2,5)))
  p[[3]] <- integral(poly.from.roots(c(-1,-1,1,1,2,2,3,3,5)))

  
  for(i in 1:length(p)){
    
    expect_false(is_monotone(p[[i]], region = c(-Inf, 5.01)), info = paste("The polynomial:", p[[i]],"should not be monotonic over -Inf, 5.01."))
    expect_false(is_monotone(p[[i]], region = c(0, 5.01)), info = paste("The polynomial:", p[[i]],"should not be monotonic over 0, 5.01."))
    expect_false(is_monotone(p[[i]], region = c(0, Inf)), info = paste("The polynomial:", p[[i]],"should not be monotonic over 0, Inf"))
    expect_false(is_monotone(p[[i]], region = c(-Inf, Inf)), info = paste("The polynomial:", p[[i]],"should not be monotonic over the real line."))
    
  }
  
})