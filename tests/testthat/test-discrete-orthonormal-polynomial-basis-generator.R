context("testing functions in R/discrete-orthonormal-polynomial-basis-generator.R script")

test_that("make_disc_orthonormal_basis() makes a discrete orthonormal basis", {
  
  for(num in c(10,100,1000,10000)){
  
    x <- runif(n = num, min = -1, max = 1)
    
    for(d in sample(x = 2:min(num-2,20), size = 3, replace = F)){
    
      p_basis <- make_disc_orthonormal_basis(x = x, deg = d)
      
      Xo <- sapply(p_basis, predict, newdata = x)
      
      expect_equal(crossprod(Xo),diag(nrow = d+1), tolerance = 1e-08)
    
      }    
  
  }
  
  })
