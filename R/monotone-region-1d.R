library(polynom)
p <- polynomial(runif(8,min = -1,max = 1))
p_i <- polynomial(c(rep(0,7),1))

a <- -Inf
b <- Inf

is_monotone(p)

is_monotone(p_i)


D_p <- deriv(p)
D_p_i <- deriv(p_i)

roots_D_p <- polyroot(D_p)
roots_D_p_i <- polyroot(D_p_i)

reroots_D_p <- Re(roots_D_p)[abs(Im(roots_D_p)) < EPS & a - EPS <  Re(roots_D_p) &  Re(roots_D_p) < b + EPS]

reroots_D_p_i <- unique(Re(roots_D_p_i)[abs(Im(roots_D_p_i)) < EPS & a - EPS <  Re(roots_D_p_i) &  Re(roots_D_p_i) < b + EPS])

crits <- c(a,sort(reroots_D_p_i),b)

max_region_bounds <- numeric()
min_region_bounds <- numeric()

for(i in 2:length(crits)){
  
  if(is.finite(crits[i-1]) & is.finite(crits[i])){
    test_point <- mean(crits[c(i-1,i)])
  } else if(is.finite(crits[i-1])){
    test_point <- crits[i-1] + 100
  } else if(is.finite(crits[i])) {
    test_point <- crits[i] - 100
  } else {
    test_point <- 0
  }

  test_val <- predict(D_p_i,test_point)
  
  if(test_val > 0){
    max_region_bounds <- c(max_region_bounds,crits[c(i-1,i)])
  } else {
    min_region_bounds <- c(min_region_bounds,crits[c(i-1,i)])
  }

}



# si

