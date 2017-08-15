# testing for functions in R/is-monotone.R script

test_that("is_monotone correct over entire real line", {
  
  p1 <- polynomial(1:6)
  
  expect_true(is_monotone(p1))

})
