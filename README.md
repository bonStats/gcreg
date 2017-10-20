# gcreg
General Constraint Regression models in `R`

This package is currently being developed. It's aim is to provide methods for fitting regression models with:
* Functional and shape constraints, e.g. monotonicity
* Parameter inequality constraints
* Joint constraints, e.g. combinations of the above
* Other constraints that create closed and convex parameter spaces

The current focus of development is on monotonicity in polynomial fixed and mixed effects models but will be extended over time to more general models and constraints. 

To get started, install this package from GitHub using the `devtools` package:

```r
devtools::install_github("bonStats/gcreg")
library(gcreg)
```

You can start fitting constrained polynomial models with the `gcreg::cpm()` function. For example

```r 
library(fda)
data(onechild)
cpm(height~day, data = onechild, degree = 5, constraint = "monotone", c_region = c(1,312))
```
See the package vignettes for more examples.