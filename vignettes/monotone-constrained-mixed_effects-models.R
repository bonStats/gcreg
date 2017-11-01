## ---- fig.show='hold', fig.cap='fda::growth data', cache=TRUE------------
library(fda)
library(ggplot2)
library(dplyr)

# format data
boys <- fda::growth$hgtm %>% melt(value.name = "height", varnames = c("age","id"))

# take a random sample of 10 boys
set.seed(2017)
ids <- boys$id %>% as.character() %>% unique() %>% sample(size = 10, replace = F)
boys_s <- boys %>% filter(id %in% ids)

#scale data
b_age_sc <- gen_scale_data_funs(boys_s$age)
b_height_sc <- gen_scale_data_funs(boys_s$height)
  
boys_s <- boys_s %>% mutate(age2 = b_age_sc$scale(age), height2 = b_height_sc$scale(height))

# set up model we wish to fit
md <- gcreg:::make_em_model_specs(height2~age2, data = boys_s,
                                    p_degree = 12, r_degree = 5,
                                    r_constrained = F, # should random effects be constrained
                                    mcontr_region = c(-1,1),
                                    group_name = "id"
                                    )

# verbose = T shows the steps of the quasi-log-likelihood
boys_em <- gcreg:::constrained_lmm_em(model = md, tol = 1e-02, verbose = T, save_steps = F)



