library(ROI)
library(ROI.plugin.glpk)
library(ompr.roi)
library(ompr)
library(magrittr)
library(Rglpk)
library(dplyr)

n <- 5
t <- 10
a <- matrix(c(0,1,1,
         1,0,1,
         1,1,0),ncol=3)

set.seed(1)

weights <- matrix(rpois(n * n, 5), ncol = n, nrow = n)
demand  <- c(rnorm(n,0,1))
demand  <- matrix(sapply(demand,function(x) ifelse(x<0,0,x)),nrow=1)
       


result <- MIPModel() %>% 
  add_variable(x[i, j, k], i = 1:n, j = 1:n, k=1:t, type = "continuous") %>% ## Amount of electricity transferred from i to j at time k
  add_vaiable(g[i,k],i=1:n,k=1:t,type = "continuous") %>% ## Amount of electricity generated at i at time k  
  
  set_objective(sum_expr(weights[i, j] * x[i, j, k], i = 1:n, j = 1:n,k=1:n)) %>% 
  add_constraint( sum_expr(x[i, j, k], k = 1:n) <= 1, i = 1:n,j=1:n) %>% 
  add_constraint(x[i,j,k]<=5,i=1:n,j=1:n,k=1:n) %>% 
  solve_model(with_ROI("glpk", verbose = TRUE))

get_solution(result, x[i, j, k]) %>% filter(value>0,variable=="x",i==1)


