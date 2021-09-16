result <- MIPModel() %>% 
  # Amount of electricity transferred from i to j at time k
  add_variable(x[i, j, k, l], i = 1:n, j = 1:n, k = 1:t, l = 1:s,  type = "continuous",lb = 0) %>% 
  
  # Amount of electricity generated at i at time k  
  add_variable(g[i, k],
               i = 1:n, k = 1:t, type = "continuous",lb = 0) %>% 
  # If line i-j selected to be the path to transfer electricity from i to j
  add_variable(r[i, j], 
               i = 1:n, j = 1:n, type = "binary") %>% 
  # Amount of electricity supplied from main grid
  add_variable(b[i, k, l], i = 1:n, k = 1:t, l = 1:s, type = "continuous",lb=0) %>% 
  
  # Line binding constraint
  add_constraint(x[i, j, k, l] <= r[i,j]* C[i,j] * A[i,j],
                 i = 1:n, j = 1:n, k = 1:t, l = 1:s) %>% 
  
  # Generation capacity constraint 
  add_constraint(g[i, k] <= G[i], 
                 i = 1:n , k = 1:t) %>%   
  
  # Balance constraint
  add_constraint(sum_expr(x[i, j, k, l], j = 1:n, i!=j) <= b[i,k,l] + sum_expr(x[j, i, k, l], j = 1:n, i!=j) + 
                   g[i, k] + rgen(i, k, l) - demand(i, k, l), i = 1:n, k = 1:t , l = 1:s) %>% 
  
  add_constraint(sum_expr(x[j, i, k, l], j = 1:n, i!=j) + 
                   g[i, k] + rgen(i, k, l) >= demand(i,k,l), i = 1:n, k = 1:t , l = 1:s)

# Objective Function   
result <- result %>% 
  set_objective(sum_expr(RHO[l,1] * (P[i,1] * g[i, k] + 100*b[i, k, l]), i = 1:n, k = 1:t, l=1:s),
                sense = "min") 


result$variables
remove(result)
result$constraints[[1]]
result$constraints[[length(result$constraints)-440]]
result$objective

result <- result %>% 
  solve_model(with_ROI("glpk", verbose = TRUE))


get_solution(result,x[i,j,k,l]) %>% filter(i==1,j==2)
get_solution(result,r[i,j]) %>% filter(value>0)
get_solution(result,b[i,k,l])
get_solution(result,g[i,k])
get_solution(result,g[i,k])

result$objective_value


get_solution(result,r[i,j]) %>% filter(value>0)