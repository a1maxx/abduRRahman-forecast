# 
# rm(list = ls())
# gc()

library(ROI)
library(ROI.plugin.glpk)
library(ompr.roi)
library(ompr)
library(magrittr)
library(Rglpk)

suppressPackageStartupMessages(library(dplyr))
n <- 5
t <- 10


###
N <- 1:n # Set of buses
K <- 1:t # Set of time segments



A <-  matrix(sapply(diag(n),function(x) ifelse(x==1,0,1)),ncol=n) # if there is connection between i and j
G <-  matrix(c(40,0,0,0,20),ncol=5) # Generation capacity of bus i
D <-  matrix(as.numeric(rep(c(10,10,10,10,10),t)),ncol=t,nrow=n) # Demand of the bus i
R <-  matrix(c(0,0,0,5,2),nrow=n,ncol=1) # Amount of renewable generation of the bus i
P <-  matrix(c(1,2,3,4,5),nrow=n,ncol=1) # Price of generating unit of electricity in bus i
C <-  matrix(sapply(diag(n),function(x)ifelse(x==1,0,rnorm(n*n,10,2))),ncol = n) # Capacity of line i-j 





result <- MIPModel() %>% 
  ## Amount of electricity transferred from i to j at time k
  add_variable(x[i, j, k], i = N, j = N, k = K, type = "continuous",lb = 0) %>%
  
  ## Amount of electricity generated at i at time k  
  add_variable(g[i, k],
               i = 1:n, k = 1:t, type = "continuous",lb = 0) %>% 
  ## If line i-j selected to be the path to transfer electricity from i to j
  add_variable(r[i, j], 
               i = 1:n, j = 1:n, type = "binary") %>% 
 
  add_variable(b[i,k], i = 1:n, k = 1:t, type = "continuous",lb=0) %>% 
  
  # Line binding constraint
  add_constraint(x[i,j,k] <= r[i,j]* C[i,j],
                 i= 1:n, j= 1:n, k = 1:t) %>% 
  
  ## Generation capacity constraint 
  add_constraint(g[i,k] <= G[i], 
                 i = 1:n , k = 1:t) %>% 
  
  # Balance constraint
  add_constraint(sum_expr(x[i, j, k], j = 1:n, j!=i) <= b[i,k] + sum_expr(x[j, i, k], j = 1:n, j!=i ) + g[i,k] + R[i,1] - D[i,k]  , 
                i = 1:n, k = 1:t) %>% 
  
  ## Path constraint
  # add_constraint(r[i,j] <= 1, i=1:n,j=1:n,i!=j)

  ## Objective Function   
  set_objective(sum_expr(P[i,1] * g[i, k], i = 1:n, k = 1:t) + sum_expr(5 * b[i, k], i = 1:n, k = 1:t),sense = "min") 
  
result <- result %>% 
  solve_model(with_ROI("glpk", verbose = TRUE))

get_solution(result,x[i,j,k]) %>% filter(i==1,j==2)
get_solution(result,r[i,j])
get_solution(result,b[i,k])


#### Scenario creation

fw <- fitdist(rweibull(2000,shape=8,scale=3),"weibull")


mean <- 3
sd <- 1
mean.est <- 2.5
sd.est <- 0.8


shape <-  8
scale <-  3
shape.est <- 7.5
scale.est <- 2.5

vci <- 1
vr <- 5
pr <-  5

scenarios <- data.frame(sn=0,load1=0,load2=0,load3=0,wspeed1=0,wpspeed2=0,prob=0)


for(i in 1:1000){
  loads <-rnorm(3,mean = mean,sd = sd)
  rweis <- rweibull(2,shape = shape,scale = scale)
  
  load.prob <- 1
  for(j in 1:length(loads))
    load.prob <- dnorm(loads[j],mean=mean.est,sd=sd.est) * load.prob
  
  rgen.prob <- 1
  for(k in 1:length(rweis))
    rgen.prob <- dweibull(rweis[k],shape = shape.est,scale=scale.est)*rgen.prob

  scenarios[i,] <- c(i,loads,rweis,rgen.prob*load.prob)
    
}

head(scenarios,5)

# ggplot(scenarios %>% filter(prob>0.005)) + geom_point(aes(x=sn,y=prob/max(prob)))

nrow(scenarios)
reducted.scenarios <- scenarios %>% filter(prob>0.005) %>% arrange(desc(prob))
reducted.scenarios <- head(reducted.scenarios,60)
nrow(reducted.scenarios)

dist.scenarios <- daisy(reducted.scenarios,"euclidean")
mat.dist.scenarios <- as.matrix(dist.scenarios)
df.dist.scenarios <- as.data.frame(mat.dist.scenarios)
melted.df <- melt(df.dist.scenarios)
melted.df <- melted.df %>% arrange(desc(value))
tail(melted.df)
threshold <- mean(as.matrix(dist.scenarios))/2

kmeans( dist.scenarios)
sil_width <- c(NA)
for(i in 2:10){
  
  pam_fit <- pam(mat.dist.scenarios,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

# plot(1:10, sil_width,
#      xlab = "Number of clusters",
#      ylab = "Silhouette Width")
# lines(1:10, sil_width)


pam_fit <- pam(mat.dist.scenarios,
               diss = TRUE,
               k = 9)

reducted.scenarios[as.numeric(pam_fit$medoids),]

##### Needs investigation and explanation
tsne_obj <- Rtsne(dist.scenarios, dims =2 ,perplexity = 10, is_distance = TRUE)


tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = reducted.scenarios$sn)

# ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster)) +theme_classic()

final.scenarios <- reducted.scenarios[pam_fit$id.med,] %>% mutate(prob = prob/sum(prob))



### Scenario based parameters-----
n <- 5
t <- 10
s <- 9

RHO <- matrix(final.scenarios$prob,nrow = 9)
PF <- matrix(seq(0.1,1,0.1),nrow = 10) 


demandSet <- 1:3
### Creates demand for bus i at time t in scenario S
demand <- function(i,t,s){
  if(i %in% demandSet){
  d <- as.matrix(final.scenarios[s,] %>%  select(load1,load2,load3),ncol=3)
  mm <- t(PF %*% d)
  return(as.numeric(mm[i,t]))
  }else
    0
    
}

#Creates renewable generation of the bus i
rgenSet <- 4:5
rgen <- function(i,t,s){
  if(i %in% rgenSet){
    d <- as.matrix(final.scenarios[s,] %>%  select(wspeed1,wpspeed2),ncol=2)
    mm <- t(PF %*% d)
    as.numeric(mm[i%%4+1,t])
  }else
    0
}



A <-  matrix(sapply(diag(n),function(x) ifelse(x==1,0,1)),ncol=n) # if there is connection between i and j
C <-  matrix(sapply(diag(n),function(x)ifelse(x==1,0,rnorm(n*n,10,2))),ncol = n) # Capacity of line i-j 

G <-  matrix(c(40,0,0,0,20),ncol=5) # Generation capacity of bus i
P <-  matrix(c(1,2,3,4,5),nrow=n,ncol=1) # Price of generating unit of electricity in bus i



result <- MIPModel() %>% 
  ## Amount of electricity transferred from i to j at time k
  add_variable(x[i, j, k, l], i = 1:n, j = 1:n, k = 1:t, l = 1:s,  type = "continuous",lb = 0) %>% 
  
  ## Amount of electricity generated at i at time k  
  add_variable(g[i, k, l],
               i = 1:n, k = 1:t, l = 1:s, type = "continuous",lb = 0) %>% 
  ## If line i-j selected to be the path to transfer electricity from i to j
  add_variable(r[i, j], 
               i = 1:n, j = 1:n, type = "binary") %>% 
  
  add_variable(b[i, k, l], i = 1:n, k = 1:t, l = 1:s, type = "continuous",lb=0) %>% 
  
  # Line binding constraint
  add_constraint(x[i, j, k, l] <= r[i,j]* C[i,j],
                 i = 1:n, j = 1:n, k = 1:t, l = 1:s) %>% 
  
  ## Generation capacity constraint 
  add_constraint(g[i, k, l] <= G[i], 
                 i = 1:n , k = 1:t, l = 1:s) %>% 
  
  # Balance constraint
  add_constraint(sum_expr(x[i, j, k, l], j = 1:n, j!=i) <= b[i,k,l] + sum_expr(x[j, i, k, l], j = 1:n, j!=i ) + 
                   g[i, k, l], 
                 # + rgen(i, k, l) - demand(i, k, l), 
                 i = 1:n, k = 1:t , l = 1:s) %>% 
  
  ## Path constraint
  # add_constraint(r[i,j] <= 1, i=1:n,j=1:n,i!=j)
  
  ## Objective Function   
  set_objective(sum_expr(RHO[l,1] * P[i,1] * g[i, k, l] + b[i, k, l], i = 1:n, k = 1:t, l=1:s),
                sense = "min") 


result$objective

result <- result %>% 
  solve_model(with_ROI("glpk", verbose = TRUE))

get_solution(result,x[i,j,k]) %>% filter(i==1,j==2)
get_solution(result,r[i,j])
get_solution(result,b[i,k])

