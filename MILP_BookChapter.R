# 
rm(list = ls())
gc()
library(ggplot2)
library(ROI)
library(ROI.plugin.glpk)
library(ompr.roi)
library(ompr)
library(magrittr)
library(Rglpk)
library(cluster)
library(reshape2)
suppressPackageStartupMessages(library(dplyr))
# detach("package:simmer", unload = TRUE)
library(Rtsne)
library(Rlab)
n <- 5
t <- 10


### Base Model ----- 
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
  add_variable(x[i, j, k], i = 1:n, j = 1:n, k = 1:t, type = "continuous",lb = 0) %>%
  
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




### Scenario creation ----
archivalDat <-  rweibull(2000,shape=8,scale=3)

fw <- fitdist(archivalDat,"weibull")

newData <-  rweibull(4000,shape=8,scale=3)


totalDat1 <- c(archivalDat,newData[1:1000])

fw1 <- fitdist(totalDat1,"weibull")
fw1$estimate


totalDat2 <- c(archivalDat,newData)


fw2 <- fitdist(totalDat2,"weibull")


fw$estimate
fw1$estimate
fw2$estimate

par(mfrow=c(1,1))
percent=50
rgb.val <- col2rgb("pink")
name="pink"
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)


hist(rweibull(2000,shape=fw$estimate[1],scale=fw$estimate[2]),col = t.col)
hist(rweibull(2000,shape=fw1$estimate[1],scale=fw1$estimate[2]),col="red",add=T,alpha=0.5)
hist(rweibull(2000,shape=fw2$estimate[1],scale=fw2$estimate[2]),col="blue",add=T,alpha=0.5)



mean <- 3
sd <- 1

mean.est <- 2.5
sd.est <- 0.8


shape <-  8
scale <-  3
shape.est <- 8
scale.est <- 3.5


scenarios <- data.frame(sn=0,load1=0,load2=0,load3=0,wspeed1=0,wpspeed2=0,prob=0)


for(i in 1:1000){
  loads <- rnorm(3,mean = mean,sd = sd)
  rweis <- rweibull(2,shape = shape,scale = scale)
  
  load.prob <- 1
  for(j in 1:length(loads))
    load.prob <- dnorm(loads[j],mean=mean.est,sd=sd.est) * load.prob
  
  rgen.prob <- 1
  for(k in 1:length(rweis))
    rgen.prob <- dweibull(rweis[k],shape = shape.est,scale=scale.est)*rgen.prob
  
  scenarios[i,] <- c(i,loads,rweis,rgen.prob*load.prob)
  
}


# ggplot(scenarios %>% filter(prob>0.005)) + geom_point(aes(x=sn,y=prob/max(prob)))

reducted.scenarios <- scenarios %>% filter(prob>0.005) %>% arrange(desc(prob))
reducted.scenarios <- head(reducted.scenarios,200)
# threshold <- mean(as.matrix(dist.scenarios))/2


dist.scenarios <- daisy(reducted.scenarios,"euclidean")
mat.dist.scenarios <- as.matrix(dist.scenarios)


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

noc <- which(sil_width==min(sil_width,na.rm=T))

pam_fit <- pam(mat.dist.scenarios,
               diss = TRUE,
               k = noc)

reducted.scenarios[as.numeric(pam_fit$medoids),]
final.scenarios <- reducted.scenarios[pam_fit$id.med,] %>% mutate(prob = prob/sum(prob))





### TSNE Needs investigation and explanation ----
tsne_obj <- Rtsne(dist.scenarios, dims =2 ,perplexity = 20, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = reducted.scenarios$sn)

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster)) +theme_classic()





### tidyr::spread example -----

library(dplyr)
library(tidyr)
stocks <- data.frame(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)
stocksm <- stocks %>% gather(stock, price, -time)
stocksm %>% spread(stock, price)





### Modification of Model-----



RHO <- matrix(final.scenarios$prob,nrow = 9)
PF <- matrix(sample(runif(100,0.9,1.1),step_size),nrow = step_size)  # Further Randomization of demand
PF2 <- matrix(sample(runif(100,0.95,1.05),step_size),nrow = step_size) # Further Randomization of renewable generation

demandSet <- 1:3
### Creates demand for bus i at time t in scenario S
demand <- function(i,t,s,PF){
  if(i %in% demandSet){
    d <- as.matrix(final.scenarios[s,] %>%  dplyr::select(load1,load2,load3),ncol=3)
    mm <- t(PF %*% d)
    return(as.numeric(mm[i,t]))
  }else
    0
  
}

#Creates renewable generation of the bus i
rgenSet <- 4:5

vin <- 3.5
vr <-  14.5
vout <- 20
pw <- 5

rgen <- function(i, t, s,PF2) {
  if (i %in% rgenSet) {
    d <-
      as.matrix(final.scenarios[s, ] %>%  dplyr::select(wspeed1, wpspeed2),
                ncol = 2)
    mm <- t(PF2 %*% d)
    ws <- as.numeric(mm[i %% 4 + 1, t])
    if (ws < vin  || ws > vout)
      return(0)
    else if (ws >= vin && ws <= vout)
      return((pw * (ws - vin)) / (vr - vin))
    else
      return(pw)
    
  } else
    return(0)
}

for(i in 1:5)
  for(t in 1:10)
    for(s in 1:9)
      print(paste(i,t,s,rgen(i,t,s),sep = "  "))


A <-  matrix(0,ncol=n,nrow=n) # if there is connection between i and j
C <-  matrix(sapply(diag(n),function(x)ifelse(x==1,0,rnorm(n*n,10,2))),ncol = n) # Capacity of line i-j 
G <-  matrix(c(40,0,0,0,20),ncol=5) # Generation capacity of bus i
P <-  matrix(c(1,2,3,4,5),nrow=n,ncol=1) # Price of generating unit of electricity in bus i



diag(A) <- 0
A[lower.tri(A)] <- c(1,0,1,1,0,1,0,1,0,1)
A[upper.tri(A)] = t(A)[upper.tri(A)]
isSymmetric(A)
A

n <- 5
t <- 10
s <- 9

remove("i","j","k")
which(objs %in% ls(envir=.GlobalEnv))


result <- MIPModel() %>% 
  # Amount of electricity transferred from i to j at time k
  add_variable(x[i, j, k], i = 1:n, j = 1:n, k = 1:t,  type = "continuous",lb = 0) %>% 
  
  # Amount of electricity generated at i at time k  
  add_variable(g[i, k],
               i = 1:n, k = 1:t, type = "continuous",lb = 0) %>% 
  # If line i-j selected to be the path to transfer electricity from i to j
  add_variable(r[i, j], 
               i = 1:n, j = 1:n, type = "binary") %>% 
  # Amount of electricity supplied from main grid
  add_variable(b[i, k, l], i = 1:n, k = 1:t, l = 1:s, type = "continuous",lb=0) %>% 
  
  # Line binding constraint
  add_constraint(x[i, j, k] <= r[i,j]* C[i,j] * A[i,j],
                 i = 1:n, j = 1:n, k = 1:t, l = 1:s) %>% 
  
  # Generation capacity constraint 
  add_constraint(g[i, k] <= G[1,i], 
                 i = 1:n , k = 1:t) %>%   
  
  # Balance constraint
  add_constraint(sum_expr(x[i, j, k], j = 1:n, i!=j) <= b[i,k,l] + sum_expr(x[j, i, k], j = 1:n, i!=j) + 
                   g[i, k] + rgen(i, k, l) - demand(i, k, l), i = 1:n, k = 1:t , l = 1:s)

  # Objective Function   
  result <- result %>% 
    set_objective(sum_expr(RHO[l,1] * (P[i,1] * g[i, k] + 100*b[i, k, l]), i = 1:n, k = 1:t, l=1:s),
                  sense = "min") 

  result <- result %>% 
   solve_model(with_ROI("glpk", verbose = TRUE))
  

  result$constraints[[length(result$constraints)-440]]



dd <- c(rnorm(3,3,1),0,0)
rg <- c(0,0,0,rweibull(2,shape = 8, scale = 3))

df <- get_solution(result,x[i,j,k])

xxf <- function(df,i,j,k){
  return( df[which(1==df$i&&1==df$j && 1 == df$k),"value"] )
  
}

params <- list()

df <- as.data.frame(matrix(ncol = 8 ))
colnames(df) <- c("node","step","arrivals","departures","generation","demand","rgen","excessORdeficit")
p <- 1
y <- 1

calculateB <- function(result){
  gg <- get_solution(result,g[i,k])
  xx <- get_solution(result,x[i,j,k])
  
  N <- max(get_solution(result,x[i,j,k])$i)
  P <-  max(get_solution(result,x[i,j,k])$k)
  for (p in 1:P) {
    for (y in 1:N) {
      arrivals <- sum(xx %>% filter(j == y, k == p) %>% select(value))
      departures <- sum(xx %>% filter(i == y, k == p) %>% select(value))
      generation <- as.numeric(gg %>% filter(i==y,k==p) %>% select(value))
      dd <- rnorm(1,3,1)
      rr <- rweibull(1,shape=8,scale=3)
      bb <- arrivals + generation + rr - dd - departures
      asb <- ifelse()
      df <-  rbind(df,c(y,p,arrivals,departures,generation,dd,rr,bb))
      
    }
 
  }
  return(df[-1,])
}



cal <- calculateB(result)

cal %>% group_by(step) %>% summarise(SUM = sum(excessORdeficit)*100)

result$objective_value



sim <- function(dd,rg,result){
  astep <- max(get_solution(result,x[i,j,k])$k) 
  gg <- get_solution(result,g[i,k])
  rr <- get_solution(result,r[i,j]) %>% filter(value>0)
  mm <- matrix(0,ncol=5,nrow=5)
  
  
  for(i in 1:10){
    for(f in 1:5){
      
      
    }
    
  }
  
  
}

### Functions ----
solveOPT <- function(final.scenarios,step_size,PF,PF2){
  line_cap <- 7.5
  A <-  matrix(0,ncol=n,nrow=n) # if there is connection between i and j
  C <-  matrix(sapply(diag(n),function(x)ifelse(x==1,0,line_cap)),ncol = n) # Capacity of line i-j 
  G <-  matrix(c(40,0,0,0,20),ncol=5) # Generation capacity of bus i
  P <-  matrix(c(4,6,3,5,8),nrow=n,ncol=1) # Price of generating unit of electricity in bus i
  RHO <- matrix(final.scenarios$prob,nrow = 9) # Scenario Probabilities
  PF <- matrix(sample(runif(100,0.9,1.1),step_size),nrow = step_size)  # Further Randomization of demand
  PF2 <- matrix(sample(runif(100,0.95,1.05),step_size),nrow = step_size) # Further Randomization of renewable generation
  result <- MIPModel() %>% 
    # Amount of electricity transferred from i to j at time k
    add_variable(x[i, j, k], i = 1:n, j = 1:n, k = 1:t,  type = "continuous",lb = 0) %>% 
    
    # Amount of electricity generated at i at time k  
    add_variable(g[i, k],
                 i = 1:n, k = 1:t, type = "continuous",lb = 0) %>% 
    # If line i-j selected to be the path to transfer electricity from i to j
    add_variable(r[i, j], 
                 i = 1:n, j = 1:n, type = "binary") %>% 
    # Amount of electricity supplied from main grid
    add_variable(b[i, k, l], i = 1:n, k = 1:t, l = 1:s, type = "continuous",lb=0) %>% 
    
    # Line binding constraint
    add_constraint(x[i, j, k] <= r[i,j]* C[i,j] * A[i,j],
                   i = 1:n, j = 1:n, k = 1:t, l = 1:s) %>% 
    
    # Generation capacity constraint 
    add_constraint(g[i, k] <= G[1,i], 
                   i = 1:n , k = 1:t) %>%   
    
    # Balance constraint
    add_constraint(sum_expr(x[i, j, k], j = 1:n, i!=j) <= b[i,k,l] + sum_expr(x[j, i, k], j = 1:n, i!=j) + 
                     g[i, k] + rgen(i, k, l) - demand(i, k, l), i = 1:n, k = 1:t , l = 1:s)
  
    # Objective Function   
    result <- result %>% 
      set_objective(sum_expr(RHO[l,1] * (P[i,1] * g[i, k] + 100*b[i, k, l]), i = 1:n, k = 1:t, l=1:s),
                    sense = "min") 
  
    result <- result %>% 
      solve_model(with_ROI("glpk", verbose = TRUE))
    
    return(result)
}




fresult <- solveOPT(final.scenarios,10,PF,PF2)

mean.est <- c(3,2.5,4)
sd.est <- c(1,1.5,0.5)


demandL = vector(mode="list",length = length(mean.est))
for(i in 1:length(mean.est)){
  demandL[[i]]$mean.est <- mean.est[1]
  demandL[[i]]$sd.est <- sd.est[1]
  
}

shape.est <- c(8,7,6.5)
scale.est <- c(3.5,4,5)

rgenL = vector(mode="list",length = length(shape.est))
for(i in 1:length(mean.est)){
  rgenL[[i]]$shape.est <- shape.est[1]
  rgenL[[i]]$scale.est <- scale.est[1]
  
}


generateScenarios <- function(demandL,rgenL){
 
  scenarios <- data.frame(sn=0,load1=0,load2=0,load3=0,wspeed1=0,wpspeed2=0,prob=0)
  
  for (i in 1:1000) {
    
    loads <- c()
    rweis <- c()
    
    for (j in 1:length(demandL))
      loads <- append(loads, rnorm(1, mean = demandL[[i]]$mean.est, sd = demandL[[i]]$sd.est))
    
    
    for (k in 1:length(rgenL)) 
      rweis <- append(rweis,rweibull(1, shape = rgenL[[k]]$shape.est , scale = rgenL[[k]]$scale.est))
    
    
    load.prob <- 1
    
    for (r in 1:length(loads))
      load.prob <- dnorm(loads[r], mean = demandL[[r]]$mean.est, sd = demandL[[r]]$sd.est) * load.prob
    
    
    rgen.prob <- 1
    
    for (p in 1:length(rweis))
      rgen.prob <- dweibull(rweis[p], shape = rgenL[[p]]$shape.est, scale = rgenL[[p]]$scale.est) * rgen.prob
    
    scenarios[i, ] <- c(i, loads, rweis, rgen.prob * load.prob)
    
  }
  
}



reduceScenarios <- function(scenarios){
  

  reducted.scenarios <- scenarios %>% filter(prob>0.005) %>% arrange(desc(prob))
  reducted.scenarios <- head(reducted.scenarios,200)
  for(i in 2:10){
    
    pam_fit <- pam(mat.dist.scenarios,
                   diss = TRUE,
                   k = i)
    
    sil_width[i] <- pam_fit$silinfo$avg.width
    
  }

  noc <- which(sil_width==min(sil_width,na.rm=T))
  
  pam_fit <- pam(mat.dist.scenarios,
                 diss = TRUE,
                 k = noc)
  
  reducted.scenarios[as.numeric(pam_fit$medoids),]
  final.scenarios <- reducted.scenarios[pam_fit$id.med,] %>% mutate(prob = prob/sum(prob))
  
  
  return(final.scenarios)
  
}





### Network visualization ----
library('igraph')
nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
head(nodes)
head(links)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net


nodes2 <- read.csv("Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
links2 <- read.csv("Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)

dim(links2)
dim(nodes2)
net2 <- graph_from_incidence_matrix(A)

net2
l <- do.call(layout_on_grid, list(net2)) 

layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
colors <- ifelse(get_solution(result,r[i,j])$value >0, "green","grey")
plot(net2, edge.arrow.size=.4,main = layout_on_grid,layout = l)

df <- get_solution(result,r[i,j])
df <- df %>% filter(value>0)
colnames(df) <- c("variable","from","to","value")
df <-  df %>% select(from)

result$objective_value








