rm(list = ls())
gc()

library(fitdistrplus)
load1 <- rnorm(250,3,1)
load2 <- rnorm(250,2.5,1.5)
load3 <- rnorm(250,4,0.5)
load1 <- ifelse(load1<0,0,load1)
load2 <- ifelse(load2<0,0,load2)
load3 <- ifelse(load3<0,0,load3)

wspeed1 <- rweibull(250,shape= 8,scale=3.5) 
wspeed2 <- rweibull(250,shape= 7,scale=4)
df <- data.frame(load1,load2,load3,wspeed1,wspeed2)




rgenL <-  list()
demandL <- list()


m <- generateEstimates(df)



assignEstimates(m)

upd.m <-  updateEstimates(archival,new)

assignEstimates(upd.m)
scenarios <- generateScenarios(demandL,rgenL)
reduced.scenes <- reduceScenarios(scenarios)
result <- solveOPT(final.scenarios = reduced.scenes,step_size = 10)



## Scenario Reduction Function  ----

assign(paste("load",1,sep = ""),rnorm(1,2,1))

d <- 3
r <- 2

createData <-  function(l,d,r){
  
  dl <- l[1:d]  
  rl <- l[(d+1):(d+r)]
  dl.names <- paste("load",1:length(dl),sep="")
  for(i in 1:length(dl.names)){
    assign(dl.names[i],rnorm(250,dl[[i]]$mean.est,dl[[i]]$sd.est))
    
  }
  
  rl.names <- paste("wspeed",1:length(rl),sep="")
  for(i in 1:length(rl.names)){
    assign(rl.names[i],rnorm(250,rl[[i]]$shape.est,rl[[i]]$scale.est))
    
  }
  
  ls <- ifelse(sapply(dl.names,function(x) get(x,envir = parent.frame()))<0,0,sapply(dl.names,function(x) get(x,envir = parent.frame())))
  rs <- ifelse(sapply(rl.names,function(x) get(x,envir = parent.frame()))<0,0,sapply(rl.names,function(x) get(x,envir = parent.frame())))
  
  df <- data.frame(cbind(ls,rs))
  
  
}



generateEstimates <- function(df){
  loads <- colnames(df)[which(grepl("load",colnames(df)))]
  wss <- colnames(df)[which(grepl("wspe",colnames(df)))]
  
  demandL <- vector(mode="list",length=length(loads))
  for(i in 1:length(loads)){
    fit <- fitdist(df[,loads[i]],"norm")
    demandL[[i]]$mean.est <- as.numeric(fit$estimate[1])
    demandL[[i]]$sd.est <- as.numeric(fit$estimate[2])
    
  }
  rgenL <- vector(mode="list",length=length(wss))
  for(i in 1:length(wss)){
    fit <- fitdist(df[,wss[i]],"weibull")
    rgenL[[i]]$shape.est <- as.numeric(fit$estimate[1])
    rgenL[[i]]$scale.est <- as.numeric(fit$estimate[2])
    
  }
  
  return(c(demandL,rgenL))
  
}



assignEstimates <-function(m){
  mean.est <- sapply(m[1:3],function(x) x$mean.est)
  sd.est <- sapply(m[1:3],function(x) x$sd.est)
  demandL <<-  m[1:length(mean.est)]

  shape.est <- sapply(m[(length(mean.est)+1):length(m)],function(x) x$shape.est)
  scale.est <- sapply(m[(length(mean.est)+1):length(m)],function(x) x$scale.est)
  rgenL <<- m[(length(mean.est)+1):length(m)]
}

updateEstimates <- function(archival,new){
  df <-  rbind(archival,new)
  upd.m <- generateEstimates(df)
  return(upd.m)
}






mean.est <- sapply(m[1:3],function(x) x$mean.est)
sd.est <- sapply(m[1:3],function(x) x$sd.est)
demandL <-  m[1:length(mean.est)]

shape.est <- sapply(m[(length(mean.est)+1):length(m)],function(x) x$shape.est)
scale.est <- sapply(m[(length(mean.est)+1):length(m)],function(x) x$scale.est)
rgenL <- m[(length(mean.est)+1):length(m)]


generateScenarios <- function(demandL,rgenL){
  
  scenarios <- data.frame(sn=0,load1=0,load2=0,load3=0,wspeed1=0,wpspeed2=0,prob=0)
  
  for (i in 1:1000) {
    
    loads <- c()
    rweis <- c()
    
    for (j in 1:length(demandL))
      loads <- append(loads, rnorm(1, mean = demandL[[j]]$mean.est, sd = demandL[[j]]$sd.est))
    
    
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
  return(scenarios)
}


reduceScenarios <- function(scenarios){
  reducted.scenarios <- scenarios %>% filter(prob>0.005) %>% arrange(desc(prob))
  reducted.scenarios <- head(reducted.scenarios,200)
  dist.scenarios <- daisy(reducted.scenarios,"euclidean")
  mat.dist.scenarios <- as.matrix(dist.scenarios)
  df.dist.scenarios <- as.data.frame(mat.dist.scenarios)
 
  sil_width <- c() 
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
step_size <- 10
final.scenarios <- reduced.scenes
#### Optimization Function
solveOPT <- function(final.scenarios,step_size){
  n <- 5
  t <- 10
  s <- nrow(final.scenarios)
  line_cap <- 7.5
  A <-  matrix(0,ncol=n,nrow=n) # if there is connection between i and j
  A[lower.tri(A)] <- c(1,0,1,1,0,1,0,1,0,1)
  A[upper.tri(A)] = t(A)[upper.tri(A)]
  C <-  matrix(sapply(diag(n),function(x)ifelse(x==1,0,line_cap)),ncol = n) # Capacity of line i-j 
  G <-  matrix(c(40,0,0,0,20),ncol=5) # Generation capacity of bus i
  P <-  matrix(c(4,6,3,5,8),nrow=n,ncol=1) # Price of generating unit of electricity in bus i
  RHO <- matrix(final.scenarios$prob,nrow = nrow(final.scenarios)) # Scenario Probabilities
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
                     g[i, k] + rgen(i, k, l,PF2) - dgen(i, k, l,PF), i = 1:n, k = 1:t , l = 1:s)
  
  
   result$constraints[[length(result$constraints)-250]]
    
    # Objective Function
    result <- result %>% 
      set_objective(sum_expr(RHO[l,1] * (P[i,1] * g[i, k] + 100*b[i, k, l]), i = 1:n, k = 1:t, l=1:s),
                    sense = "min") 
  

    result$constraints[[length(result$constraints)-200]]
    
    result <- result %>% 
      solve_model(with_ROI("glpk", verbose = TRUE))
    
  return(result)
}


rgen(5,10,2,PF2)

### Simulation Functions

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

rgen <- function(i, t, s,PF2) {
  vin <- 3.5
  vr <-  14.5
  vout <- 20
  pw <- 5
  rgenSet <- 4:5
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

dgen <- function(i,t,s,PF){
  demandSet <- 1:3
  if(i %in% demandSet){
    d <- as.matrix(final.scenarios[s,] %>%  dplyr::select(load1,load2,load3),ncol=3)
    mm <- t(PF %*% d)
    return(as.numeric(mm[i,t]))
  }else
    0
  
}

