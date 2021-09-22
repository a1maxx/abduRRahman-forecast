rm(list = ls())
gc()
rm(list = setdiff(ls(), lsf.str()))


library(ggplot2)
library(ROI)
library(ROI.plugin.glpk)
library(ompr.roi)
library(ompr)
library(magrittr)
library(Rglpk)
library(fitdistrplus)
library(cluster)
library(reshape2)
# detach("package:simmer", unload = TRUE)
# library(Rtsne)
library(Rlab)
library(arpr)
suppressPackageStartupMessages(library(dplyr))
options(dplyr.summarise.inform = FALSE)

createRealTimeData<- function(means, sds,shapes,scales,n){
  dl <- length(ms)
  rl <- length(shapes)
  dl.names <- paste("load",1:dl,sep="")
  for(i in 1:length(dl.names)){
    assign(dl.names[i],rnorm(n,means[i],sds[i]))
    
  }
  
  rl.names <- paste("wspeed",1:rl,sep="")
  for(i in 1:length(rl.names)){
    assign(rl.names[i],rweibull(n,shapes[i],scales[i]))
    
  }
  
  ls <- ifelse(sapply(dl.names,function(x) get(x,envir = parent.frame()))<0,0,sapply(dl.names,function(x) get(x,envir = parent.frame())))
  rs <- ifelse(sapply(rl.names,function(x) get(x,envir = parent.frame()))<0,0,sapply(rl.names,function(x) get(x,envir = parent.frame())))
  
  df <- data.frame(cbind(ls,rs))
  
}

assignActualParams <- function(means, sds,shapes,scales,env=parent.env()){
  nofD <- length(means)
  demand.mean.actual <- means
  demand.sd.actual <- sds
  demand.actual.params <- vector(mode="list",length=nofD)
  for(i in 1:nofD) {
    demand.actual.params[[i]]$mean <- demand.mean.actual[i]
    demand.actual.params[[i]]$sd <- demand.sd.actual[i]
    
  }
  nofR <- length(shapes)
  rgen.shape.actual <- shapes
  rgen.scale.actual <- scales
  rgen.actual.params <- vector(mode="list",length=nofR)
  for( i in 1:nofR){
    rgen.actual.params[[i]]$shape <- rgen.shape.actual[i]
    rgen.actual.params[[i]]$scale <- rgen.scale.actual[i]
    
  }
  
  actual.params <<- append(demand.actual.params,rgen.actual.params)
}

createDatawithSampleParamEst <-  function(l,d,r,n){
  
  dl <- l[1:d]  
  rl <- l[(d+1):(d+r)]
  dl.names <- paste("load",1:length(dl),sep="")
  for(i in 1:length(dl.names)){
    assign(dl.names[i],rnorm(n,dl[[i]]$mean.est,dl[[i]]$sd.est))
    
  }
  
  rl.names <- paste("wspeed",1:length(rl),sep="")
  for(i in 1:length(rl.names)){
    assign(rl.names[i],rweibull(n,shape = rl[[i]]$shape.est,rl[[i]]$scale.est))
    
  }
  
  ls <- ifelse(sapply(dl.names,function(x) get(x,envir = parent.frame()))<0,0,sapply(dl.names,function(x) get(x,envir = parent.frame())))
  rs <- ifelse(sapply(rl.names,function(x) get(x,envir = parent.frame()))<0,0,sapply(rl.names,function(x) get(x,envir = parent.frame())))
  
  df <- data.frame(cbind(ls,rs))
  
  
}

generateEstimates <- function(df){
  sink(nullfile())
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
  sink()
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

updatedEstimates <- function(archival,new){
  df <-  rbind(archival,new)
  upd.m <- generateEstimates(df)
  
  return(upd.m)
}

generateScenarios <- function(demandL,rgenL){
  
  scenarios <- data.frame(sn=0,load1=0,load2=0,load3=0,wspeed1=0,wspeed2=0,prob=0)
  
  for (i in 1:1000) {
    
    loads <- c()
    rweis <- c()
    
    for (j in 1:length(demandL))
      loads <- append(loads, rnorm(1, mean = as.numeric(demandL[[j]][1]), sd = as.numeric(demandL[[j]][2])))
    
    
    for (k in 1:length(rgenL)) 
      rweis <- append(rweis,rweibull(1, shape = as.numeric(rgenL[[k]][1]) , scale = as.numeric(rgenL[[k]][2])))
    
    
    load.prob <- 1
    
    for (r in 1:length(loads))
      load.prob <- dnorm(loads[r], mean = as.numeric(demandL[[r]][1]), sd = as.numeric(demandL[[r]][2])) * load.prob
    
    
    rgen.prob <- 1
    
    for (p in 1:length(rweis))
      rgen.prob <- dweibull(rweis[p], shape = as.numeric(rgenL[[p]][1]), scale = as.numeric(rgenL[[p]][2])) * rgen.prob
    
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

solveOPT <- function(final.scenarios,step_size,env=emptyenv()){
  n <<- 5
  
  rgen <- function(i, t, s, PF2) {
    vin <- 3.5
    vr <-  14.5
    vout <- 20
    pw <- 5
    rgenSet <- 4:5
    if (i %in% rgenSet) {
      d <-
        as.matrix(final.scenarios[s,] %>%  dplyr::select(wspeed1, wspeed2),
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
  
  dgen <- function(i, t, s, PF) {
    demandSet <- 1:3
    if (i %in% demandSet) {
      d <-
        as.matrix(final.scenarios[s, ] %>%  dplyr::select(load1, load2, load3),
                  ncol = 3)
      mm <- t(PF %*% d)
      return(as.numeric(mm[i, t]))
    } else
      0
    
  }


  t <- step_size
  s <- nrow(final.scenarios)
  line_cap <- 7.5
  A <-  matrix(0,ncol=n,nrow=n) # if there is connection between i and j
  A[lower.tri(A)] <- c(1,0,1,1,0,1,0,1,0,1)
  A[upper.tri(A)] = t(A)[upper.tri(A)]
  C <-  matrix(sapply(diag(n),function(x)ifelse(x==1,0,line_cap)),ncol = n) # Capacity of line i-j 
  G <-  matrix(c(40,0,0,0,20),ncol=n) # Generation capacity of bus i
  P <-  matrix(c(4,0,0,0,8),nrow=n,ncol=1) # Price of generating unit of electricity in bus i
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
  
   
    # Objective Function
    result <- result %>% 
      set_objective(sum_expr(RHO[l,1] * (P[i,1] * g[i, k] + 10*b[i, k, l]), i = 1:n, k = 1:t, l=1:s),
                    sense = "min")
  
    result <- result %>% 
      solve_model(with_ROI("glpk", verbose = TRUE))
  
  return(result)
}

calculateB <- function(result,env = parent.env()){
  df <- as.data.frame(matrix(ncol = 8))
  colnames(df) <- c("node","step","arrivals","departures","generation","demand","rgen","excessORdeficit")
  gg <- get_solution(result,g[i,k])
  xx <- get_solution(result,x[i,j,k])
  N <- max(get_solution(result,x[i,j,k])$i)
  P <-  max(get_solution(result,x[i,j,k])$k)
  rij <- get_solution(result,r[i,j]) %>% filter(value>0)
  prices <- matrix(c(4,0,0,0,8),nrow=N,ncol=1) 
  
  for (p in 1:P) {
    for (y in 1:N) {
      arrivals <- sum(xx %>% filter(j == y, k == p) %>% dplyr::select(value))
      departures <- sum(xx %>% filter(i == y, k == p) %>% dplyr::select(value))
      generation <- as.numeric(gg %>% filter(i==y,k==p) %>% dplyr::select(value))
      
      dd <-ifelse(y %in% demandSet,
                  rnorm(1, actual.params[[y]]$mean, actual.params[[y]]$sd) * runif(1, 0.9, 1.1),0)
      
      rr <- ifelse(y %in% rgenSet,
                   rweibull(1, actual.params[[y]]$shape, actual.params[[y]]$scale) * runif(100, 0.95, 1.05),0) 
      
      bb <- arrivals + generation + rr - dd - departures
      
      
      df <-  rbind(df,c(y,p,arrivals,departures,generation,dd,rr,bb))
      
    }
    
  }

  
  df <- df[-1,] %>% mutate(genCost = generation * prices[node,])
  
  return(df)
}

furtherSimulate <- function(result, result.df) {
  
  defs <-
    as.numeric(unlist(result.df %>% filter(excessORdeficit < 0) %>% dplyr::select(node)))
  excs <-
    as.numeric(unlist(result.df %>% filter(excessORdeficit > 0) %>% dplyr::select(node)))
  
  xijk <- get_solution(result,x[i,j,k]) %>% 
    filter(value>0) %>% dplyr::select(-variable)
    
  t.df <- xijk %>% group_by(across(all_of(c("i","j")))) %>% 
    summarise(sumTrans= sum(value)) %>% mutate(ij = paste(i,j,sep=""))

  
  my_df <- result.df
  i <- 4
  j <- 2
  
  for (i in excs) {
    for (j in defs) {
      line <- paste(i,j,sep="")
      if(line %in% t.df$ij){
        row.n <- which(t.df$ij == line)
        if(as.logical(t.df[row.n,"sumTrans"] < 6)){
          nj.df <- which(my_df$node==j)
          ni.df <- which(my_df$node==i)
          possible_transfers <- c()
          possible_transfers[1] <- -my_df[nj.df,"excessORdeficit"]
          possible_transfers[2] <- as.numeric(6 - t.df[row.n,"sumTrans"])
          possible_transfers[3] <- my_df[ni.df,"excessORdeficit"]
          final_transfer <- min(possible_transfers)
          
          my_df[nj.df,"arrivals"] = my_df[nj.df,"arrivals"] + final_transfer
          my_df[nj.df,"excessORdeficit"] = my_df[nj.df,"excessORdeficit"] + final_transfer
          
          my_df[ni.df,"departures"] = my_df[ni.df,"departures"] + final_transfer
          my_df[ni.df,"excessORdeficit"] = my_df[ni.df,"excessORdeficit"] - final_transfer
  
        }
        
      }
       
    }
  
  }
  
  my_df <- my_df %>% mutate(cost = ifelse(excessORdeficit<0,-excessORdeficit*40+genCost,excessORdeficit*(2.5) + genCost))
  

 return(my_df)  
}

solveOPT2 <- function(final.scenarios,step_size,env=emptyenv()){
  
  rgen <- function(i, t, s, PF2) {
    vin <- 3.5
    vr <-  14.5
    vout <- 20
    pw <- 5
    rgenSet <- 4:5
    if (i %in% rgenSet) {
      d <-
        as.matrix(final.scenarios[s,] %>%  dplyr::select(wspeed1, wspeed2),
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
  
  dgen <- function(i, t, s, PF) {
    demandSet <- 1:3
    if (i %in% demandSet) {
      d <-
        as.matrix(final.scenarios[s, ] %>%  dplyr::select(load1, load2, load3),
                  ncol = 3)
      mm <- t(PF %*% d)
      return(as.numeric(mm[i, t]))
    } else
      0
    
  }

  n <- 5
  t <- step_size
  s <- nrow(final.scenarios)
  line_cap <- 7.5
  A <-  matrix(0,ncol=n,nrow=n) # if there is connection between i and j
  A[lower.tri(A)] <- c(1,0,1,1,0,1,0,1,0,1)
  A[upper.tri(A)] = t(A)[upper.tri(A)]
  C <-  matrix(sapply(diag(n),function(x)ifelse(x==1,0,line_cap)),ncol = n) # Capacity of line i-j 
  G <-  matrix(c(40,0,0,0,20),ncol=n) # Generation capacity of bus i
  P <-  matrix(c(4,0,0,0,8),nrow=n,ncol=1) # Price of generating unit of electricity in bus i
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
  
  # Objective Function
  result <- result %>%
    set_objective(sum_expr(
      RHO[l, 1] * (P[i, 1] * g[i, k] + 10 * b[i, k, l]),
      i = 1:n,
      k = 1:t,
      l = 1:s
    ) -  5 * (
      sum_expr(
        b[i, k, l] + sum_expr(x[j, i, k], j = 1:n, i != j) +
          g[i, k] + rgen(i, k, l, PF2) - dgen(i, k, l, PF) - sum_expr(x[i, j, k], j = 1:n, 
                                                                      i != j),
        i = 1:n,
        k = 1:t,
        l = 1:s
      )
    ),
    sense = "min")
  
  result <- result %>% 
    solve_model(with_ROI("glpk", verbose = TRUE))
  
  return(result)
}

solveOPT3 <- function(final.scenarios,step_size,env=emptyenv()){
  n <- 5
  rgen <- function(i, t, s, PF2) {
    vin <- 3.5
    vr <-  14.5
    vout <- 20
    pw <- 5
    rgenSet <- 4:5
    if (i %in% rgenSet) {
      d <-
        as.matrix(final.scenarios[s,] %>%  dplyr::select(wspeed1, wspeed2),
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
  
  dgen <- function(i, t, s, PF) {
    demandSet <- 1:3
    if (i %in% demandSet) {
      d <-
        as.matrix(final.scenarios[s, ] %>%  dplyr::select(load1, load2, load3),
                  ncol = 3)
      mm <- t(PF %*% d)
      return(as.numeric(mm[i, t]))
    } else
      0
    
  }
  nn <- 5
  t <- step_size
  s <- nrow(final.scenarios)
  line_cap <- 6.0
  A <-  matrix(0,ncol=nn,nrow=nn) # if there is connection between i and j
  A[lower.tri(A)] <- c(1,0,1,1,0,1,0,1,0,1)
  A[upper.tri(A)] = t(A)[upper.tri(A)]
  C <-  matrix(sapply(diag(nn),function(x)ifelse(x==1,0,line_cap)),ncol = nn) # Capacity of line i-j 
  G <-  matrix(c(5,0,0,0,2.5),ncol=nn) # Generation capacity of bus i
  P <-  matrix(c(20,0,0,0,40),nrow=nn,ncol=1) # Price of generating unit of electricity in bus i
  RHO <- matrix(final.scenarios$prob,nrow = nrow(final.scenarios)) # Scenario Probabilities
  PF <- matrix(sample(runif(100,0.9,1.1),step_size),nrow = step_size)  # Further Randomization of demand
  PF2 <- matrix(sample(runif(100,0.95,1.05),step_size),nrow = step_size) # Further Randomization of renewable generation
  
  result <- MIPModel() %>% 
    # Amount of electricity transferred from i to j at time k
    add_variable(x[i, j, k], i = 1:nn, j = 1:nn, k = 1:t,  type = "continuous",lb = 0) %>% 
    
    # Amount of electricity generated at i at time k  
    add_variable(g[i, k],
                 i = 1:nn, k = 1:t, type = "continuous",lb = 0) %>% 
    
    # If line i-j selected to be the path to transfer electricity from i to j
    add_variable(r[i, j], 
                 i = 1:nn, j = 1:nn, type = "binary") %>% 
    
    # Amount of electricity supplied from main grid
    add_variable(b[i, k], i = 1:nn, k = 1:t, type = "continuous",lb=0) %>% 
    
    # Line binding constraint
    add_constraint(x[i, j, k] <= r[i,j]* C[i,j] * A[i,j],
                   i = 1:nn, j = 1:nn, k = 1:t, l = 1:s) %>% 
    
    # Generation capacity constraint 
    add_constraint(g[i, k] <= G[1,i], 
                   i = 1:nn , k = 1:t)  
    
    # Balance constraint
    result <- result %>% add_constraint(sum_expr(x[i, j, k], j = 1:nn, i!=j) <= b[i,k] + sum_expr(x[j, i, k], j = 1:nn, i!=j) + 
                     g[i, k] + rgen(i, k, l,PF2) - dgen(i, k, l,PF), i = 1:nn, k = 1:t , l = 1:s) %>% 
    
    add_constraint(b[i, k] <= 10, i = 1:nn, k = 1:t)
    
    
    # Objective Function
    result <- result %>%
      set_objective(
        sum_expr(
          RHO[l, 1] * (P[i, 1] * g[i, k] + 40 * b[i, k]) - 
            2.5 * ( b[i, k] + sum_expr(x[j, i, k], j = 1:nn, i != j) + g[i, k] + rgen(i, k, l, PF2) - dgen(i, k, l, PF) - sum_expr(x[i, j, k], j = 1:nn, i != j) ) ,
              i = 1:nn,
              k = 1:t,
              l = 1:s
        ),
        sense = "min"
      )
     
    
    sink(nullfile())
    result <- result %>% solve_model(with_ROI("glpk", verbose = TRUE))
    
    sink()
  return(result)
}



result$objective_value
get_solution(result,b[i,k])
get_solution(result,g[i,k])


triel <- get_solution(result,x[i,j,k])
triel <- triel  %>% select(-variable,-k)
spread(triel,key=-j,value =value,drop = T) %>% select(-i) 
