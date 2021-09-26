# 
rm(list = setdiff(ls(), lsf.str()))
gc()




actual.params <- list()
ms <-  c(3,2,4)+1
sds <- c(1,1.5,0.5)
shapes <- c(8,7)
scales <- c(3.5,4)
assignActualParams(ms,sds,shapes,scales)

nn <- 5
load1 <- rnorm(250,ms[1],sds[1])
load2 <- rnorm(250,ms[2],sds[2])
load3 <- rnorm(250,ms[3],sds[3])
wspeed1 <- rweibull(250,shape= shapes[1],scale=scales[1]) 
wspeed2 <- rweibull(250,shape= shapes[1],scale=scales[2])
init.df <- data.frame(load1,load2,load3,wspeed1,wspeed2)

demandSet <- 1:3
rgenSet <- 4:5
step_size <- 1




maxIter <- 10
f.result.final <-  as.data.frame(matrix(ncol = 3))
colnames(f.result.final) <- c("optimization","simulation","lat")
lats <-  c(seq(1000,50,-250))

for(lat in lats){
  
  demandL <- list()
  rgenL <-  list()
  f.result <- as.data.frame(matrix(0,ncol=3))
  colnames(f.result) <- c("optimization","simulation","lat")
  
  for(iter in 1:maxIter){
    
    ms.n <- ms * runif(1,1,1.35)
    sds.n <- sds * runif(1,1,1.35)
    shapes.n <- shapes * runif(1,0.95,1.05)
    scales.n <- scales * runif(1,0.95,1.05)
    assignActualParams(ms.n,sds.n,shapes.n,scales.n)
    
    r.df <- createRealTimeData(means = ms.n,sds = sds.n,
                               shapes = shapes.n, scales = scales.n,
                               lat)
    
    # r.df.bs <- do.call("rbind", replicate(2, r.df, simplify = FALSE))
    
    upd.est <- updatedEstimates2(r.df)
    assignEstimates(upd.est)
    
    scenarios <- generateScenarios(demandL,rgenL) 
    reduced.scenes <- reduceScenarios(scenarios)
    
    result <- solveOPT3(reduced.scenes,step_size = 1)
    result.df <- calculateB(result)
    result.df2 <- furtherSimulate(result.df, result = result)
    
    
    f.result[iter,] <- c(result$objective_value,sum(result.df2$cost),lat)
     
  }
  
  f.result.final <- rbind(f.result.final,f.result)
}


