# 
rm(list = setdiff(ls(), lsf.str()))
gc()




actual.params <- list()
ms <-  c(5,4,6)
sds <- c(1,1.5,0.5)
shapes <- c(8,7)
scales <- c(3.5,4)
assignActualParams(ms,sds,shapes,scales)


nn <- 5
demandSet <- 1:3
rgenSet <- 4:5
step_size <- 1


maxIter <- 20
f.result.final <-  as.data.frame(matrix(ncol = 4))
colnames(f.result.final) <- c("iter","optimization","simulation","lat")

lats <- c(5,100,200)
repS <- 1:10

for(iter in 1:maxIter){
  
  ms.n <- ms * runif(1,1,2)
  sds.n <- sds * runif(1,1,2)
  shapes.n <- shapes * runif(1,0.995,1.005)
  scales.n <- scales * runif(1,0.995,1.005)
  assignActualParams(ms.n,sds.n,shapes.n,scales.n)
  
  r.df <- createRealTimeData(means = ms.n,sds = sds.n,
                             shapes = shapes.n, scales = scales.n,
                             max(lats))
  
  f.result <- as.data.frame(matrix(0,ncol=4))
  colnames(f.result) <- c("iter","optimization","simulation","lat")
  
  
  
  for(lat in lats){
    
    demandL <- list()
    rgenL <-  list()

    # r.df.bs <- do.call("rbind", replicate(2, r.df, simplify = FALSE))
    
    r.df2 <- r.df[1:lat,] 
    
    upd.est <- updatedEstimates2(r.df2)
    
    assignEstimates(upd.est)
    
    scenarios <- generateScenarios(demandL,rgenL)
    
    reduced.scenes <- reduceScenarios(scenarios)
    
    result <- solveOPT4(reduced.scenes,step_size = 1)
    
    
    
    for(reps in repS){
      
      result.df <- calculateB2(result,reduced.scenes$prob)
      result.df2 <- furtherSimulate2(result.df, result = result, reduced.scenes$prob)
      f.result <- rbind(f.result,c(iter,result$objective_value,sum(result.df2$cost),
                           lat))
    }
    
    
  
  }
  
  f.result.final <- rbind(f.result.final,f.result[-1,])
  
}

