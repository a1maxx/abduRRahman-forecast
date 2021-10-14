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


maxIter <- 96
f.result.final <-  as.data.frame(matrix(ncol = 5))
colnames(f.result.final) <- c("iter","optimization","simulation","lat","distmethod")

meths <- c("euclidean","manhattan")
lats <- c(3,5,50)
for(meth in meths){

  for(iter in 1:maxIter){
    
    ms.n <- ms * runif(1,1,2)
    sds.n <- sds * runif(1,1,2)
    shapes.n <- shapes * runif(1,0.9,1.1)
    scales.n <- scales 
    assignActualParams(ms.n,sds.n,shapes.n,scales.n)
    
    r.df <- createRealTimeData(means = ms.n,sds = sds.n,
                               shapes = shapes.n, scales = scales.n,
                               max(lats))
    
    f.result <- as.data.frame(matrix(0,ncol=5))
    colnames(f.result) <- c("iter","optimization","simulation","lat","distmethod")
    
    dds <- sapply(1:5, function(x) if(x %in% demandSet) 
      rnorm(1,actual.params[[x]]$mean,actual.params[[x]]$sd) else 0 )
    
    rrs <- sapply(1:5, function(x) if(x %in% rgenSet) 
      rweibull(1,actual.params[[x]]$shape,actual.params[[x]]$scale) else 0 )
    
    
    for(lat in lats){
      
      demandL <- list()
      rgenL <-  list()
      
      # r.df.bs <- do.call("rbind", replicate(2, r.df, simplify = FALSE))
      
      r.df2 <- r.df[1:lat,] 
      
      upd.est <- updatedEstimates2(r.df2)
      
      demandL <- upd.est[demandSet]
      rgenL <- upd.est[rgenSet]
        
      scenarios <- generateScenarios(demandL,rgenL)
      
      reduced.scenes <- reduceScenarios(scenarios,meth)
      
      result <- solveOPT4(reduced.scenes,step_size = 1)
      
  
      result.df <- calculateB3(result,reduced.scenes$prob,dds,rrs)
      result.df2 <- furtherSimulate2(result.df, result = result, reduced.scenes$prob)
      f.result <- rbind(f.result,c(iter,result$objective_value,sum(result.df2$cost),
                                   lat,meth))
  
      
      
      
    }
    
    f.result.final <- rbind(f.result.final,f.result[-1,])
   
    
  }

}
