
rm(list = setdiff(ls(), lsf.str()))
gc()

load1 <- rnorm(250,3,1)
load2 <- rnorm(250,2.5,1.5)
load3 <- rnorm(250,4,0.5)
load1 <- ifelse(load1<0,0,load1)
load2 <- ifelse(load2<0,0,load2)
load3 <- ifelse(load3<0,0,load3)
wspeed1 <- rweibull(250,shape= 8,scale=3.5) 
wspeed2 <- rweibull(250,shape= 7,scale=4)
demandSet <- 1:3
rgenSet <- 4:5
step_size <- 1


actual.params <- list()
ms <-  c(3,2,4)
sds <- c(1,1.5,0.5)
shapes <- c(8,7)
scales <- c(3.5,4)
assignActualParams(ms,sds,shapes,scales)
init.df <- data.frame(load1,load2,load3,wspeed1,wspeed2)


demandL <- list()
rgenL <-  list()
maxIter <- 5

for(iter in 1:maxIter){
  
  ms <- ms * runif(1,1,1.5)
  sds <- sds * runif(1,1,1.5)
  shapes <- shapes * runif(1,1,1.5)
  scales <- scales * runif(1,1,1.5)
  assignActualParams(ms,sds,shapes,scales)
  
  r.df <- createRealTimeData(ms,sds,shapes,scales,1000)
  
  upd.est <- updatedEstimates(init.df,r.df)
  
  assignEstimates(upd.est)
  scenarios <- generateScenarios(demandL,rgenL)
  reduced.scenes <- reduceScenarios(scenarios)
  result <- solveOPT2(reduced.scenes,step_size = 1)
  
  result.df <- calculateB(result)
  
  result.df
  
  
  get_solution(result,r[i,j]) %>% dplyr::filter(value>0)
  
  result$objective_value
  agg.df <-
    result.df %>% group_by(step) %>% dplyr::summarise(
      utilityCostbyStep = sum(utilityCost),
      generationCostbyStep = sum(genCost),
      totalSum = sum(utilityCostbyStep, generationCostbyStep)
    ) 
  print(result$objective_value)
  print(sum(agg.df$totalSum))
  
  
}


