library(fitdistrplus)
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


actual.params <- list()
ms <-  c(3,2,4)
sds <- c(1,1.5,0.5)
shapes <- c(8,7)
scales <- c(3.5,4)
## Changes actual parameter estimations to be used in the simulation
assignActualParams(ms,sds,shapes,scales)


init.df <- data.frame(load1,load2,load3,wspeed1,wspeed2)
init.ests <- generateEstimates(init.df)
rgenL <-  list()
demandL <- list()
assignEstimates(init.ests)
scenarios <- generateScenarios(demandL,rgenL)
reduced.scenes <- reduceScenarios(scenarios)
result <- solveOPT(final.scenarios = reduced.scenes,step_size = 10)
result.df <- calculateB(result)
