library(tsoutliers)
library(TSclust)
library(factoextra)
library(NbClust)
library(reshape2)
library(dplyr)
library(minpack.lm)
library(TSdist)
library(ggthemes)
library(rminer)
library(ggplot2)
library(ggsci)

 # split <- holdout(sogoods, 0.9)
 split2 <- holdout(sogoods,0.85)

f = function(p, q, t) {
  res = (exp((p + q) * t) * p * (p + q) ^ 2) / (p * exp((p + q) * t) + q) ^
    2
  
}
FF <- expression(p * (exp((p + q) * t) - 1) / (p * exp((p + q) * t) + q))
ff <- function(p, q, t) {
  res <-  D(FF, "t")
}
PLC.ts <- ts(deltat = 1 / 52)





for (i in sogoods[split2$tr]) {
  data <- Product[[i]]@sales
  
  startY <- as.numeric(format(data$week[1], "%y"))
  startM <- as.numeric(format(data$week[1], "%m"))
  data.ts <-
    ts(
      data = data$quant,
      start = c(startY, startM),
      deltat = 1 / 52,
      names = c("week", "sales")
    )
  
  
  
 #### OUTLIERS DETECTION and REPLACEMENT
  
  
  outliers <- tso(data.ts,
                  types = "AO",
                  maxit.iloop = 20,
                  cval = 10.0)
  
  while (length(outliers$outliers$ind) != 0) {
    outliers <- tso(data.ts,
                    types = "AO",
                    maxit.iloop = 20,
                    cval = 10.0)
    
    if (length(outliers$outliers$ind) != 0) {
      k1 <- 1:nrow(data.ts)
      t <- outliers$outliers$ind[1]
      k1 <- k1[-(t:nrow(data.ts))]
      term1  <- k1 * data.ts[k1]
      k2 <- (nrow(data.ts) - t):1
      indice <- (t + 1):nrow(data.ts)
      term2 <- k2 * data.ts[indice]
      correction <- (sum(term1) + sum(term2)) / (sum(k1) + sum(k2))
      data.ts[t] <- correction
    }
    
  }
  #### Normalization of data
  normalized.ts <- data.ts / sum(data.ts)
  
  ### FITTING BASS MODEL (NLS)
  mnlS <- sum(normalized.ts) * 0.75
  pnlS <- 0.01
  qnlS <- 0.1
  t <- 1:length(normalized.ts)
  control <-
    nls.control(
      maxiter = 1024,
      minFactor = 1 / 4096 * 2,
      printEval = FALSE,
      warnOnly = TRUE
    )
  # Bass.nls<-nls(normalized.ts ~ M*(((P+Q)^2/P)*exp(-(P+Q)*t))/(1+(Q/P)*exp(-(P+Q)*t))^2,
  #               start=c(list(M=mnlS,P=pnlS,Q=qnlS)),control = control,algorithm = "port")
  Bass.nls <-
    nlsLM(
      normalized.ts ~ M * (((P + Q) ^ 2 / P) * exp(-(P + Q) * t)) / (1 + (Q /
                                                                            P) * exp(-(P + Q) * t)) ^ 2,
      start = c(list(
        M = mnlS, P = pnlS, Q = qnlS
      )),
      control = control,
      algorithm = "LM",
      lower = c(0, 0, 0),
      upper = c(Inf, 1, 1)
    )
  
  m <- coef(Bass.nls)[1]
  p <- coef(Bass.nls)[2]
  q <- coef(Bass.nls)[3]
  
  
  ####
  
  t <- 1:30
  PLC.data <- f(p, q, t) * m
  
  PLCs <-
    ts(
      data = PLC.data,
      start = c(1, 1),
      deltat = 1 / 52,
      names = c("week", "sales")
    )
  
  # PLCs <- PLCs / sum(PLCs)
  
  PLC.ts <- rbind(PLC.ts, PLCs)
  
  
}

PLC.ts <- PLC.ts[-1, ]
PLC.ts <- ts(t(PLC.ts), start = c(1, 1), deltat = 1 / 52)

colnames(PLC.ts) <- sogoods[split2$tr]
# PLC.ts<-scale(PLC.ts)
distanze <- diss(t(PLC.ts), "CORT",k=3,deltamethod="Euclid")
# distanze <- diss(PLC.ts, "CORT")

# distanze2 <- diss(PLC.ts, "CORT")
# distanze2 <- diss(PLC.ts, "DTWARP")
# distanze3 <- diss(PLC.ts, "DWT")
# distanze4 <- TSDatabaseDistances(PLC.ts, distance = "cort")


# wss <- 0
# for (i in 2:15) {
# km <- pam(x = distanze,
# k = i,
# )
# wss[i] <- km$silinfo$avg.width
# }
# 
# plot(
#   1:15,
#   wss,
#   type = "b",
#   col = "grey20",
#   main = "Scree Plot",
#   ylab = "Toplam Küme İçi Uzaklık",
#   xlab = "Küme Sayısı"
# )


##PAM clustering
pam.res <- pam(x = distanze, k = 4, diss = inherits(distanze,"dist"))
PLC.df <- as.data.frame(PLC.ts)
PLC.df[31, ] <- pam.res$clustering
#PAM END

 plot(pam.res)


plot(PLC.ts[,pam.res$medoids[1]],ylim=c(0,0.1),xlab="Time Periods",ylab="Normalized Values",main="Representative Curves")
lines(PLC.ts[,pam.res$medoids[2]],col="blue")
lines(PLC.ts[,pam.res$medoids[3]],col="green")
lines(PLC.ts[,pam.res$medoids[4]],col="red")
legend("topright", legend = c(paste("Center",1:4)),
       col = c("black", "blue","green","red"), lty = 1, cex = 0.8) 


# lines(PLC.ts[,pam.res$medoids[5]],col="darkgreen")
# lines(PLC.ts[,pam.res$medoids[6]],col="grey")

PLC.mat <- as.matrix(PLC.df)
PLC.mat <- t(PLC.mat)
PLC.df <- as.data.frame(PLC.mat)
colnames(PLC.df)[31] <- "clusters"

plcdf2 <-
  PLC.df %>% mutate(clusters = as.factor(paste("Cluster", clusters)))


df <- melt(plcdf2)  #the function melt reshapes it from wide to long
df$rowid <- 1:length(sogoods[split2$tr])
df
ggplot(df, aes(variable, value, group = factor(rowid),color=clusters)) +
  geom_line( alpha = 0.2) + 
  facet_wrap( ~ clusters) + theme_Publication() + scale_color_aaas() +
  labs(x="Time Periods",y="Normalized Values") + ggtitle("Clusters") + coord_cartesian(ylim = c(0,0.2))




##KMeans Clustering
# set.seed(12324)
# km.res <- kmeans(x = distanze, 6, nstart = 20)
# PLC.df <- as.data.frame(PLC.ts)
# PLC.df[31, ] <- km.res$cluster
##KmeansEND

# ##Kmedoids Start
# kmed.res <- KMedoids(PLC.ts, 6, "cort")
# PLC.df <- as.data.frame(PLC.ts)
# PLC.df[31, ] <- kmed.res
# #Kmedoids END







silhoutte.err <-  c()
for (k in 2:15) {
  kmed.res <- pam(x = distanze, k , diss = inherits(distanze,"dist"))
  sil.err <- silhouette(kmed.res, distanze)
  silhoutte.err <- append(silhoutte.err, mean(sil.err[, 3]))
}
a <- c(0,silhoutte.err)
ggplot()+geom_line(aes(x=1:15,y=a)) +theme_Publication() + scale_color_aaas() + 
  labs(y="Silüet genişliği",x="Küme sayısı",size="Ortalama genişlik") + 
  geom_point(aes(x=2:15,y=a[-1],size=a[-1])) 



# datam <- as.matrix(t(PLC.ts))
# # indizes <- c("frey", "mcclain", "cindex", "silhouette", "dunn")
# indizes <- c(
#   "kl",
#   "ch",
#   "hartigan",
#   "cindex",
#   "db",
#   "silhouette",
#   "pseudot2",
#   "duda",
#   "ratkowsky",
#   "ball",
#   "ptbiserial",
#   "tau",
#   "dunn",
#   "hubert",
#   "sdindex",
#   "sdbw"
# )
# 
# for (i in indizes) {
#   nb <-
#     NbClust(
#       data = datam,
#       diss = distanze,
#       distance = NULL,
#       min.nc = 3,
#       max.nc = 20,
#       method = "centroid",
#       alphaBeale = 0.1,
#       index = i
#     )
#   
#   print(i)
#   print(nb$Best.nc)
#   
# }
# 
# 
# fviz_nbclust(nb, method = "silhouette", pam, diss = distanze4)
# fviz_nbclust(
#   t(PLC.ts),
#   method = "silhouette",
#   diss = distanze4,
#   kmeans,
#   k.max = 10,
#   nboot = 5
# )
# 
# 
fviz_nbclust(
  PLC.ts,
  method = "wss",
  diss = distanze,
  pam,
  k.max = 10,
  nboot = 5
)

i<- sogoods[380]
iteras <-1
for(i in sogoods)
  
{print(paste(iteras,Product[[i]]@evaluationNLS[Product[[i]]@nrow,2]))
  iteras <- iteras +1
}

nrowa <- Product[[i]]@nrow
pa <- Product[[i]]@evaluation2[nrowa,2]
qa <- Product[[i]]@evaluation2[nrowa,3]
ma <- Product[[i]]@evaluation2[nrowa,4]
ta <- 1:Product[[i]]@nrow  

sum(Product[[i]]@sales$quant)
plot(1:nrowa,f(pa,qa,ta)*ma,ylim=c(0,100))
lines(1:nrowa,Product[[i]]@sales$quant)  
