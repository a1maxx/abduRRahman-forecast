library(ROI)
library(ROI.plugin.glpk)
library(ompr.roi)
library(ompr)
library(magrittr)
library(Rglpk)
library(dplyr)
#
n <- 5
t <- 10
a <- matrix(c(0,1,1,
         1,0,1,
         1,1,0),ncol=3)

set.seed(1)

weights <- matrix(rpois(n * n, 5), ncol = n, nrow = n)
demand  <- c(rnorm(n,0,1))
demand  <- matrix(sapply(demand,function(x) ifelse(x<0,0,x)),nrow=1)
       


result <- MIPModel() %>% 
  add_variable(x[i, j, k], i = 1:n, j = 1:n, k=1:t, type = "continuous") %>% ## Amount of electricity transferred from i to j at time k
  add_vaiable(g[i,k],i=1:n,k=1:t,type = "continuous") %>% ## Amount of electricity generated at i at time k  
  
  set_objective(sum_expr(weights[i, j] * x[i, j, k], i = 1:n, j = 1:n,k=1:n)) %>% 
  add_constraint( sum_expr(x[i, j, k], k = 1:n) <= 1, i = 1:n,j=1:n) %>% 
  add_constraint(x[i,j,k]<=5,i=1:n,j=1:n,k=1:n) %>% ß
  solve_model(with_ROI("glpk", verbose = TRUE))

get_solution(result, x[i, j, k]) %>% filter(value>0,variable=="x",i==1)

grepl('livd','alive')

sort(c('a','g','b'))

a <- function(a1, a2) {
  temp <- c()
  for (i in a1) {
    for (j in a2) {
      if (grepl(i,j)) {
        temp <- append(temp, i)
        a1 <- a1[-which(a1==i)]
        break
      }
    }

  }
  if(length(temp)>0)
    return(sort(temp))
  else 
    return(c())
}



a1 = c("arp", "arp", "bull") 
a2 = c("lively", "alive", "harp", "sharp", "armstrong")
a(a1,a2)
rm(list=ls())
lng <-5
wd <- 3
squaresInRect <- function(lng, wd) {
  totalSquares <- lng*wd
  maxit <- min(lng,wd)
  vec <- c()
  if(lng!=wd){
    while(sqrt(totalSquares)>1){
      item <- min(floor(sqrt(totalSquares)),min(lng,wd))
      vec <- append(vec,item)
      totalSquares <- totalSquares - item**2
      lng <- ifelse(lng-item==0,lng,lng-item)
      wd <- ifelse(wd-item==0,wd,wd-item)
    }
    
    return(append(vec,rep(1,totalSquares)))
  }else{
    return(NULL)
  }
}
sqrt(280)
sentence <- "is2 Thi1s T4est 3a"

fsplit <-  strsplit(sentence," ")

unlist(strsplit(fsplit[[1]][1],""))
as.numeric(unlist(strsplit(fsplit[[1]][1],""))[grep("[0-9]+",unlist(strsplit(fsplit[[1]][1],"")))])

sapply(fsplit, function(x) length(strsplit(x,"")) )

sentence <- ""
is.null(sentence)
length(unlist(strsplit(sentence,""))) ==0 

a <- sapply(fsplit[[1]],function(x) as.numeric(unlist(strsplit(x,"")) [ grep("[0-9]+",unlist(strsplit(x,""))) ]) ) 
paste(names(sort(a)),collapse = " ")

lapply(unlist(fsplit),function(x) split(x,""))
  

grep("[0-9]+",unlist(strsplit("th1s",split="")))

as.character(strsplit("this",split=""))

totalSquares <- 14*20
squaresInRect(5,3)
squaresInRect(3,5)
squaresInRect(5,5)
squaresInRect(20,14)


lng <-3
wd <- 5
totalSquares <- lng*wd
maxit <- min(lng,wd)
vec <- c()
if(lng!=wd){
  while(sqrt(totalSquares)>1){
    item <- min(floor(sqrt(totalSquares)),min(lng,wd))
    vec <- append(vec,item)
    totalSquares <- totalSquares - item**2
    lng <- ifelse(lng-item==0,lng,lng-item)
    wd <- ifelse(wd-item==0,wd,wd-item)
    
    # maxit <- min(max(lng-item,wd-item),maxit)
  }
  final <- append(vec,rep(1,totalSquares))
}else{
  return(NULL)
}
final
sum(final**2)     
