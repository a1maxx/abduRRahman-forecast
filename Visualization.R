library(ggplot2)
library(ggthemes)
library(ggsci)



# result %>% get_solution(x[i,j,k]) %>% filter(value>0)
# result %>% get_solution(g[i,k]) %>% filter(value>0)
# result %>% get_solution(b[i,k])


 # write.csv(f.result.final2,"C:\\Users\\Administrator\\Desktop\\Datas\\newSetup4_10_48.csv", row.names = FALSE) 


f.result.final <-  read.csv("C:\\Users\\Administrator\\Desktop\\Datas\\newSetup3_5_100_24.csv") ### Used, created plots with



dist.names <- as_labeller(
  c(`euclidean` = "Euclidean Distance", `manhattan` = "Manhattan Distance"))




f.result.final2 <- f.result.final %>% mutate(optimization = as.numeric(optimization),simulation = as.numeric(simulation),iter = as.numeric(iter)) %>% 
  mutate(p.diff =  ((simulation-optimization) / optimization) * 100, lat = factor(lat,levels= c(as.character(unique(f.result.final$lat)))))


#### Aggregate plots

zz <- f.result.final2 %>% group_by(iter) %>% summarise(Meaniter = mean(simulation),SDiter = sd(simulation))
f.result.final3 <- f.result.final2 
f.result.final3$meanIter <- rep(zz$Meaniter,each=length(c(as.character(unique(f.result.final$lat)))),2) 
f.result.final3$sdIter <- rep(zz$SDiter,each=length(c(as.character(unique(f.result.final$lat)))),2)
f.result.final3 <-  f.result.final3 %>% mutate(normPerf = (simulation-meanIter)/sdIter)

font_add_google("Montserrat", "Montserrat")
font_add_google("Roboto", "Roboto")
font_add("Palatino", "palab.ttf")
showtext_auto()
showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)
myFont1 <- "Montserrat"
myFont2 <- "Roboto"
myFont3 <- "Palatino"

figure1 <- ggplot(f.result.final3,aes(x=lat,y=normPerf,color=lat)) + geom_boxplot() + ylim(c(-2.5,2.5)) +
  facet_wrap(~distmethod,scales = "free", labeller = dist.names) + theme_classic() + scale_color_aaas(name = "Latency",labels= c("High","Moderate","Low")) + 
  ylab("Normalized Cost") + theme(strip.text = element_text( face= "bold",size = unit(14,"pt"),family = myFont3),
                                  legend.text = element_text(face="bold",family=myFont3),
                                  legend.title = element_text(face="bold",family=myFont3),
                                  axis.title = element_text(face="bold",family= myFont3,size = unit(14,"pt")),
                                  axis.text = element_text(face="bold",family=myFont3,size = unit(12,"pt")),
                                  legend.position = c(0.9,0.85)) +scale_x_discrete(labels=c("High","Moderate","Low")) + xlab("Latency")


figure1

# 
# ggsave(figure1, filename = "C:\\Users\\Administrator\\Dropbox\\BookChapter-vol3\\Figures\\figure1_cairo.png", dpi = 300, type ="cairo",
#        width =8, height = 6, units = "in")
dev.off()

showtext_opts(dpi = 300)
figure2 <- ggplot(f.result.final3,
       aes(
         x = rep(rep(1:24, 3), 2),
         y = normPerf,
         color = lat,
         linetype = lat
       )) + theme_minimal() + geom_smooth(method="gam",formula = y~poly(x,6), se = FALSE) + theme_classic() + 
  # facet_wrap(~distmethod,labeller = dist.names) +od="gam",formula = y~poly(x,5), se = FALSE) + theme_classic() + facet_wrap(~distmethod,labeller = dist.names)+
  theme(strip.text = element_text( face= "bold",size = unit(14,"pt"),family = myFont3),
        legend.text = element_text(face="bold",family=myFont3),
        legend.title = element_text(face="bold",family=myFont3),
        axis.title = element_text(face="bold",family= myFont3,size = unit(14,"pt")),
        axis.text = element_text(face="bold",family=myFont3,size = unit(12,"pt")),
        legend.position = "top",
        legend.key = element_rect(colour ="transparent", fill = "white")) + xlab("Time Steps") +
  labs(color = NULL, linetype = NULL, shape = NULL) + 
  xlab("Time Steps") + ylab("Normalized Objective Value") + 
  scale_linetype_discrete(name="Latency", labels=c("High","Moderate","Low")) + scale_color_aaas(name="Latency", labels=c("High","Moderate","Low")) 

figure2
ggsave(figure2, filename = "C:\\Users\\Administrator\\Dropbox\\BookChapter-vol3\\Figures\\figure2_cairo.png", dpi = 300, type ="cairo",
        width =8, height = 6, units = "in")


#### Plots for manhattan distance


aa <- f.result.final2 %>% filter(distmethod == "manhattan") %>% group_by(iter) %>% summarise(Meaniter = mean(simulation),SDiter = sd(simulation))

f.result.final3 <- f.result.final2 %>% filter(distmethod == "manhattan")
f.result.final3$meanIter <- rep(aa$Meaniter,each=length(c(as.character(unique(f.result.final$lat))))) 
f.result.final3$sdIter <- rep(aa$SDiter,each=length(c(as.character(unique(f.result.final$lat)))))
f.result.final3 <-  f.result.final3 %>% mutate(normPerf = (simulation-meanIter)/sdIter)
  
ggplot(f.result.final3,aes(x=iter,y=normPerf,fill=lat)) + geom_boxplot() + ylim(c(-3,3))

figure3 <- ggplot(f.result.final3,aes(x=lat,y=p.diff,fill = lat)) + geom_boxplot(outlier.size = -1,alpha= 0.65) +
  scale_fill_aaas(name="Latency") + ylim(c(min(f.result.final3$p.diff),80)) + ylab("( % )Deviation") + 
  guides(fill="none") + scale_x_discrete(name="Latency", labels = c("High","Moderate","Low"))  +
  theme(axis.title.y = element_text(face = "bold",margin = margin(t = 1,r = 5, b = 1,l =1)),
        strip.text = element_text( face= "bold",size = unit(14,"pt"),family = myFont3),
       legend.text = element_text(face="bold",family=myFont3),
       legend.title = element_text(face="bold",family=myFont3),
       axis.title = element_text(face="bold",family= myFont3,size = unit(14,"pt")),
       axis.text = element_text(face="bold",family=myFont3,size = unit(12,"pt")),
       legend.position = c(0.9,0.85)) + xlab("Latency")




ggsave(figure3, filename = "C:\\Users\\Administrator\\Dropbox\\BookChapter-vol3\\Figures\\figure3_cairo.png", dpi = 300, type ="cairo",
       width =4, height = 3, units = "in")



#### Plots for euclidean distance


bb <- f.result.final2 %>% filter(distmethod == "euclidean") %>% group_by(iter) %>% summarise(Meaniter = mean(simulation),SDiter = sd(simulation))

f.result.final3 <- f.result.final2 %>% filter(distmethod == "euclidean")
f.result.final3$meanIter <- rep(aa$Meaniter,each=length(c(as.character(unique(f.result.final$lat))))) 
f.result.final3$sdIter <- rep(aa$SDiter,each=length(c(as.character(unique(f.result.final$lat)))))
f.result.final3 <-  f.result.final3 %>% mutate(normPerf = (simulation-meanIter)/sdIter)

ggplot(f.result.final3,aes(x=iter,y=normPerf,fill=lat)) + geom_boxplot() + ylim(c(-5,5))

ggplot(f.result.final3,aes(x=lat,y=p.diff,fill = lat)) + geom_boxplot(outlier.size = -1,alpha= 0.65) +
  scale_fill_aaas(name="Latency") + ylim(c(min(f.result.final3$p.diff),80)) + ylab("( % )Deviation") + 
  theme(axis.title =  element_text(face = "bold",size = unit(14,"pt")),
        axis.title.y = element_text(face = "bold",margin = margin(t = 1,r = 5, b = 1,l =1),family = 'sans'  ),
        axis.text = element_text(face = "bold")) + guides(fill="none") + scale_x_discrete(name="Latency", labels = c("High","Moderate","Low"))



grouped.df <- f.result.final2 %>% filter(distmethod == "manhattan") %>%  group_by(iter,lat) %>% 
  summarise(meanSim = mean(simulation),meanDif = mean(p.diff),meanOpt= mean(optimization)) %>% as.data.frame()


ggplot(grouped.df,aes(x=iter,y=meanSim,fill=lat)) + geom_boxplot()   
  

grouped.df <- f.result.final2 %>%  group_by(iter,lat,distmethod) %>% summarise(meanSim = mean(simulation),meanDif = mean(p.diff),meanOpt= mean(optimization)) %>% 
  as.data.frame()
  
  ggplot(grouped.df,aes(x=lat,y=meanSim,fill=distmethod)) + geom_boxplot() 
  
  ggplot(grouped.df,aes(x=lat,y=meanDif,fill=lat)) + geom_boxplot() 
  
  ggplot(grouped.df,aes(x=iter,y=meanOpt,fill=lat)) + geom_boxplot() 

  
  ggplot(grouped.df,aes(x=iter,y=meanSim,color=lat,linetype=lat)) + geom_line(lwd=1) + scale_color_futurama() +ylab("Objective Function Value")
  

  ggplot(f.result.final2,aes(x=lat,y=simulation)) + geom_boxplot()  + theme_clean()
  
  
  v1 <- f.result.final3[f.result.final3$lat==100,"normPerf"]
  v2 <- f.result.final3[f.result.final3$lat==3,"normPerf"]
  
  t.test(v1,v2,conf.level = 0.95)


p <-  ggplot(f.result.final2 %>% filter(distmethod == "euclidean"),aes(x=rep(rep(1:24,3),1),y=simulation,color=lat,linetype=lat)) +theme_minimal()

p + geom_line()




plt2 <- p + geom_smooth(method="lm", formula = y~poly(x,6), se = FALSE) + labs(y = "")+ theme_base() +
  theme(strip.text = element_text(face = "bold", size = unit(14, "pt")), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        legend.position = c(.7, .7), 
        axis.title = element_text(size = unit(16, "pt"), face = "bold"), 
        legend.text = element_text(face = "bold"))  + 
        labs(color = NULL, linetype = NULL, shape = NULL) + 
        xlab("Time Steps") + ylab("Objective Value") +
  scale_linetype_discrete(name="Latency", labels=c("High","Moderate","Low")) + scale_color_aaas(name="Latency", labels=c("High","Moderate","Low")) 

plt2



p + geom_smooth(method="auto",formula = y~x, se = FALSE) + labs(y = "")+ theme_base() +
  theme(strip.text = element_text(face = "bold", size = unit(14, "pt")), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        legend.position = "top", 
        axis.title = element_text(size = unit(16, "pt"), face = "bold"), 
        legend.text = element_text(face = "bold"))  + 
  labs(color = NULL, linetype = NULL, shape = NULL) + 
  xlab("Time Steps") + ylab("Objective Value") +
  scale_linetype_discrete(name="Latency", labels=c("High","Moderate","Low")) + scale_color_aaas(name="Latency", labels=c("High","Moderate","Low")) 




# ggsave(plt2, filename = paste(as.character(Sys.time()),".png"), dpi = 800, type = "cairo",
#        width =8, height = 6, units = "in")


a <- ggplot(f.result.final2,aes(x=lat,y=p.diff*100,fill=distmethod)) + geom_boxplot(outlier.size = -1) +
    xlab("Latency") + ylab("( % )Deviation") + 
  theme(axis.title = element_text(face="bold",size = unit(16,"pt")),
        axis.text = element_text(face = "bold",size = unit(14,"pt"))) +
  scale_x_discrete(labels=c("High","Moderate","Low")) +scale_fill_aaas()
  
a

# 
# ggsave(a, filename = "diff.png", dpi = 800, type = "cairo",
#        width =8, height = 6, units = "in")
  


library(Rtsne)
scenarios <- generateScenarios(demandL,rgenL)

reducted.scenarios0 <- scenarios %>% filter(prob>fivenum(prob)[3]) %>% arrange(desc(prob))

reducted.scenarios <- head(reducted.scenarios0,500) %>% select(-prob)

scaled.reducted <- reducted.scenarios %>% scale()

dist.scenarios <- daisy(scaled.reducted,meth)

mat.dist.scenarios <- as.matrix(dist.scenarios)


sil_width <- c(-Inf)
for(i in 2:13){
  
  pam_fit <- pam(mat.dist.scenarios,
                 diss = TRUE,
                 k = i,nstart = 5)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}


plot(sil_width,type = "l")
noc <- which(sil_width == max(sil_width))
pam_res <- pam(scaled.reducted,noc)


tsne_fit <- Rtsne(scaled.reducted,pca_center = F,dims=2) 
tsne_df <- tsne_fit$Y %>% 
  as.data.frame() %>% cbind(cluster = as.factor(pam_res$clustering)) %>% 
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())


pam.cent = tsne_df[pam_res$id.med,] %>% select(tSNE1, tSNE2,cluster) 
km.cent = tsne_df %>% group_by(cluster) %>% select(tSNE1, 
                                                    tSNE2) %>% summarize_all(mean)

figure6 <- tsne_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             colour = cluster))+
  geom_point(alpha = 0.5,size=4 ) + geom_label_repel(aes(label=cluster),data = km.cent) +
  theme(legend.position="bottom") + guides(color="none")


ggsave(figure6, filename = "C:\\Users\\Administrator\\Dropbox\\BookChapter-vol3\\Figures\\figure6_cairo.png", dpi = 300, type ="cairo",
       width =4, height = 3, units = "in")




