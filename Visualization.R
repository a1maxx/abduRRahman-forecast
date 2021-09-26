library(ggplot2)
library(ggthemes)
library(ggsci)



f.result.final2 <- f.result.final[-1,] %>% 
  mutate(p.dif =  (simulation-optimization) / optimization, lat = as.factor(lat) )


ggplot(f.result.final2,aes(x=lat,y=simulation)) + geom_boxplot()  

ggplot(f.result.final2,aes(x=lat,y=p.dif)) + geom_boxplot() 

f.result.final3  <- f.result.final2 %>% group_by(lat) %>%  summarise(fivenum(p.dif))


f.result4 <- f.result.final2

f.result4.melted <- f.result4 %>% select(-p.dif,-simulation) %>% melt()

ggplot(f.result4.melted,aes(x=lat,y=value,color=variable)) + geom_boxplot()

f.result.final2 %>% select(lat,p.dif) %>% ggplot(aes(x=lat,y=p.dif)) + 
  geom_boxplot()



ggplot(f.result4) + geom_point(aes(x, y = value, shape = variable, color =
                                     variable))  +
  geom_smooth(aes(group = variable))+ 
  theme_bw() + facet_wrap(~ lat) +
  theme(
    axis.text.x = element_text(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  scale_color_manual(name = 'legend',
                     labels = c("optimization","simulation","moderate","none"),
                     values = c("firebrick1", "deepskyblue1","green","red")) +
  scale_shape_manual(name ='legend',
                     labels = c("optimization","simulation","moderate","none"),
                     values = c(1,2,3,4)) 
  


theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

p <- ggplot(f.result4,aes(x=x,y=value,shape=variable, color=latency)) +facet_wrap(~latency)  

pg <- p + geom_point()

pg + stat_smooth(
  method = "gam",
  formula = y ~
    s(x),
  size = 1,
  se = FALSE) + 
  theme_clean() +scale_color_aaas() + labs(y="")+
  theme(strip.text = element_text(face="bold",size = unit(14,"pt")),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.9,.8),
        axis.title=element_text(size=unit(12,"pt"),face="bold"),
        legend.text = element_text(face="bold")) 


pg + stat_smooth(
  method = "loess",
  formula = y ~
    x,
  size = 1,
  se = FALSE) + 
  theme_clean() +scale_color_aaas() + labs(y="",colour=NULL)+
  theme(strip.text = element_text(face="bold",size = unit(14,"pt")),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.5,.5),
        axis.title=element_text(size=unit(12,"pt"),face="bold"),
        legend.text = element_text(face="bold",size=unit(8,"pt"))) +  
  labs(color = NULL, linetype = NULL, shape = NULL)

themeF <- function(p){
  print(p + theme_clean() +scale_color_aaas() + labs(y = "")+
          theme(strip.text = element_text(face = "bold", size = unit(14, "pt")), 
                strip.background = element_rect(colour = "white", fill = "white"), 
                legend.position = c(.9, .8), 
                axis.title = element_text(size = unit(12, "pt"), face = "bold"), legend.text = element_text(face =
                                                                                                              "bold")) )
}

pg + stat_smooth(
  method = "gam",
  formula = y ~
    s(x),
  size = 1,
  se = FALSE) + 
  theme_clean() +scale_color_aaas() + labs(y = "")+
  theme(strip.text = element_text(face = "bold", size = unit(14, "pt")), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        legend.position = c(.9, .8), 
        axis.title = element_text(size = unit(12, "pt"), face = "bold"), legend.text = element_text(face = "bold")) +labs(shape="a",color="a")

pg + geom_smooth(method="lm", formula = y~poly(x,10), se = FALSE) +  theme_clean() +scale_color_aaas() + labs(y = "")+
  theme(strip.text = element_text(face = "bold", size = unit(14, "pt")), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        legend.position = c(.9, .8), 
        axis.title = element_text(size = unit(12, "pt"), face = "bold"), legend.text = element_text(face =
                                                                                                      "bold"))  + 
  labs(color = NULL, linetype = NULL, shape = NULL)

