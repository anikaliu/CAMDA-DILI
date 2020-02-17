#https://stackoverflow.com/questions/47542849/marginal-plots-using-axis-canvas-in-cowplot-how-to-insert-gap-between-main-pane

# http://www.lreding.com/nonstandard_deviations/2017/08/19/cowmarg/

library(tidyverse) # for ggplot2

library(magrittr) # for pipes and %<>%

library(ggpubr) # for theme_pubr()

library(cowplot)

library(ggsci)


df <- read_csv("./Structural_alerts_with_DrugBank.csv")
df <- df[df$p_value <= 0.05, ]


df$Source<-as.factor(df$Source)

levels(df$Source)<-c("Liu et al. (2015) [n=5]", "SARpy [n=11]", "MOSS [n=23]")

#levels(df$Source)<-c("Hewitt et al. (2013) [n=16]", "Liu et al. (2015) [n=5]", "SARpy [n=11]", "MOSS [n=21]")
#df<-df[!(df$Source=="Hewitt et al. (2013) [n=16]"),]


original_plot <- ggplot(data = df, mapping = aes(x = PPV, y = perc_hits, color=Source)) +
  
  geom_point(alpha = 0.6, size=1.6) +
  
  xlab("Precision") +
  
  ylab("Compounds with positive DILI labels \n containing the respective substructure (%)")  +
  
  scale_color_futurama()+
  
  xlim(0, 1) +
  
  ylim(0, 100) +
  
  #  scale_color_brewer(palette='Set2')+
  
  theme(legend.position="right") + theme_bw()+ labs(color = "Substructure Source")


original_plot <- original_plot + theme(axis.text = element_text(size = 13))
original_plot <- original_plot + theme(axis.title = element_text(size = 13))
original_plot <- original_plot + theme(legend.text = element_text(size = 13))
original_plot <- original_plot + theme(legend.title = element_text(size = 13))
original_plot

y_box <- axis_canvas(original_plot, axis = "y") +
  
  geom_boxplot(data = df, aes(x=factor(Source),y = perc_hits, color = factor(Source))) +scale_color_futurama()+ scale_x_discrete() 


x_box <- axis_canvas(original_plot,axis="x", coord_flip = TRUE) + 
  
  geom_boxplot(data=df, aes(y=PPV, x = factor(Source), color = factor(Source))) + scale_color_futurama()+ scale_x_discrete() +coord_flip() 



#l <- ggplot(data = df, aes(x=factor(Source),y = perc_hits, color = factor(Source))) +

#  geom_boxplot()

#ggdraw(l)

p1 <- insert_xaxis_grob(original_plot,x_box,grid::unit(0.7, "in"),position = "top")

p2 <- insert_yaxis_grob(p1, y_box, grid::unit(0.7, "in"),position = "right")

ggdraw(p2)


ggsave(p2, filename = './Plots/Figure5.jpeg', height = 95, width=170, units="mm")
ggsave(p2, filename = './Plots/Figure5.pdf', height = 95, width=170, units="mm")

#p3 <- p2 + 

#    scale_fill_discrete(name = "Substructure Source", labels = c("MOSS", "SARpy", "Hewitt", "Liu"))

