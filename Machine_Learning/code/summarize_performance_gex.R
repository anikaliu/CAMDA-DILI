####Load librariers####
library(tidyverse)
library(ggplot2)
library(viridis)

####Aggregate individual performance files####
cv<-read.csv('~/Downloads/OneDrive_1_10-11-2020/MCNC_cv_scores_GEX_yscr.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset','rest','csv','data',
                  'method', 'true_or_y_scrambled_labels_used',
                  'train_test_split','cv_split'))%>%
  separate(dataset, into=c('cell_line', 'time', 'dose'), sep=', ')%>%
  mutate(cell_line=str_extract(pattern='[\\d\\w]+', cell_line))%>%
  mutate('testset'='LOCO-CV', 'descriptor'='L1000')
ts<-read.csv('~/Downloads/OneDrive_1_10-11-2020/MCNC_ts_scores_GEX_yscr.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset','rest','csv','data', 'method', 
                  'true_or_y_scrambled_labels_used',
                  'train_test_split','cv_split'))%>%
  separate(dataset, into=c('cell_line', 'time', 'dose'), sep=', ')%>%
  mutate(cell_line=str_extract(pattern='[\\d\\w]+', cell_line))%>%
  mutate('testset'='External Test Set', 'descriptor'='L1000')

df_performance<-rbind(cv,ts)%>%
  mutate(dataset='MCNC')%>%
  mutate(condition=ifelse(cell_line=='All_Cell_Lines','All cell lines',paste0(cell_line,' (',time,'h, ', dose,' uM)')))%>%
  mutate(condition=factor(condition, levels=c('All cell lines', unique(df_performance$condition))))%>%
  dplyr::select(condition, balanced_accuracy,dataset,testset,true_or_y_scrambled_labels_used)%>%
  mutate(balanced_accuracy=as.numeric(balanced_accuracy))

df_performance$dataset_readable=NA
df_performance$dataset_readable[which(df_performance$dataset=='all')]<-'DILIrank \n (+SIDER)'
df_performance$dataset_readable[which(df_performance$dataset=='MCLCNC')]<-'DILIrank'
df_performance$dataset_readable[which(df_performance$dataset=='MCNC')]<-'DILIrank \n (-vLessConcern)'
#write.csv(df_performance,'../data/Model_Results_Parameters/summarized_performance.csv', row.names = F)

####Generate figure####
#Convert to factors and assing orders
df_performance$dataset_readable<-factor(df_performance$dataset_readable, levels=c('DILIrank \n (-vLessConcern)','DILIrank','DILIrank \n (+SIDER)'))
df_performance$testset<-factor(df_performance$testset, levels=c('LOCO-CV','External Test Set','FDA Validation Set'))
#true_or_y_scrambled_labels_used =1 means standard model, everything else is scrambled

df_models<-df_performance%>%
  filter(true_or_y_scrambled_labels_used ==1)
df_scrambledmedians<-df_performance%>%
  filter(true_or_y_scrambled_labels_used !=1)%>%
  group_by(condition, dataset, testset, dataset_readable)%>%summarise("median_balacc"=median(balanced_accuracy))%>%mutate('scrambled'="Scrambled performance")
df_standardmeans<-df_performance%>%
  filter(true_or_y_scrambled_labels_used ==1)%>%
  group_by(condition, dataset, testset, dataset_readable)%>%summarise("median_balacc"=mean(balanced_accuracy, na.rm=T))%>%mutate('mean'="Mean performance")


gg_gex<- ggplot(data = df_models,aes(x=dataset, y=as.numeric(balanced_accuracy ))) +
    geom_boxplot(data = df_models, aes(fill=condition))+
    geom_point(data=df_scrambledmedians,
               aes(y=median_balacc,fill=condition, shape=scrambled, group=condition),
               position=position_dodge(width=0.75), size=2)+
    geom_point(data=df_standardmeans,
               aes(y=median_balacc,fill=condition,  shape=mean, group=condition),
               position=position_dodge(width=0.75), size=2)+
    scale_shape_manual(values=c(23,24))+
    scale_fill_viridis_d(option='plasma', begin=0.3, end=1)+
    facet_wrap(~testset) + theme_bw()+
    ylab("Balanced Accuracy")+
  guides(guide_legend(ncol=2))+
    theme(text=element_text(size=13), axis.text.x=element_text(angle=60, hjust=1),legend.title = element_blank(), axis.title.x = element_blank())



ggsave(gg_gex,filename = '../plots/gex_performance.pdf', height = 150, width=300, units="mm")
ggsave(gg_gex,filename = '../plots/gex_performance.jpeg', height = 150, width=300, units="mm")
