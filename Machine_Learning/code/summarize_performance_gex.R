####Load librariers####
library(tidyverse)
library(ggplot2)
library(viridis)

####Aggregate individual performance files####
cv<-read.csv('~/Downloads/OneDrive_1_11-11-2020/MCNC_cv_scores_GEX_yscr.csv')%>%
  separate(splits, sep='\\.csv\\.', into=c('cell', 'modelinfo'))%>%
  separate(cell, into=c('cell_line', 'time', 'dose'), sep=', ')%>%
  separate(modelinfo, sep='\\.',
           into=c('dataset',
                  'method', 'true_or_y_scrambled_labels_used',
                  'train_test_split','cv_split'))%>%
  mutate(cell_line=str_extract(pattern='[\\d\\w]+', cell_line))%>%
  mutate('testset'='LOCO-CV', 'descriptor'='L1000')
ts<-read.csv('~/Downloads/OneDrive_1_11-11-2020/MCNCts_scores_GEX_yscr.csv')%>%
  separate(splits, sep='\\.csv\\.', into=c('cell', 'modelinfo'))%>%
  separate(cell, into=c('cell_line', 'time', 'dose'), sep=', ')%>%
  separate(modelinfo, sep='\\.',
           into=c('dataset',
                  'method', 'true_or_y_scrambled_labels_used',
                  'train_test_split','cv_split'))%>%
  mutate(cell_line=str_extract(pattern='[\\d\\w]+', cell_line))%>%
  mutate('testset'='External Test Set', 'descriptor'='L1000')

cv_mclcnc<-read.csv('~/Downloads/cv_Updated_MCLCNC_scores_GEX_yscr.csv')%>%
  separate(splits, sep='\\.csv\\.', into=c('cell', 'modelinfo'))%>%
  separate(cell, into=c('cell_line', 'time', 'dose'), sep=', ')%>%
  separate(modelinfo, sep='\\.',
           into=c('dataset',
                  'method', 'true_or_y_scrambled_labels_used',
                  'train_test_split','cv_split'))%>%
  mutate(cell_line=str_extract(pattern='[\\d\\w]+', cell_line))%>%
  mutate('testset'='LOCO-CV', 'descriptor'='L1000')
ts_mclcnc<-read.csv('~/Downloads/OneDrive_1_12-11-2020/MCLCNC_ts_scores_GEX_yscr.csv')%>%
  separate(splits, sep='\\.csv\\.', into=c('cell', 'modelinfo'))%>%
  separate(cell, into=c('cell_line', 'time', 'dose'), sep=', ')%>%
  separate(modelinfo, sep='\\.',
           into=c('dataset',
                  'method', 'true_or_y_scrambled_labels_used',
                  'train_test_split','cv_split'))%>%
  mutate(cell_line=str_extract(pattern='[\\d\\w]+', cell_line))%>%
  mutate('testset'='External Test Set', 'descriptor'='L1000')

df_performance<-rbind(cv,ts)%>%
  mutate(dose=gsub(dose,pattern='\\.0\\)', replacement = ''))%>%
  mutate(condition=ifelse(cell_line=='All_Cell_Lines','All cell lines',paste0(cell_line,' (',time,' h, ', dose,' uM)')))%>%
  mutate(condition=factor(condition, levels=unique(c('All cell lines', unique(condition)))))%>%
  dplyr::select(condition, balanced_accuracy,dataset,testset,true_or_y_scrambled_labels_used,method)%>%
  mutate(balanced_accuracy=as.numeric(balanced_accuracy))

df_performance$dataset_readable=NA
df_performance$dataset_readable[which(df_performance$dataset=='all')]<-'DILIrank'
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
  group_by(condition, dataset, testset, dataset_readable, method)%>%summarise("median_balacc"=median(balanced_accuracy, na.rm=T))%>%mutate('scrambled'="Scrambled performance")
df_standardmeans<-df_performance%>%
  filter(true_or_y_scrambled_labels_used ==1)%>%
  group_by(condition, dataset, testset, dataset_readable, method)%>%summarise("median_balacc"=mean(balanced_accuracy, na.rm=T))%>%mutate('mean'="Mean performance")
df_allmedians<-df_performance%>%
  filter(true_or_y_scrambled_labels_used ==1)%>%
  group_by(true_or_y_scrambled_labels_used, condition, dataset, testset, dataset_readable, method)%>%summarise("median_balacc"=median(balanced_accuracy, na.rm=T))

write.csv(df_allmedians, '../data/LOOCV_Results/median_performance_gex.csv')

gg_gex<- ggplot(data = df_models,aes(x=testset, y=as.numeric(balanced_accuracy ))) +
    geom_boxplot(data = df_models, aes(fill=condition))+
    geom_point(data=df_scrambledmedians,
               aes(y=median_balacc,fill=condition, shape=scrambled, group=condition),
               position=position_dodge(width=0.75), size=2)+
    geom_point(data=df_standardmeans,
               aes(y=median_balacc,fill=condition,  shape=mean, group=condition),
               position=position_dodge(width=0.75), size=2)+
    scale_shape_manual(values=c(23,24))+
    scale_fill_viridis_d(option='plasma', begin=0.3, end=1)+
    facet_wrap(~method) + theme_bw()+
    ylab("Balanced accuracy")+
  guides(guide_legend(ncol=2))+
    theme(text=element_text(size=13),
          #axis.text.x=element_text(angle=60, hjust=1),
          legend.title = element_blank(), axis.title.x = element_blank())



ggsave(gg_gex,filename = '../plots/gex_performance_balanced_accuracy.pdf', height = 130, width=250, units="mm")
ggsave(gg_gex,filename = '../plots/gex_performance_balanced_accuracy.jpeg',height = 130, width=250, units="mm")

gg_gex<- ggplot(data = df_models,aes(x=condition, y=as.numeric(balanced_accuracy ))) +
  geom_boxplot(data = df_models, aes(fill=testset))+
  geom_point(data=df_scrambledmedians,
             aes(y=median_balacc,fill=testset, shape=scrambled, group=testset),
             position=position_dodge(width=0.75), size=2)+
  geom_point(data=df_standardmeans,
             aes(y=median_balacc,fill=testset,  shape=mean, group=testset),
             position=position_dodge(width=0.75), size=2)+
  scale_shape_manual(values=c(23,24))+
 # scale_fill_viridis_d(option='plasma', begin=0.3, end=1)+
  facet_wrap(~method) + theme_bw()+
  ylab("Balanced Accuracy")+
  guides(guide_legend(ncol=2))+
  theme(text=element_text(size=13), axis.text.x=element_text(angle=60, hjust=1),legend.title = element_blank(), axis.title.x = element_blank())



ggsave(gg_gex,filename = '../plots/gex_performance.pdf', height = 250, width=300, units="mm")
ggsave(gg_gex,filename = '../plots/gex_performance.jpeg', height = 250, width=300, units="mm")

