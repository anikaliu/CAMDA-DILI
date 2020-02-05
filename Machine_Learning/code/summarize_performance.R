####Load librariers####
library(tidyverse)
library(ggplot2)

####Aggregate individual performance files####
cv_ECFP<-read.csv('../data/Model_Results_Parameters/ECFP/cv_scores_ECFP.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='CV', 'descriptor'='ECFP')
ts_ECFP<-read.csv('../data/Model_Results_Parameters/ECFP/test_scores_ECFP.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='External Test Set', 'descriptor'='ECFP')
ex_ECFP<-read.csv('../data/FDA_Validation_Set_Results/Predictions_ECFP4/ambiguous_results.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='FDA Validation Set', 'descriptor'='ECFP')
cv_MD<-read.csv('../data/Model_Results_Parameters/MD/cv_scores_MD.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='CV', 'descriptor'='MD')
ts_MD<-read.csv('../data/Model_Results_Parameters/MD/test_scores_MD.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='External Test Set', 'descriptor'='MD')
ex_MD<-read.csv('../data/FDA_Validation_Set_Results/Predictions_MD/ambiguous_results.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='FDA Validation Set', 'descriptor'='MD')
cv_PT<-read.csv('../data/Model_Results_Parameters/PT/cv_scores_PT.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='CV', 'descriptor'='PT')
ts_PT<-read.csv('../data/Model_Results_Parameters/PT/test_scores_PT.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='External Test Set', 'descriptor'='PT')
ex_PT<-read.csv('../data/FDA_Validation_Set_Results/Predictions_PT/ambiguous_results.csv')%>%
  separate(splits, sep='\\.',
           into=c('dataset', 'method', 'true_or_y_scrambled_labels_used','train_test_split','cv_split'))%>%
  mutate('testset'='FDA Validation Set', 'descriptor'='PT')
df_performance<-cv_ECFP%>%
  bind_rows(ts_ECFP)%>%
  bind_rows(cv_MD)%>%
  bind_rows(ts_MD)%>%
  bind_rows(cv_PT)%>%
  bind_rows(ts_PT)%>%
  bind_rows(ex_ECFP)%>%
  bind_rows(ex_MD)%>%
  bind_rows(ex_PT)
df_performance<- df_performance%>%select(dataset, method, descriptor, testset, everything())
df_performance$dataset_readable=NA
df_performance$dataset_readable[which(df_performance$dataset=='all')]<-'DILIrank \n (+SIDER)'
df_performance$dataset_readable[which(df_performance$dataset=='MCLCNC')]<-'DILIrank'
df_performance$dataset_readable[which(df_performance$dataset=='MCNC')]<-'DILIrank \n (-vLessConcern)'
write.csv(df_performance,'../data/Model_Results_Parameters/summarized_performance.csv', row.names = F)

####Generate figure####
#Convert to factors and assing orders
df_performance$dataset_readable<-factor(df_performance$dataset_readable, levels=c('DILIrank \n (-vLessConcern)','DILIrank','DILIrank \n (+SIDER)'))
df_performance$testset<-factor(df_performance$testset, levels=c('CV','External Test Set','FDA Validation Set'))
#true_or_y_scrambled_labels_used =1 means standard model, everything else is scrambled
df_models<-df_performance%>%
  filter(true_or_y_scrambled_labels_used ==1)
df_scrambledmedians<-df_performance%>%
  filter(true_or_y_scrambled_labels_used !=1)%>%
  group_by(dataset_readable, method, testset)%>%summarise("median_balacc"=median(balanced_accuracy))%>%mutate('scrambled'="Scrambled performance")
df_standardmeans<-df_performance%>%
  filter(true_or_y_scrambled_labels_used ==1)%>%
  group_by(dataset_readable, method, testset)%>%summarise("median_balacc"=mean(balanced_accuracy))%>%mutate('mean'="Mean performance")

#Plotting function
plot_performance<-function(df_models){
  g <- ggplot(data = df_models,aes(x=dataset_readable, y=balanced_accuracy )) +
    geom_boxplot(data = df_models, aes(fill=method))+
    geom_point(data=df_scrambledmedians,
               aes(y=median_balacc,fill=method, shape=scrambled, group=method),
               position=position_dodge(width=0.75), size=2)+
    geom_point(data=df_standardmeans,
               aes(y=median_balacc,fill=method,  shape=mean, group=method),
               position=position_dodge(width=0.75), size=2)+
    scale_shape_manual(values=c(23,24))+
    scale_fill_brewer(palette='Set2')+
    facet_wrap(~testset) + theme_bw()+
    ylab("Balanced Accuracy")+
    theme(text=element_text(size=15), axis.text.x=element_text(angle=60, hjust=1),legend.title = element_blank(), axis.title.x = element_blank())
  return(g)
}

#Generate and save
gg_ECFP<-plot_performance(df_models%>%filter(descriptor=='ECFP'))
ggsave(gg_ECFP,filename = '../plots/ECFP_performance.pdf', height = 4.5, width=8)

gg_MD<-plot_performance(df_models%>%filter(descriptor=='MD'))
ggsave(gg_MD,filename = '../plots/MD_performance.pdf', height = 4.5, width=8)

gg_PT<-plot_performance(df_models%>%filter(descriptor=='PT'))
ggsave(gg_PT,filename = '../plots/PT_performance.pdf', height = 4.5, width=8)
