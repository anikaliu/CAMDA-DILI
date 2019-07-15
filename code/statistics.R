runttest<-function(df_DILI, df_NoDILI){
  result<-as.data.frame(matrix(NA, nrow = ncol(df_DILI), ncol=12))
  colnames(result)<-c("statistic","parameter","p.value",'estimate1','estimate2',"conf.int1","conf.int2","null.value","alternative","method", "median1", "median2")
  rownames(result)<-colnames(df_DILI)
  for (index in 1:ncol(df_DILI)){
    feature=colnames(df_DILI)[index]
    contrast<-t.test(df_NoDILI[,feature], df_DILI[,feature], na.action='na.exclude', alternative = "less")
    result$statistic[index]=contrast$statistic
    result$parameter[index]=contrast$parameter
    result$p.value[index]=contrast$p.value
    result$null.value[index]=contrast$null.value
    result$conf.int1[index]=contrast$conf.int[1]
    result$conf.int2[index]=contrast$conf.int[2]
    result$estimate1[index]=contrast$estimate[1]
    result$estimate2[index]=contrast$estimate[2]
    result$median1<-median(df_NoDILI[,feature], na.rm = T)
    result$median2<-median(df_DILI[,feature], na.rm = T)
    result$alternative[index]=contrast$alternative
    result$method[index]=contrast$method
    result$null.value[index]=contrast$null.value
 
  }
  return(result)
}
  
  runttest_threshold<-function(df_DILI, df_NoDILI, threshold=2){
    result<-as.data.frame(matrix(NA, nrow = ncol(df_DILI), ncol=8))
    colnames(result)<-c("estimate","p.value",'freqDILI','freqNoDILI','DILI_inact','DILI_act','NoDILI_inact','NoDILI_act')
    rownames(result)<-colnames(df_DILI)
    for (index in 1:ncol(df_DILI)){
      feature=colnames(df_DILI)[index]
      NoDILI_act=sum(ifelse(df_NoDILI[,feature]>threshold,1,0)==1)
      NoDILI_inact=sum(ifelse(df_NoDILI[,feature]>threshold,1,0)==0)
      DILI_act=sum(ifelse(df_DILI[,feature]>threshold,1,0)==1)
      DILI_inact=sum(ifelse(df_DILI[,feature]>threshold,1,0)==0)
      TeaTasting <-
        matrix(c(DILI_act, NoDILI_act, DILI_inact,  NoDILI_inact),
      #         matrix(c(fp, tp, tn, fn),
               nrow = 2,
               dimnames = list(Guess = c('DILI', "NoDILI"),
                               Truth = c("1", "0")))
      fisher<-fisher.test(TeaTasting,, alternative = "greater")
      result$estimate[index]=fisher$estimate
      result$p.value[index]=fisher$p.value
      result$freqDILI[index]<-mean(ifelse(df_DILI[,feature]>threshold,1,0))
      result$freqNoDILI[index]<-mean(ifelse(df_NoDILI[,feature]>threshold,1,0))
      result$NoDILI_inact[index]<-NoDILI_inact
      result$NoDILI_act[index]<-NoDILI_act
      result$DILI_act[index]<-DILI_act
      result$DILI_inact[index]<-DILI_inact
      result$p.adj<-p.adjust(result$p.value,, method="fdr")
 
    }
    return(result)
  }
  runwilcoxon<-function(df_DILI, df_NoDILI, alternative="two.sided"){
    result<-as.data.frame(matrix(NA, nrow = ncol(df_DILI), ncol=2))
    colnames(result)<-c("statistic","p.value")
    rownames(result)<-colnames(df_DILI)
    for (index in 1:ncol(df_DILI)){
      feature=colnames(df_DILI)[index]
      contrast<-wilcox.test(x=df_DILI[,index],y=df_NoDILI[,index], alternative=alternative)
      result$statistic[index]=contrast$statistic
      result$p.value[index]=contrast$p.value
      
    }
    result$p.adj<-p.adjust(result$p.value, method="fdr")
    return(result)
  }
  
    
  #Example
  
  #df_DILI<-as.data.frame(t(read.csv('../processed_data/gex/gex-VDILI.csv', row.names = 'X')))
  #df_NoDILI<-as.data.frame(t(read.csv('../processed_data/gex/gex-vNO_VDILI.csv', row.names = 'X')))
  #result<-runwilcoxon(df_DILI, df_NoDILI)
  