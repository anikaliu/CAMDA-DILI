library("cmapR")
library("tidyverse")
library("dplyr")

"+" = function(x,y) { # useful function for string concatenation.
  if(is.character(x) || is.character(y)) {
    return(paste(x , y, sep=""))
  } else {
    .Primitive("+")(x,y)
  }
}

#load(file = "../data/l1000_1314pert_inames.872matchedBROAD_ids-dilirank.rda")
#broad_id_map = drank.sel
#write_csv("../data/l1000_1314pert_inames.872matchedBROAD_ids.csv", x = broad_id_map)


load(file = "../data/CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda")
load(file = "../data/CAMDA_l1000_1314compounds-dilirank.v2.rda")

Comp_GC = drank.sel

cdesc<-gct@cdesc
#write_csv("../data/GEXP_cdesc.csv", x = cdesc)
#cdesc = read.csv("../data/GEXP_cdesc.csv")


# (1) Check Cell line data for comps with other data
sum_cp<-cdesc%>%
  select(pert_iname,cell_id,pert_dose,pert_time)%>%
  unique()

data_GC = sum_cp[sum_cp$pert_iname %in% Comp_GC$Compound.Name,] # unique comp-cell-dose-time with GC data


lol_GC<-data_GC%>%
  select(pert_iname,cell_id,pert_dose,pert_time)%>%
  unique()%>%
  group_by(cell_id, pert_dose, pert_time)%>%
  summarise(n=n())%>%
  arrange(n)%>%
  ungroup()%>%
  filter(n==max(n))


# (2) Parse dataset for cell+dose+time conditions with data for all compounds

GE_data_GC = cdesc %>% inner_join(lol_GC, by=c("cell_id","pert_dose","pert_time")) # get full mat cell lines
GE_data_GC = GE_data_GC[GE_data_GC$pert_iname %in% Comp_GC$Compound.Name,] # get data only for comps with clinical data

metadata=inner_join(GE_data_GC, drank.sel)
write.csv(metadata,'../data/LINCS_DILI_metadata.csv')
#length(unique(GE_data_GC$pert_iname))

#GE_data_GC %>% count(pert_iname,cell_id,pert_dose,pert_time)

# (3) Weighted co-correlation of replicate signatures
# https://stackoverflow.com/questions/45532058/subsetting-a-matrix-on-the-basis-of-a-list-of-some-of-the-column-names
GE_data_GC_mat = as.data.frame(gct@mat) %>% select(one_of(GE_data_GC$id))


#Map gene ID to gene
gene_info <- read.csv("../data/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt", sep = '\t') ## PW moved
gene_info <- gene_info[gene_info$pr_is_lm == 1,] ## PW moved

# Get only landmark genes
GE_data_GC_landmark_non_collapsed <- GE_data_GC_mat[rownames(GE_data_GC_mat) %in% gene_info$pr_gene_id, ] # Landmark genes only)


genesymbols = gene_info[match(rownames(GE_data_GC_landmark_non_collapsed),gene_info$pr_gene_id ),2] # mapping to gene symbol
rownames(GE_data_GC_landmark_non_collapsed) = genesymbols
GE_data_GC_landmark_non_collapsed$id = genesymbols
write_csv("../data/L1000_Data_uncollapsed_all_conditionsv2_new.csv", x = GE_data_GC_landmark_non_collapsed)


# # Find replicate columns - loop over do weighted sum, save into new column and replace old replicates
# 
# modzs <- function(m) {
#   cc = cor(m, method="spearman") # pair-wise correlations
#   cc[cc<0] = 0
#   wt = 0.5*(colSums(cc)-1) # Per-sample values
#   wt[wt<0.01] = 0.01
#   sumWT = sum(wt)
#   normWT = wt / sumWT  # Normalized weights
#   # Return the scaled input values
#   return(m %*% normWT)
# }
# 
# conditions<-GE_data_GC%>%select(cell_id,pert_time, pert_dose, pert_iname)%>%unique()
# df<-data.frame(matrix(NA, nrow=nrow(GE_data_GC_mat), ncol=nrow(conditions)))
# 
# for (index in 1:nrow(conditions)){
#   print(index)
#  # print(nrow(conditions))
#   cond<-conditions[index,]
#   id_interest<-GE_data_GC%>%right_join(cond) # should be GE_data_GC
#   id_col<-id_interest$id
#   replicates<-GE_data_GC_mat%>%select(all_of(id_col)) # should be GE_data_GC_mat 
#   
#   #If there are multiple replicates aggregate using modzs
#   if (ncol(replicates)>1){
#     signature<-modzs(as.matrix(replicates))
#   }else{
#     signature<-replicates
#   }
#   colname<-paste0(cond[,1],'_', cond[,2],'_', cond[,3],'_', cond[,4])
#   df[,index]<-signature
#   colnames(df)[index]<-colname
# }
# 
# # Collapsed to only landmark genes
# dftest = df
# #rownames(dftest) = rownames(df$A375_6_10.0_alaproclate) ####AL column shouldn't have rownames
# rownames(dftest) = rownames(GE_data_GC_mat)
# dftest <- dftest[rownames(dftest) %in% gene_info$pr_gene_id, ] # Landmark genes only 
# genesymbols = gene_info[match(rownames(dftest),gene_info$pr_gene_id ),2] # mapping to gene symbol
# 
# rownames(dftest)= genesymbols
# dftest$id = genesymbols
# 
# 
# 
# write_csv("../data/L1000_Data_collapsed_all_conditionsv2_new.csv", x = dftest)
