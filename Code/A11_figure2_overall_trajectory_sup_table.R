library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(survival)
library(fastcluster)
library(cluster)
library(survminer)
library(umap)
library(Rtsne)
library(enrichR)
library(DESeq2)
library(lme4)
library(lmerTest)
library(emmeans)
library(gplots)
library(stringr)
library(glmnet)
library(randomForest)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))]
RNAseq=RNAseq%>%mutate(Info_id = as.character(Info_id))

colnames(olink)[-(1:2)]=paste0('olink_',colnames(olink)[-(1:2)])
RNAseq = left_join(RNAseq,olink,by=c('Info_id'='id','Info_day'='day'))

colnames(ab1)[-(1:2)]=paste0('IgG_',colnames(ab1)[-(1:2)])
ab1 = ab1 %>%mutate(id = as.character(id))
RNAseq = left_join(RNAseq,ab1,by=c('Info_id'='id','Info_day'='day'))

colnames(lab1)[-(1:2)]=paste0('Lab_',colnames(lab1)[-(1:2)])
lab1 = lab1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
RNAseq = left_join(RNAseq,lab1,by=c('Info_id'='id','Info_day'='day'))

#colnames(eve1)[-(1:2)]=paste0('EVE_',colnames(eve1)[-(1:2)])
#eve1[,-(1:2)] = apply(eve1[,-(1:2)],2,function(x){log2(as.numeric(x))})
#eve1 = eve1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
#RNAseq = left_join(RNAseq,eve1,by=c('Info_id'='id','Info_day'='day'))

##### correlate with onset #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  filter(Info_severitycat!='Asymptomatic')

pLM = read.csv('Result/A11_trajectory.csv')
pLM = pLM%>%
  dplyr::select(gene,p,p2)%>%
  unique()%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  #filter(pAdj<0.05)%>%
  filter(nchar(gene)<40)



# remove random effect
dfTractory = data.frame()
days = data.frame(Info_onset=0:15)
V = c()
for (i in grep('^GO_|^IgG_|^olink_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c(15,i,17)]%>%na.omit()
  colnames(t1)[2]='y'
  if(var(t1$y)==0){next}
  LM1 =lmer(y ~ poly(Info_onset,2) + (1 | Info_id), data = t1)
  LM2 =lmer(y ~ 1 + (1 | Info_id), data = t1)
  p = anova(LM1, LM2)
  
  pValue = p$`Pr(>Chisq)`[2]
  t1 = data.frame(gene = colnames(df1)[i], 
                  p = pValue)
  dfTractory = rbind(dfTractory,t1)
  
}
pMix = dfTractory%>%
  dplyr::select(gene,p)%>%
  unique()%>%
  mutate(type = gsub('_.*','',gene))%>%
  group_by(type)%>%
  mutate(pAdj = p.adjust(p,'fdr'))%>%
  #filter(pAdj<0.05)%>%
  filter(nchar(gene)<40)

combined = inner_join(pLM,pMix,by='gene')%>%select(-type.x,-type.y)
colnames(combined)=c("gene","overall_p_value_linear_regression",
                 "quadratic_p_value_linear_regression",
                 "overall_FDR_linear_regression",
                 "overall_p_value_mixed_model",
                 "overall_FDR_mixed_model")
write.csv(combined,"Result/A11_sup_table.csv",row.names = F)
