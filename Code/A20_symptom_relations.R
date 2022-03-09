##### load packages #####
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
library(xCell)
library(glmnet)
library(randomForest)
library(PRROC)
library(pROC)
##### load data #####
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
RNAseq = RNAseq[,grepl('olink|Info',colnames(RNAseq))]
RNAseq=RNAseq%>%mutate(Info_id = as.character(Info_id))

colnames(olink)[-(1:2)]=paste0('olink_',colnames(olink)[-(1:2)])
RNAseq = inner_join(RNAseq,olink,by=c('Info_id'='id','Info_day'='day'))

##### join with clinical data #####
neut_ab = read.csv('Data/Antibodies/D28_M7 Neut IgG Binding.csv')
df_combine = neut_ab%>%
  dplyr::select(Info_id=ID,D28_IgG=d28_IgG_AUC,D28_Neut,
                M7_IgG = M7_IgG_AUC,M7_Neut)%>%
  mutate(D28_IgG = log2(D28_IgG),
         D28_Neut=log2(D28_Neut),
         M7_IgG=log2(M7_IgG) ,
         M7_Neut=log2(M7_Neut))%>%
  mutate(Info_id = as.character(Info_id))
RNAseq = right_join(df_combine,RNAseq,by='Info_id')

tcell = read.csv('Data/Tcell_Data/Lambda ICS D28_combination cytokines_bgsub_freq of parent.csv')
RNAseq = tcell%>%
  filter(Stim%in%c('MN','SS1'))%>%
  dplyr::select(Info_id=id,mCD4)%>%
  group_by(Info_id)%>%summarise(mCD4=sum(mCD4))%>%
  right_join(RNAseq,by='Info_id')

RNAseq = opAUC%>%mutate(id = as.character(id))%>%
  right_join(RNAseq,by=c('id'='Info_id'))

colnames(RNAseq) = gsub('\\-| ','',colnames(RNAseq))
txt1='olink|mCD4|^auc|M7_Neut|D28_Neut|Info_age|Info_sex|Info_severitycat'
RNAseq = RNAseq[,grepl(txt1,colnames(RNAseq))]

##### plot #####
RNAseq = RNAseq[,!grepl('olink|age|sex',colnames(RNAseq))]
RNAseq = RNAseq%>%
  gather(key='clin',value = 'value',-Info_severitycat)%>%
  na.omit()%>%filter(is.finite(value))
p = ggplot(RNAseq,aes(x =Info_severitycat,y=value,
                      color = Info_severitycat))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(width = 0.2)+
  facet_wrap(~clin,nrow = 1,scales = 'free')+
  theme_bw()
pdf('Result/A20_severity_vs_all.pdf',width = 12, height = 3)
plot(p)
dev.off()

f1 = function(value,Info_severitycat){
  L1 = lm(value~Info_severitycat)
  L1 = summary.aov(L1)
  return(L1[[1]]['Info_severitycat','Pr(>F)'])
}
pvalue = RNAseq%>%
  group_by(clin)%>%
  summarise(p = f1(value,Info_severitycat))





