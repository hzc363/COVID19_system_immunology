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
#library(xCell)
library(glmnet)
library(randomForest)
library(PRROC)

##### load data #####
# ab1,demo1,eve1,symp1,olink,RNAseq,timeTo1,tcell,go_df
load('Result/A01_compiled_data.rda')
RNAseq = fread('Result/A02_RNA_seq_derived.csv',data.table = F)
#t1 = c('Gene_CCR2','Gene_CSF1R','Gene_LIFR',
#       'Gene_OSMR','Gene_IL6ST','Gene_PDCD1','Gene_CXCR3')
geneDf = read.csv('A27_top_go_genes.csv')
t1 = paste0('Gene_',geneDf$gene_symbol)
RNAseq = RNAseq[,!grepl('^Gene',colnames(RNAseq))|colnames(RNAseq)%in%t1]

RNAseq=RNAseq%>%mutate(Info_id = as.character(Info_id))

colnames(olink)[-(1:2)]=paste0('olink_',colnames(olink)[-(1:2)])
RNAseq = left_join(RNAseq,olink,by=c('Info_id'='id','Info_day'='day'))

colnames(ab1)[-(1:2)]=paste0('IgG_',colnames(ab1)[-(1:2)])
ab1 = ab1 %>%mutate(id = as.character(id))
RNAseq = left_join(RNAseq,ab1,by=c('Info_id'='id','Info_day'='day'))

colnames(lab1)[-(1:2)]=paste0('Lab_',colnames(lab1)[-(1:2)])
lab1 = lab1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
RNAseq = left_join(RNAseq,lab1,by=c('Info_id'='id','Info_day'='day'))

colnames(eve1)[-(1:2)]=paste0('EVE_',colnames(eve1)[-(1:2)])
eve1[,-(1:2)] = apply(eve1[,-(1:2)],2,function(x){log2(as.numeric(x))})
eve1 = eve1 %>%mutate(id = as.character(id))%>%mutate(day = as.numeric(day))
RNAseq = left_join(RNAseq,eve1,by=c('Info_id'='id','Info_day'='day'))

##### correlate with severity #####
df1 = RNAseq%>%
  filter(Info_tx=='placebo')%>%
  mutate(Info_severitycat = factor(Info_severitycat,
                                   levels =c('Asymptomatic',
                                             "Moderate Symptomatic",
                                             "Severe")))%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(!is.na(Info_onset))%>%
  mutate(Info_severitycat = as.integer(Info_severitycat))
  

dfResult = data.frame()
for (i in grep('^GO_|^IgG_|^Lab_|^olink_|^xcell_|Gene_|EVE_',colnames(df1))) { #18:ncol(df1)
  t1 = df1[,c(15,16,i)]%>%na.omit()
  colnames(t1)[3]='y'
  LM1 = lm(y~poly(Info_onset,2)+Info_severitycat,t1)
  LM1 = summary(LM1)
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = LM1$coefficients['Info_severitycat','Estimate'], 
                  t= LM1$coefficients['Info_severitycat','t value'],
                  p = LM1$coefficients['Info_severitycat','Pr(>|t|)'])
  dfResult = rbind(dfResult,t1)
}

write.csv(dfResult,'Result/A12_severity_results_control_only.csv',row.names = F)

dfResultAll = read.csv('Result/A12_severity_results.csv')
dfResultAll = dfResultAll%>%inner_join(dfResult,by='gene')%>%
  dplyr::select(gene,tAll = t.x,tControl=t.y)%>%
  mutate(type = gsub('_.*','',gene))%>%
  filter(type %in%c('GO','olink'))

p = ggplot(dfResultAll,aes(x=tAll,y = tControl))+
  geom_point()+theme_bw()+
  facet_wrap(~type,scales = 'free')+
  geom_smooth(method = 'lm')+
  xlab('t stat using full dataset')+
  ylab('t stat using control arm')
pdf('Result/A12_severity_correlation_arms.pdf',width = 8,height = 5)  
print(p)
dev.off()

dfResultAll = dfResultAll%>%group_by(type)%>%
  summarise(cor = cor(tAll,tControl,method = 'spearman'))










