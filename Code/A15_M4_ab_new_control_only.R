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
library(corrplot)

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

##### combined with Ab data ######
neut_ab = read.csv('Data/Antibodies/D28_M7 Neut IgG Binding.csv')
df_combine = neut_ab%>%
  dplyr::select(Info_id=ID,D28_IgG=d28_IgG_AUC,D28_Neut,
                M7_IgG = M7_IgG_AUC,M7_Neut)%>%
  mutate(D28_IgG = log2(D28_IgG+1),
         D28_Neut=log2(D28_Neut+1),
         M7_IgG=log2(M7_IgG+1) ,
         M7_Neut=log2(M7_Neut+1))%>%
  mutate(Info_id = as.character(Info_id))
RNAseq = inner_join(df_combine,RNAseq,by='Info_id')

##### correlate with Ab #####
df1 = RNAseq%>%
  mutate(Info_onset = Info_day-Info_onset)%>%
  filter(Info_severitycat!='Asymptomatic')%>%
  filter(Info_tx=='placebo')

plotDf = df1%>%dplyr::select(Info_id,Info_tx,D28_IgG)%>%unique()
p = ggplot(plotDf,aes(x=Info_tx,y =D28_IgG ))+
  geom_boxplot()+theme_bw()
pdf('Result/A14_tx_vs_D28_IgG.pdf',3,3)
plot(p)
dev.off()

dfResult = data.frame()
for (i in grep('^GO_|^olink_|Gene_|EVE_',colnames(df1))) { 
  t1 = df1[,c("D28_IgG",'Info_onset','olink_DDX58_O95786_OID01018',
              colnames(df1)[i])]%>%na.omit()
  colnames(t1)[4]='y'
  LM1 = lm(D28_IgG~poly(Info_onset,2)+y,t1)
  LM1 = summary(LM1)
  
  LM2 = lm(olink_DDX58_O95786_OID01018~poly(Info_onset,2)+y
           ,t1)
  LM2 = summary(LM2)
  t1 = data.frame(gene = colnames(df1)[i],
                  beta = LM1$coefficients['y','Estimate'], 
                  t= LM1$coefficients['y','t value'],
                  p = LM1$coefficients['y','Pr(>|t|)'],
                  p_rig = LM2$coefficients['y','Pr(>|t|)'],
                  t_rig = LM2$coefficients['y','t value'])
  dfResult = rbind(dfResult,t1)
}

write.csv(dfResult,'Result/A15_Ab_M4_results_new_control_only.csv',row.names = F)


dfResultAll = read.csv('Result/A15_Ab_M4_results_new.csv')
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
pdf('Result/A15_severity_correlation_arms.pdf',width = 8,height = 5)  
print(p)
dev.off()

dfResultAll = dfResultAll%>%group_by(type)%>%
  summarise(cor = cor(tAll,tControl,method = 'spearman'))

